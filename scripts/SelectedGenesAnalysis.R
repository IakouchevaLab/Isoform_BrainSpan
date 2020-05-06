library(tidyverse)
library(cowplot)
library(ggrepel)
library(biomaRt)
source("scripts/utility/plot_expression.R")

################################################################################
# GLOBAL                                                                       #
################################################################################

genes <- c("BTRC", "DYRK1A", "SCN2A", "DLG2", "CELF2")
SNV_locations <- c(
    "BTRC" = 103221816,
    "DYRK1A" = 38865466,
    "SCN2A" = 166187838,
    "DLG2" = 83194295,
    "CELF2" = 11356223
)

################################################################################
# DATA                                                                         #
################################################################################

annotations <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]
tt_combined <- read_tsv(
    "data/limmaResultsAll_LoFTargets.tsv",
    col_types = cols(.default = "c")
)
gene_modules <- lapply(
    setNames(
        nm = readxl::excel_sheets(
            "data/SupplementaryTables/Supplementary Table 7.xlsx"
        )[-1]
    ),
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 7.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()
isoform_modules <- lapply(
    setNames(
        nm = readxl::excel_sheets(
            "data/SupplementaryTables/Supplementary Table 8.xlsx"
        )[-1]
    ),
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 8.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()
variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 6.xlsx",
    sheet = 2
)
control_targets <- c(
    filter(variants, `Affected status` == 1)$`Ensembl Gene ID`,
    filter(variants, `Affected status` == 1)$`Ensembl Transcript ID`
)
asd_targets <- c(
    filter(variants, `Affected status` == 2)$`Ensembl Gene ID`,
    filter(variants, `Affected status` == 2)$`Ensembl Transcript ID`
)
asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]
asd_genes_ensembl <- as.character(annotations$ensembl_gene_id[match(
    asd_genes, annotations$external_gene_id
)])
feature_metadata <- bind_rows(
    gene_modules %>%
        mutate(name = `Ensembl Gene ID`) %>%
        mutate(Data = "Gene"),
    isoform_modules %>%
        mutate(name = `Ensembl Transcript ID`) %>%
        mutate(Data = "Isoform")
) %>%
    dplyr::select(name, everything()) %>%
    mutate(ASD = `Ensembl Gene ID` %in% asd_genes_ensembl) %>%
    mutate(
        Target = case_when(
            name %in% control_targets & name %in% asd_targets ~ "Both",
            ! name %in% control_targets & name %in% asd_targets ~ "ASD",
            name %in% control_targets & ! name %in% asd_targets ~ "Control",
            ! name %in% control_targets & ! name %in% asd_targets ~ "Neither"
        )
    )
metadata <- read_tsv("data/metadata.tsv")

iexpr <- readRDS("data/RegressIsoformCounts.rds")
gexpr <- readRDS("data/RegressGeneCounts.rds")

################################################################################
# ANALYSIS                                                                     #
################################################################################

# Expression profiles
impacted_isoforms <- c()
expression_plots <- vector(mode = "list", length = length(genes)) %>%
    setNames(nm = genes)
for (gene in genes) {
    features <- feature_metadata %>%
        filter(Data == "Isoform" & `Gene Symbol` == gene) %>%
        pull(`Ensembl Transcript ID`)
    # Get exons and check which isoforms are impacted by given SNV
    impacted_exons <- getBM(
        attributes = c(
            "ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"
        ),
        filters = c("ensembl_transcript_id"),
        values = features, 
        mart = useMart(
            biomart = "ensembl", 
            dataset = "hsapiens_gene_ensembl",
            host = "GRCh37.ensembl.org"
        )
    ) %>%
        mutate(
            SNV_dist = mapply(
                function(start, end) {
                    min(
                        abs(SNV_locations[[gene]] - start),
                        abs(SNV_locations[[gene]] - end)
                    )
                }, exon_chrom_start, exon_chrom_end
            )
        ) %>%
        filter(SNV_dist <= 2) %>%
        pull(ensembl_transcript_id)
    impacted_isoforms <- c(impacted_isoforms, impacted_exons)
    
    feature_expressions <- iexpr[features, ] %>%
        as.matrix() %>%
        expr_to_z() %>%
        as.data.frame() %>%
        mutate(Feature = rownames(.)) %>%
        gather(
            key = "Sample", value = "Expression Z-Score", 
            starts_with("BrainSpan")
        ) %>%
        left_join(
            filter(feature_metadata, `Gene Symbol` == gene & Data == "Isoform"),
            by = c("Feature" = "name")
        ) %>%
        left_join(metadata, by = "Sample") %>%
        mutate(isImpacted = Feature %in% impacted_exons)
    final_loess_points <- feature_expressions %>%
        distinct(
            Feature, `Transcript Symbol`, `Module Colour`, `Module Label`
        ) %>%
        mutate(
            FinalLoessPrediction = sapply(
                Feature, function(feat) {
                    this_expr <- feature_expressions %>%
                        filter(Feature == feat) %>%
                        mutate(Feature = as.character(Feature))
                    loess_predict <- loess(
                        `Expression Z-Score` ~ Days, data = this_expr
                    ) %>%
                        predict(
                            newdata = data.frame(Days = max(this_expr$Days))
                        )
                    return(loess_predict)
                }
            )
        ) %>%
        mutate(Label = paste0(`Transcript Symbol`, "_iM", `Module Label`))
    
    expression_plots[[gene]] <- expression_plot_fill(ggplot(), metadata) +
        facet_grid(. ~ `Gene Symbol`) +
        geom_line(
            data = feature_expressions,
            mapping = aes(
                x = Days, y = `Expression Z-Score`,
                group = Feature, colour = isImpacted
            ),
            stat = "smooth", method = "loess", size = 0.5
        ) +
        geom_text_repel(
            data = final_loess_points,
            mapping = aes(
                x = max(feature_expressions$Days), y = FinalLoessPrediction,
                label = Label
            ),
            size = 2
        ) +
        labs(
            x = "Days",
            y = "Expression Z-Score"
        ) +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
        theme_bw() +
        theme(
            text = element_text(size = 10),
        )
}
# pdf("data/figures/ExperimentalGenes_ExpressionProfiles.pdf", width = 12, height = 9)
# invisible(lapply(expression_plots, print))
# dev.off()

four_plot <- plot_grid(
    plotlist = lapply(
        list(
            expression_plots$SCN2A, expression_plots$DYRK1A, 
            expression_plots$DLG2, expression_plots$CELF2
        ),
        function(plt) plt + 
            theme(
                text = element_text(size = 7),
                legend.position = "none"
            )
    ),
    ncol = 2, nrow = 2, align = "hv", axis = "lrtb"
)

ggsave(
    filename = "data/figures/ExperimentalGenes_ExpressionProfiles_4genes.pdf",
    plot = four_plot, device = "pdf", width = 12, height = 6, units = "cm"
) 


################################################################################
# PPI
################################################################################

gene_ppi <- read.table(
    "data/source/CuratedLists/allVidal.lit17.mvp.hprdInvivo.hintBinary.txt",
    header = TRUE, stringsAsFactors = FALSE
)[, 1:2] %>%
    setNames(nm = c("ENSG1", "ENSG2")) %>% 
    as_tibble() %>%
    left_join(
        annotations[, c(
            "ensembl_gene_id", "external_gene_id", 
            "ensembl_transcript_id", "external_transcript_id"
        )],
        by = c("ENSG1" = "ensembl_gene_id")
    ) %>%
    rename(
        SYMBOL1 = external_gene_id,
        ENST1 = ensembl_transcript_id,
        TXSYMBOL1 = external_transcript_id
    ) %>% 
    left_join(
        annotations[, c(
            "ensembl_gene_id", "external_gene_id", 
            "ensembl_transcript_id", "external_transcript_id"
        )],
        by = c("ENSG2" = "ensembl_gene_id")
    ) %>%
    rename(
        SYMBOL2 = external_gene_id,
        ENST2 = ensembl_transcript_id,
        TXSYMBOL2 = external_transcript_id
    ) %>%
    filter(SYMBOL1 %in% genes | SYMBOL2 %in% genes) %>%
    filter(ENST1 %in% rownames(iexpr) & ENST2 %in% rownames(iexpr)) %>%
    mutate(
        PCC = mapply(
            function(x, y) {
                cor(iexpr[x, ], iexpr[y, ])
            }, 
            ENST1, ENST2
        )
    ) %>%
    left_join(
        isoform_modules[, c("Ensembl Transcript ID", "Module Label")],
        by = c("ENST1" = "Ensembl Transcript ID")
    ) %>%
    rename(
        ML1 = `Module Label`
    ) %>%
    left_join(
        isoform_modules[, c("Ensembl Transcript ID", "Module Label")],
        by = c("ENST2" = "Ensembl Transcript ID")
    ) %>%
    rename(
        ML2 = `Module Label`
    ) %>%
    filter(PCC != 1) %>%
    filter(ML1 == ML2) %>%
    filter(ML1 != 0 & ML2 != 0) %>%
    mutate(IMPACT1 = ENST1 %in% impacted_isoforms) %>%
    mutate(IMPACT2 = ENST2 %in% impacted_isoforms) %>%
    filter(ENSG1 != ENSG2)
gene_ppi

lapply(
    setNames(nm = genes),
    function(g) {
        print(nrow(filter(gene_ppi, SYMBOL1 == g | SYMBOL2 == g)))
        print(summary(filter(gene_ppi, SYMBOL1 == g | SYMBOL2 == g)$PCC))
    }
)

edges <- gene_ppi %>% {
    lapply(
        genes,
        function(g) {
            e <- filter(., SYMBOL1 == g | SYMBOL2 == g)
            q <- quantile(e$PCC, probs = 0.8)
            filter(e, PCC >= q)
        }
    ) %>%
        bind_rows()
}
vertices <- gene_ppi %>% {
    df1 <- dplyr::select(., ends_with("1"), -PCC) %>%
        setNames(nm = gsub("1", "", colnames(.)))
    df2 <- dplyr::select(., ends_with("2"), -PCC) %>%
        setNames(nm = gsub("2", "", colnames(.)))
    bind_rows(df1, df2)
} %>%
    distinct() %>%
    mutate(Seed = SYMBOL %in% genes)

net_list = lapply(
    setNames(nm = genes),
    function(g) { filter(edges, SYMBOL1 == g | SYMBOL2 == g)}
)
net_list[["vertices"]] = vertices
writexl::write_xlsx(
    net_list,
    "data/ppi/ExperimentalGenes_Networks.xlsx"
)
