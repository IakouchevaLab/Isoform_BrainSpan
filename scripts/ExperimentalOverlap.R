library(tidyverse)
library(biomaRt)

ensembl <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)

experimental_genes <- c("SCN2A", "DYRK1A", "DLG2", "CELF2")
metadata <- read_tsv("data/metadata.tsv")
feature_anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]
anno_expand <- feature_anno %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "start_position", "end_position"),
            filters = "ensembl_gene_id",
            values = unique(.$ensembl_gene_id),
            mart = ensembl
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(gene_length = abs(start_position - end_position))
sv_gn <- 13
sv_tx <- 16
tt_genes <- readRDS(
    paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_gene_id")
tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")

module_assigns <- list(
    gene = read_tsv(
        "data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv"
    ) %>%
        left_join(
            feature_anno %>% 
                dplyr::select(-contains("transcript")) %>%
                distinct(),
            by = c("Feature" = "ensembl_gene_id")
        ),
    isoform = read_tsv(
        "data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv"
    ) %>%
        left_join(
            feature_anno, by = c("Feature" = "ensembl_transcript_id")
        )
)

module_assigns_experimental <- lapply(
    module_assigns,
    function(x) {
        x %>% filter(external_gene_id %in% experimental_genes)
    }
)

mutation_positions <- c(
    "SCN2A" = 166187838,
    "DYRK1A" = 38865466,
    "DLG2" = 83194295,
    "CELF2" = 11356223
)

experimental_isoform_exons <- module_assigns_experimental$isoform %>%
    left_join(
        getBM(
            attributes = c(
                "exon_chrom_start", "exon_chrom_end", "ensembl_transcript_id"
            ),
            filters = c("ensembl_transcript_id"),
            values = .$Feature,
            mart = ensembl
        ),
        by = c("Feature" = "ensembl_transcript_id")
    ) %>%
    mutate(
        dist_exon_start = mapply(
            function(gsym, start) {
                abs(start - mutation_positions[[gsym]])
            },
            external_gene_id, exon_chrom_start
        ),
        dist_exon_end = mapply(
            function(gsym, end) {
                abs(end - mutation_positions[[gsym]])
            },
            external_gene_id, exon_chrom_end
        )
    ) %>%
    mutate(isTargetedByExperimentalSpliceMutation = mapply(
        function(start, end) {
            if(start <= 3 | end <= 3) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        },
        dist_exon_start, dist_exon_end
    )) %>%
    dplyr::select(
        -contains("exon"), -starts_with("kME")
    ) %>%
    distinct()
target_features <- experimental_isoform_exons %>%
    filter(isTargetedByExperimentalSpliceMutation) %>%
    pull(Feature)
experimental_isoform_exons <- experimental_isoform_exons %>%
    mutate(isTarget = Feature %in% target_features) %>%
    filter(isTargetedByExperimentalSpliceMutation == isTarget) %>%
    dplyr::select(-isTarget) %>%
    left_join(
        tt_iso %>%
            dplyr::select(
                ensembl_transcript_id, logFC, AveExpr, t, P.Value, adj.P.Val, B,
                Contrast
            ),
        by = c("Feature" = "ensembl_transcript_id")
    )

write_tsv(
    experimental_isoform_exons, "data/ModuleAssign_ExpSpliceSiteTargets.tsv"
)

gene_expr <- readRDS("data/RegressGeneCounts.rds")
isoform_expr <- readRDS("data/RegressIsoformCounts.rds")
period_breakpoints <- sapply(
    2:length(unique(metadata$Period)),
    function(i) {
        p1 <- sort(unique(metadata$Period))[[i - 1]]
        p2 <- sort(unique(metadata$Period))[[i]]
        max_p1 <- max(pull(filter(metadata, Period == p1), Days))
        min_p2 <- min(pull(filter(metadata, Period == p2), Days))
        mean(c(max_p1, min_p2))
    }
)
all_expr <- bind_rows(
    gene_expr[unique(experimental_isoform_exons$ensembl_gene_id), ] %>%
        reshape2::melt() %>%
        rename(Feature = Var1, Sample = Var2, Expression = value) %>%
        left_join(
            metadata, by = "Sample"
        ) %>%
        mutate(Data = "Gene") %>%
        mutate(isTarget = TRUE) %>%
        left_join(
            dplyr::select(anno_expand, ensembl_gene_id, external_gene_id),
            by = c("Feature" = "ensembl_gene_id")
        ),
    isoform_expr[unique(experimental_isoform_exons$Feature), ] %>%
        reshape2::melt() %>%
        rename(
            Feature = Var1, Sample = Var2, Expression = value
        ) %>%
        left_join(
            metadata, by = "Sample"
        ) %>%
        mutate(Data = "Isoform") %>%
        mutate(
            isTarget = Feature %in% target_features
        ) %>%
        left_join(
            dplyr::select(anno_expand, ensembl_transcript_id, external_gene_id),
            by = c("Feature" = "ensembl_transcript_id")
        )
) %>%
    distinct()

impacted_expressions <- ggplot() +
    facet_grid(external_gene_id ~ .) +
    geom_point(
        data = all_expr %>%
            mutate(
                external_gene_id = factor(
                    as.character(external_gene_id), 
                    levels = c("SCN2A", "DYRK1A", "DLG2", "CELF2")
                )
            ) %>%
            filter(isTarget),
        mapping = aes(
            x = Days, y = Expression, group = Feature, 
            # size = Data, 
            # alpha = isTarget,
            colour = Data
        ),
        alpha = 0.1, size = 0.5
    ) +
    geom_line(
        data = all_expr %>%
            mutate(
                external_gene_id = factor(
                    as.character(external_gene_id), 
                    levels = c("SCN2A", "DYRK1A", "DLG2", "CELF2")
                )
            ) %>%
            filter(isTarget),
        mapping = aes(
            x = Days, y = Expression, group = Feature, 
            # size = Data, 
            # alpha = isTarget,
            colour = Data
        ),
        stat = "smooth", method = "loess", size = 1
    ) +
    geom_line(
        data = all_expr %>%
            mutate(
                external_gene_id = factor(
                    as.character(external_gene_id), 
                    levels = c("SCN2A", "DYRK1A", "DLG2", "CELF2")
                )
            ) %>%
            filter(!isTarget),
        mapping = aes(
            x = Days, y = Expression, group = Feature, 
        ),
        alpha = 0.3, colour = "grey",
        stat = "smooth", method = "loess", size = 1
    ) +
    geom_vline(xintercept = period_breakpoints) +
    labs(
        y = "Normalized Expression"
    ) +
    guides(
        alpha = FALSE,
        colour = guide_legend(title = "")
    ) +
    scale_x_continuous(
        trans = "log1p", breaks = c(100, 500, 5000, 10000)
    ) +
    scale_y_continuous(
        trans = "log1p",
        breaks = c(10, 100, 1000, 10000)
    ) +
    scale_alpha_discrete(range = c(0.3, 1)) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.text = element_text(size = 24),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing = unit(0, "lines")
    )

ggsave(
    filename = "data/figures/ExperimentalSpliceSiteImpactExpressions.pdf",
    plot = impacted_expressions,
    device = "pdf", width = 12, height = 12
)
