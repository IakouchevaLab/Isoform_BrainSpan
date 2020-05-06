library(tidyverse)
library(biomaRt)
library(ggrepel)
library(ggsignif)
library(doParallel)
source("scripts/utility/plot_expression.R")

registerDoParallel(cores = detectCores() - 1)

lof <- c(
    "splice_donor_variant" = "Splice Donor",
    "splice_acceptor_variant" = "Splice Acceptor",
    "frameshift_variant" = "Frameshift",
    "stop_gained" = "Stop Gain",
    "start_lost" = "Start Loss"
)
mutations <- read_tsv("data/SNVs/SatterstromProcessedVEP.txt")
lof_mutations <- mutations %>%
    filter(LoF)

anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]
lof_mutations %>%
    left_join(
        anno %>%
            dplyr::select(
                ensembl_gene_id, ensembl_transcript_id,
                external_gene_id, external_transcript_id
            ),
        by = c("ensembl_gene_id", "ensembl_transcript_id")
    ) %>%
    write_tsv("data/SNVs/LOF_MUTATIONS.tsv")
# # Add extra information to annotation table
# ensembl <- useMart(
#     biomart = "ensembl",
#     dataset = "hsapiens_gene_ensembl",
#     host = "GRCh37.ensembl.org"
# )
# anno_expand <- anno %>%
#     left_join(
#         getBM(
#             attributes = c(
#                 "ensembl_gene_id", "gene_biotype",
#                 "start_position", "end_position"
#             ),
#             filters = "ensembl_gene_id",
#             values = unique(.$ensembl_gene_id),
#             mart = ensembl
#         ),
#         by = "ensembl_gene_id"
#     ) %>%
#     mutate(gene_length = abs(start_position - end_position))
# 
# sv_gn <- 13
# sv_tx <- 16
# 
# tt_genes <- readRDS(
#     paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
# ) %>%
#     left_join(anno_expand, by = "ensembl_gene_id") %>%
#     dplyr::select(-contains("transcript")) %>%
#     distinct() %>%
#     rename(biotype = gene_biotype) %>%
#     mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
#     mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
#     left_join(
#         lof_mutations %>%
#             dplyr::select(
#                 Chromosome, Variant_start, Variant_end, Ref, Alt,
#                 Affected_status, ensembl_gene_id, ensembl_transcript_id
#             ),
#         by = "ensembl_gene_id"
#     ) %>%
#     mutate(Data = "Gene")
# tt_iso <- readRDS(
#     paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
# ) %>%
#     left_join(anno_expand, by = "ensembl_transcript_id") %>%
#     rename(biotype = transcript_biotype) %>%
#     mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
#     mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
#     left_join(
#         lof_mutations %>%
#             dplyr::select(
#                 Chromosome, Variant_start, Variant_end, Ref, Alt,
#                 Affected_status, ensembl_gene_id, ensembl_transcript_id
#             ),
#         by = c("ensembl_transcript_id", "ensembl_gene_id")
#     ) %>%
#     mutate(Data = "Isoform")
# 
# de_genes_list <- tt_genes %>%
#     filter(Significant) %>%
#     pull(selector)
# 
# de_iso_list <- tt_iso %>%
#     filter(Significant) %>%
#     pull(selector)
# 
# tt_combined <- bind_rows(
#     tt_genes %>%
#         mutate(Specific = ! selector %in% de_iso_list),
#     tt_iso %>%
#         mutate(Specific = ! selector %in% de_genes_list)
# )
# 
# write_tsv(tt_combined, "data/limmaResultsAll_LoFTargets.tsv")

tt_combined <- read_tsv(
    "data/limmaResultsAll_LoFTargets.tsv",
    col_types = cols(.default = "c")
)

metadata <- read_tsv("data/metadata.tsv")
# iexpr <- readRDS("data/RegressIsoformCounts.rds")
# gexpr <- readRDS("data/RegressGeneCounts.rds")

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

################################################################################
# Check expression profiles of targets vs nontargets                           #
################################################################################

selectors <- tt_combined %>%
    filter(Data == "Isoform") %>%
    mutate(
        LoF = ensembl_transcript_id %in% (
            lof_mutations %>%
                filter(Affected_status == 2) %>%
                pull(ensembl_transcript_id)
        )
    ) %>%
    filter(as.logical(Specific) & as.logical(Significant) & LoF) %>%
    group_by(ensembl_gene_id) %>%
    mutate(n_iso = n()) %>%
    ungroup() %>%
    filter(n_iso > 1) %>%
    mutate(sel = paste(ensembl_gene_id, Contrast, sep = "_")) %>%
    distinct(ensembl_gene_id, Contrast, sel)

lof_targets <- lof_mutations %>%
    filter(Affected_status == 2) %>%
    pull(ensembl_transcript_id)

gtpm <- readRDS("data/gene_tpm_filter.rds")
itpm <- readRDS("data/iso_tpm_filter.rds")
# gexpr <- readRDS("data/RegressGeneCounts.rds")
# iexpr <- readRDS("data/RegressIsoformCounts.rds")

# selected_expressions <- mclapply(
#     setNames(nm = c(
#         "P02P03", "P03P04", "P04P05", "P05P06", "P06P07", "P07P08", "P08P09",
#         "P09P10", "P10P11", "P11P12", "P12P13", "PrePost"
#     )),
#     function(contrast) {
#         plts <- apply(
#             filter(selectors, Contrast == contrast), 1, function(param) {
#                 ensg <- param[["ensembl_gene_id"]]
#                 ext <- tt_combined %>%
#                     filter(ensembl_gene_id == ensg) %>%
#                     pull(external_gene_id) %>%
#                     unique()
#                 enst <- tt_combined %>%
#                     filter(Data == "Isoform") %>%
#                     filter(Contrast == contrast) %>%
#                     filter(ensembl_gene_id == ensg) %>%
#                     pull(ensembl_transcript_id) %>%
#                     unique()
#                 significant_tx <- tt_combined %>%
#                     filter(Data == "Isoform") %>%
#                     filter(Contrast == contrast) %>%
#                     filter(ensembl_gene_id == ensg) %>%
#                     filter(as.logical(Significant)) %>%
#                     pull(ensembl_transcript_id)
#                 tpm <- rbind(
#                     gexpr[ensg, , drop = FALSE],
#                     iexpr[enst, , drop = FALSE]
#                 ) %>%
#                     expr_to_z() %>%
#                     as.data.frame() %>%
#                     mutate(Feature = rownames(.)) %>%
#                     mutate(
#                         Data = ifelse(
#                             str_detect(Feature, "ENSG"), "Gene", "Isoform"
#                         )
#                     ) %>%
#                     mutate(
#                         Significant = factor(ifelse(
#                             Data == "Gene", "Gene", 
#                             ifelse(
#                                 Feature %in% significant_tx, 
#                                 "DE Isoform", "Not DE Isoform"
#                             )
#                         ), levels = c("Gene", "DE Isoform", "Not DE Isoform")),
#                         LoF = factor(ifelse(
#                             Data == "Gene", "LoF Target", 
#                             ifelse(
#                                 Feature %in% lof_targets, 
#                                 "LoF Target", "Not LoF Target"
#                             )
#                         ), levels = c("LoF Target", "Not LoF Target"))
#                     ) %>%
#                     gather(
#                         key = "Sample", value = "TPM Z-Score", 
#                         starts_with("BrainSpan")
#                     ) %>%
#                     left_join(metadata, by = "Sample") %>%
#                     mutate(GeneSymbol = ext)
#                 p <- ggplot(
#                     data = tpm,
#                     mapping = aes(
#                         x = Days, y = `TPM Z-Score`, group = Feature,
#                         colour = Significant, alpha = LoF
#                     )
#                 ) +
#                     facet_wrap( ~ GeneSymbol) +
#                     geom_line(
#                         stat = "smooth", method = "loess", size = 2
#                     ) +
#                     labs(
#                         title = paste(contrast, ext)
#                     ) +
#                     scale_colour_manual(
#                         values = c(
#                             "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
#                             "DE Isoform" = rgb(84 / 255, 188 / 255, 194 / 255),
#                             "Not DE Isoform" = "grey"
#                         ),
#                         drop = FALSE
#                     ) +
#                     scale_alpha_manual(
#                         values = c(
#                             "Gene" = 1,
#                             "LoF Target" = 1,
#                             "Not LoF Target" = 0.25
#                         ),
#                         drop = FALSE
#                     )
#                 p <- expression_plot_fill(p, metadata) +
#                     theme_bw() +
#                     theme(text = element_text(size = 14))
#                 return(p)
#             }
#         )
#     }
# )
# 
# dir.create("data/figures/ExpressionProfiles")
# for (contrast in names(selected_expressions)) {
#     fn <- paste0(
#         "data/figures/ExpressionProfiles/SpecificDELoFTargetsAndSiblings",
#         contrast, ".pdf"
#     )
#     pdf(fn, width = 16, height = 9)
#     invisible(lapply(selected_expressions[[contrast]], print))
#     dev.off()
# }

################################################################################
# Generate targets and siblings for specific DE                                #
################################################################################

targets_ensg <- tt_combined %>%
    filter(Data == "Isoform") %>%
    filter(
        as.logical(Specific) & as.logical(Significant)
    ) %>%
    mutate(
        LoF = ensembl_transcript_id %in% (
            lof_mutations %>%
                filter(Affected_status == 2) %>%
                pull(ensembl_transcript_id)
        )
    ) %>%
    filter(LoF) %>%
    mutate(selector = paste(ensembl_gene_id, Contrast))
targets_and_sibs <- tt_combined %>%
    filter(Data == "Isoform") %>%
    mutate(selector = paste(ensembl_gene_id, Contrast)) %>%
    filter(selector %in% targets_ensg$selector) %>%
    mutate(
        LoF = ensembl_transcript_id %in% (
            lof_mutations %>%
                filter(Affected_status == 2) %>%
                pull(ensembl_transcript_id)
        )
    ) %>%
    dplyr::select(
        -selector
    ) %>%
    mutate(SatterstromASD = external_gene_id %in% asd_genes)
write_tsv(
    targets_and_sibs, "data/SNVs/SpecificDEIsoformTargetsAndSiblings.tsv"
)

biotype_breakdown <- targets_and_sibs %>%
    distinct(ensembl_transcript_id, biotype, LoF, Contrast) %>%
    count(biotype, LoF, Contrast) %>%
    spread(key = Contrast, value = n, fill = 0) %>%
    mutate(
        Total_nonOverlap = mapply(
            function(bt, lof) {
                targets_and_sibs %>%
                    distinct(ensembl_transcript_id, biotype, LoF) %>%
                    filter(biotype == bt & LoF == lof) %>%
                    nrow()
            },
            biotype, LoF
        )
    ) %>%
    bind_rows(
        .,
        t(data.frame(
            Total_nonOverlap = sapply(
                setNames(nm = colnames(dplyr::select(., starts_with("P")))),
                function(contrast) {
                    targets_and_sibs %>%
                        filter(Contrast == contrast) %>%
                        distinct(ensembl_transcript_id) %>%
                        nrow()
                }
            ),
            Total_LoFTarget = sapply(
                setNames(nm = colnames(dplyr::select(., starts_with("P")))),
                function(contrast) {
                    targets_and_sibs %>%
                        filter(Contrast == contrast) %>%
                        filter(LoF) %>%
                        distinct(ensembl_transcript_id) %>%
                        nrow()
                }
            ),
            Total_NotLoFTarget = sapply(
                setNames(nm = colnames(dplyr::select(., starts_with("P")))),
                function(contrast) {
                    targets_and_sibs %>%
                        filter(Contrast == contrast) %>%
                        filter(!LoF) %>%
                        distinct(ensembl_transcript_id) %>%
                        nrow()
                }
            ), 
            row.names = colnames(dplyr::select(., starts_with("P")))
        )) %>%
            as.data.frame() %>%
            mutate(biotype = rownames(.))
    )

# LoF target more pc vs nonpc?

biotype_DE_lof_fisher <- lapply(
    setNames(nm = colnames(dplyr::select(biotype_breakdown, starts_with("P")))),
    function(contrast) {
        a <- biotype_breakdown %>%
            filter(biotype == "protein_coding" & LoF) %>%
            pull(!!sym(contrast)) %>%
            sum()
        b <- biotype_breakdown %>%
            filter(biotype != "protein_coding" & LoF) %>%
            pull(!!sym(contrast)) %>%
            sum()
        c <- biotype_breakdown %>%
            filter(biotype == "protein_coding" & ! LoF) %>%
            pull(!!sym(contrast)) %>%
            sum()
        d <- biotype_breakdown %>%
            filter(biotype != "protein_coding" & ! LoF) %>%
            pull(!!sym(contrast)) %>%
            sum()
        fisher.test(matrix(c(a, b, c, d), byrow = T, nrow = 2, ncol = 2)) %>%
            broom::tidy() %>%
            mutate(Contrast = contrast)
    }
) %>%
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH"))

# Narrow down to ASD genes with at least 2 isoforms differentially impacted

asd_differentially_targeted <- targets_and_sibs %>%
    filter(
        ensembl_gene_id %in% (
            targets_and_sibs %>%
                filter(SatterstromASD) %>%
                distinct(ensembl_gene_id, LoF) %>%
                count(ensembl_gene_id) %>%
                filter(n > 1) %>%
                pull(ensembl_gene_id)
            )
    ) %>%
    arrange(ensembl_gene_id, ensembl_transcript_id, Contrast)

collapse_contrast <- function(ctr) {
    paste(ctr[!is.na(ctr)], collapse = ", ")
}

asd_differentially_targeted_summarize <- asd_differentially_targeted %>%
    distinct(
        external_gene_id, external_transcript_id, LoF, Contrast, Significant
    ) %>%
    mutate(Contrast = ifelse(Significant, Contrast, NA)) %>%
    group_by(external_gene_id, external_transcript_id, LoF) %>%
    summarise(Contrast = collapse_contrast(Contrast)) %>%
    ungroup() %>%
    mutate(Contrast = ifelse(Contrast == "", NA, Contrast))

asd_diffTar_ctImpact <- asd_differentially_targeted_summarize %>%
    count(external_gene_id, LoF) %>%
    spread(key = LoF, value = n) %>%
    dplyr::select(external_gene_id, `TRUE`, `FALSE`)
colnames(asd_diffTar_ctImpact) <- c(
    "GeneSymbol", "NumImpactedIso", "NumNotImpactedIso"
)

asd_differentially_targeted_pc <- targets_and_sibs %>%
    filter(biotype == "protein_coding") %>%
    filter(
        ensembl_gene_id %in% (
            targets_and_sibs %>%
                filter(SatterstromASD) %>%
                distinct(ensembl_gene_id, LoF) %>%
                count(ensembl_gene_id) %>%
                filter(n > 1) %>%
                pull(ensembl_gene_id)
        )
    ) %>%
    arrange(ensembl_gene_id, ensembl_transcript_id, Contrast)

asd_differentially_targeted_summarize_pc <- asd_differentially_targeted_pc %>%
    distinct(
        external_gene_id, external_transcript_id, LoF, Contrast, Significant
    ) %>%
    mutate(Contrast = ifelse(Significant, Contrast, NA)) %>%
    group_by(external_gene_id, external_transcript_id, LoF) %>%
    summarise(Contrast = collapse_contrast(Contrast)) %>%
    ungroup() %>%
    mutate(Contrast = ifelse(Contrast == "", NA, Contrast))

asd_diffTar_ctImpact_pc <- asd_differentially_targeted_summarize_pc %>%
    count(external_gene_id, LoF) %>%
    spread(key = LoF, value = n) %>%
    dplyr::select(external_gene_id, `TRUE`, `FALSE`)
colnames(asd_diffTar_ctImpact_pc) <- c(
    "GeneSymbol", "NumImpactedIso", "NumNotImpactedIso"
)

# Export to workbook

targets_and_sibs_allOutput <- list(
    SpecDEIsoTargetsAndSibs = targets_and_sibs,
    BiotypeBreakdown = biotype_breakdown,
    BiotypePCFisher = biotype_DE_lof_fisher,
    DiffTargetASDGenes = asd_differentially_targeted,
    DiffTargetASDGenesSummary = asd_differentially_targeted_summarize,
    DiffTargetASDGenesCounts = asd_diffTar_ctImpact,
    DiffTargetASDGenesPC = asd_differentially_targeted_pc,
    DiffTargetASDGenesSummaryPC = asd_differentially_targeted_summarize_pc,
    DiffTargetASDGenesCountsPC = asd_diffTar_ctImpact_pc
)
targets_and_sibs_allOutput[["README"]] <- data.frame(
    Sheet = names(targets_and_sibs_allOutput),
    Description = c(
        "Differentially expressed isoforms, 
specific to isoform-level differential expression, 
targeted by ASD LoF variants, and their sibling isoforms",
        "Breakdown of transcript biotypes",
        "Fisher-exact test for association between ASD LoF targets 
and protein_coding classification",
        "ASD Genes (from SpecDEIsoTargetsAndSibs sheet) that have 
differentially impacted isoforms",
        "Above sheet, grouped by transcripts 
and in which contrast the transcript is differentially expressed",
        "Above, count LoF targets vs non-targets",
        "DiffTargetASDGenes, protein_coding only",
        "DiffTargetASDGenesSummary, protein_coding only",
        "DiffTargetASDGenesCounts, protein_coding only"
    )
)
writexl::write_xlsx(
    targets_and_sibs_allOutput,
    path = "data/SNVs/SpecificDEIsoformTargetsAndSiblings_Extended.xlsx"
)

# Are ASD SNV targets more likely to be downregulated in PrePost

diff_asd_iso_targets_breakdown <- ggplot(
    data = asd_diffTar_ctImpact_pc %>%
        filter(!is.na(NumImpactedIso) & !is.na(NumNotImpactedIso)) %>%
        gather(
            NumImpactedIso, NumNotImpactedIso, key = "Impact", value = "Count"
        ) %>%
        group_by(GeneSymbol) %>%
        mutate(total = sum(Count)),
    mapping = aes(
        x = GeneSymbol, y = Count, fill = Impact
    )
) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    geom_text(
        mapping = aes(
            x = GeneSymbol, y = 0, 
            label = total
        ),
        vjust = 0, nudge_y = 0.01, size = 5
    ) +
    scale_fill_brewer(
        palette = "Set1",
        labels = c(
            "NumImpactedIso" = "ASD LoF Target",
            "NumNotImpactedIso" = "Untargeted by ASD LoF"
        )
    )
saveRDS(
    diff_asd_iso_targets_breakdown,
    "data/isoforms/figures/DiffASDIsoLoFTargetFreq.rds"
)

# diff_asd_iso_targets_fc <- tt_combined %>%
#     filter(Data == "Isoform") %>%
#     filter(
#         external_gene_id %in% (
#             asd_diffTar_ctImpact_pc %>%
#                 filter(!is.na(NumImpactedIso) & !is.na(NumNotImpactedIso)) %>%
#                 pull(GeneSymbol)
#         )
#     ) %>%
#     dplyr::select(
#         logFC, adj.P.Val, ensembl_gene_id, ensembl_transcript_id, Contrast,
#         external_gene_id, external_transcript_id, biotype
#     ) %>%
#     distinct() %>%
#     mutate(
#         LoF = ensembl_transcript_id %in% pull(
#             filter(lof_mutations, Affected_status == 2), ensembl_transcript_id
#         )
#     ) %>%
#     mutate(Direction = sign(as.numeric(logFC))) %>%
#     filter(adj.P.Val <= 0.05) %>%
#     filter(biotype == "protein_coding")
# diff_asd_iso_targets_fc %>% count(Contrast, LoF, Direction) %>% View()

# Heatmap

library(ggdendro)
library(dendextend)
library(cowplot)

heatmap_input <- tt_combined %>%
    filter(Data == "Isoform") %>%
    filter(
        external_gene_id %in% (
            asd_diffTar_ctImpact_pc %>%
                filter(!is.na(NumImpactedIso) & !is.na(NumNotImpactedIso)) %>%
                pull(GeneSymbol)
        )
    ) %>%
    dplyr::select(
        ensembl_gene_id, ensembl_transcript_id,
        external_gene_id, external_transcript_id, biotype
    ) %>%
    distinct() %>%
    mutate(
        LoF = ensembl_transcript_id %in% pull(
            filter(lof_mutations, Affected_status == 2), ensembl_transcript_id
        )
    ) %>%
    filter(biotype == "protein_coding") %>%
    bind_cols(
        as.data.frame(itpm[.$ensembl_transcript_id, ])
    ) %>%
    gather(key = "Sample", value = "TPM", starts_with("BrainSpan"))

dend_iso <- (1 - as.matrix(cor(
    t(itpm[unique(heatmap_input$ensembl_transcript_id), ]), method="pearson"
))) %>%
    as.dist() %>%
    hclust(method = "ward.D") %>%
    as.dendrogram()
dendata_iso <- dendro_data(dend_iso)
segment_iso <- segment(dendata_iso)
colnames(segment_iso) <- c("y", "x", "yend", "xend")
position_iso <- with(
    dendata_iso$labels,
    data.frame(
        y_center = x, 
        ensembl_transcript_id = as.character(label), 
        height = 1
    )
) %>%
    left_join(
        dplyr::select(heatmap_input, ensembl_transcript_id, LoF),
        by = "ensembl_transcript_id"
    ) %>%
    distinct()
dend_plot <- ggplot(
    data = segment_iso %>%
        left_join(
            position_iso, by = c("y" = "y_center")
        ) %>%
        mutate(colour = ifelse(! LoF | is.na(LoF), "black", "red"))
) +
    geom_segment(
        mapping = aes(x = x, y = y, xend = xend, yend = yend, colour = colour)
    ) +
    geom_tile(
        data = position_iso, 
        mapping = aes(y = y_center, x = 1),
        show.legend = FALSE, fill = NA, colour = NA
    ) +
    scale_x_reverse(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(
        values = setNames(nm = c("black", "red")),
        labels = c(
            "red" = "ASD LoF Target"
        )
    ) +
    # geom_text(
    #     data = position_iso,
    #     mapping = aes(
    #         x = x_center, y = -0.1, 
    #         label = gsub("ME", "M", module), colour = Data
    #     ),
    #     angle = -90, show.legend = FALSE, fontface = "bold", hjust = 0,
    #     size = 5
    # ) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black")
    )
dend_plot
dend_lab_plot <- ggplot() +
    geom_tile(
        data = position_iso,
        mapping = aes(y = y_center, x = 1, fill = LoF), 
        colour = NA
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(
        values = c(
            "TRUE" = rgb(210 / 255, 53 / 255, 43 / 255), 
            "FALSE" = rgb(72 / 255, 126 / 255, 179 / 255)),
        labels = c("TRUE" = "ASD LoF Target"),
        breaks = c("TRUE")
    ) +
    theme_bw()
dend_lab_plot
dend_plot_labeled <- cowplot::plot_grid(
    dend_plot +
        theme(
            axis.title = element_blank(), 
            axis.text = element_blank(),
            plot.margin = margin(),
            legend.position = "none"
        ), 
    dend_lab_plot +
        theme(
            text = element_text(size = 30),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "top",
            legend.title = element_blank()
        ),
    ncol = 2, rel_widths = c(1, 0.1), align = "h", axis = "tb"
)
dend_plot_labeled
heat_plot <- ggplot(
    data = heatmap_input %>%
        left_join(
            position_iso, by = "ensembl_transcript_id"
        ) %>%
        left_join(metadata, by = "Sample") %>%
        group_by(
            ensembl_gene_id, ensembl_transcript_id,
            external_gene_id, external_transcript_id,
            Period, LoF.x, y_center, height
        ) %>%
        summarise(MeanTPM = mean(TPM)),
    mapping = aes(
        x = Period, y = y_center, fill = MeanTPM
    )
) +
    geom_tile() +
    scale_fill_gradient(
        # low = rgb(72 / 255, 126 / 255, 179 / 255), 
        # high = rgb(210 / 255, 53 / 255, 43 / 255), 
        low = "white",
        high = "red",
        trans = "log1p"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw()
heat_plot

heat_and_dend <- plot_grid(
    dend_plot +
        theme(
            axis.title = element_blank(), 
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(),
            legend.position = "none"
        ), 
    dend_lab_plot +
        theme(
            text = element_text(size = 30),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "top",
            legend.title = element_blank()
        ),
    heat_plot +
        scale_x_continuous(
            position = "bottom", expand = c(0, 0),
            breaks = seq(2, 13)
        ) +
        theme(
            text = element_text(size = 30),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "top",
            legend.key.width = unit(5, "lines"),
            panel.border = element_rect(colour = "black"),
            rect = element_blank(),
            plot.margin = margin()
        ),
    ncol = 3, rel_widths = c(0.1, 0.05, 1), align = "h", axis = "tb"
)
heat_and_dend
saveRDS(
    heat_and_dend,
    "data/isoforms/figures/DiffASDIsoLoFTargetHeat.rds"
)
ggsave(
    filename = "data/isoforms/figures/DiffASDIsoLoFTargetHeat.pdf",
    heat_and_dend,
    device = "pdf", width = 16, height = 12
)

