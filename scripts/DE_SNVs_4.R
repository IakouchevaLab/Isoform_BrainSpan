# Focus on isoforms

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
# Add extra information to annotation table
ensembl <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)
anno_expand <- anno %>%
    left_join(
        getBM(
            attributes = c(
                "ensembl_gene_id", "gene_biotype",
                "start_position", "end_position"
            ),
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
    left_join(anno_expand, by = "ensembl_gene_id") %>%
    dplyr::select(-contains("transcript")) %>%
    distinct() %>%
    rename(biotype = gene_biotype) %>%
    mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
    mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    left_join(
        lof_mutations %>%
            dplyr::select(
                Chromosome, Variant_start, Variant_end, Ref, Alt,
                Affected_status, ensembl_gene_id, ensembl_transcript_id
            ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(Data = "Gene")
tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id") %>%
    rename(biotype = transcript_biotype) %>%
    mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
    mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    left_join(
        lof_mutations %>%
            dplyr::select(
                Chromosome, Variant_start, Variant_end, Ref, Alt,
                Affected_status, ensembl_gene_id, ensembl_transcript_id
            ),
        by = "ensembl_transcript_id"
    ) %>%
    mutate(Data = "Isoform")

de_genes_list <- tt_genes %>%
    filter(Significant) %>%
    pull(selector)

de_iso_list <- tt_iso %>%
    filter(Significant) %>%
    pull(selector)

tt_combined <- bind_rows(
    tt_genes %>%
        mutate(Specific = ! selector %in% de_iso_list),
    tt_iso %>%
        mutate(Specific = ! selector %in% de_genes_list)
)

write_tsv(tt_combined, "data/limmaResultsAll_LoFTargets.tsv")
tt_combined <- read_tsv(
    "data/limmaResultsAll_LoFTargets.tsv",
    col_types = cols(.default = "c")
)


metadata <- read_tsv("data/metadata.tsv")
# iexpr <- readRDS("data/RegressIsoformCounts.rds")
# gexpr <- readRDS("data/RegressGeneCounts.rds")

itpm <- readRDS("data/iso_tpm_filter.rds")
gtpm <- readRDS("data/gene_tpm_filter.rds")

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

################################################################################
# Volcano plots of case vs control targets                                     #
################################################################################

# lof_case_control_volcano <- ggplot(
#     data = tt_combined %>%
#         filter(Data == "Isoform"),
#     mapping = aes(
#         x = logFC, y = -log10(adj.P.Val)
#     )
# ) +
#     facet_wrap(~ Contrast) +
#     geom_point(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(!Significant),
#         alpha = 0.1
#     ) +
#     geom_point(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Significant) %>%
#             filter(is.na(Phenotype_1ctrl_2case)),
#         alpha = 0.1
#     ) +
#     geom_point(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Significant) %>%
#             filter(!is.na(Phenotype_1ctrl_2case)) %>%
#             arrange(Phenotype_1ctrl_2case) %>%
#             distinct(
#                 ensembl_transcript_id, Contrast, logFC, adj.P.Val, Specific,
#                 Phenotype_1ctrl_2case
#             ) %>%
#             group_by(
#                 ensembl_transcript_id, Contrast, logFC, adj.P.Val, Specific
#             ) %>%
#             summarise(
#                 Phenotype_1ctrl_2case = toString(Phenotype_1ctrl_2case)
#             ) %>%
#             ungroup(),
#         mapping = aes(colour = Phenotype_1ctrl_2case)
#     ) +
#     scale_y_continuous(
#         trans = "log1p", breaks = c(0, 50, 100, 200)
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 30)
#     )
# ggsave(
#     filename = "data/figures/DEI_CaseControl_Volcano.pdf",
#     plot = lof_case_control_volcano,
#     device = "pdf", width = 16, height = 12
# )

################################################################################
# Case v Control Effect Sizes                                                  #
################################################################################

# ggplot(
#     data = tt_combined %>%
#         mutate(Affected_status = factor(Affected_status)),
#     mapping = aes(fill = Affected_status)
# ) +
#     facet_grid(. ~ Contrast) +
#     geom_histogram(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Affected_status == 1) %>%
#             filter(as.numeric(logFC) > 0) %>%
#             filter(as.numeric(adj.P.Val) <= 0.05),
#         mapping = aes(
#             x = abs(as.numeric(logFC)), y = ..count..
#         ),
#         alpha = 0.75
#     ) +
#     geom_histogram(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Affected_status == 1) %>%
#             filter(as.numeric(logFC) < 0) %>%
#             filter(as.numeric(adj.P.Val) <= 0.05),
#         mapping = aes(
#             x = abs(as.numeric(logFC)), y = -..count..
#         ),
#         alpha = 0.75
#     ) +
#     geom_histogram(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Affected_status == 2) %>%
#             filter(as.numeric(logFC) > 0) %>%
#             filter(as.numeric(adj.P.Val) <= 0.05),
#         mapping = aes(
#             x = abs(as.numeric(logFC)), y = ..count..
#         ),
#         alpha = 0.75
#     ) +
#     geom_histogram(
#         data = tt_combined %>%
#             filter(Data == "Isoform") %>%
#             filter(Affected_status == 2) %>%
#             filter(as.numeric(logFC) < 0) %>%
#             filter(as.numeric(adj.P.Val) <= 0.05),
#         mapping = aes(
#             x = abs(as.numeric(logFC)), y = -..count..
#         ),
#         alpha = 0.75
#     ) +
#     theme_bw()
