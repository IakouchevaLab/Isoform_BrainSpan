library(WGCNA)
library(tidyverse)
source("scripts/utility/transform_data.R")
allowWGCNAThreads()
enableWGCNAThreads()

NETWORK_DATA <- as.character(commandArgs(trailingOnly = TRUE)[1])
i <- Sys.getenv("SLURM_ARRAY_TASK_ID")

metadata <- read_tsv("data/metadata.tsv")
model_matrix <- model.matrix(
    ~ Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
rownames(model_matrix) <- metadata$Sample

################################################################################
# Data regression                                                              #
################################################################################

gene_cts <- readRDS("data/gene_cts_filter.rds")
iso_cts <- readRDS("data/iso_cts_filter.rds")
gn_sv <- readRDS("data/genes/sva_results/13.rds")
iso_sv <- readRDS("data/genes/sva_results/16.rds")

gene_transform <- transform_data(gene_cts, model_matrix, gn_sv$sv)
iso_transform <- transform_data(iso_cts, model_matrix, iso_sv$sv)

saveRDS(gene_transform, "data/RegressGeneCounts.rds")
saveRDS(iso_transform, "data/RegressIsoformCounts.rds")

gene_transform <- readRDS("data/RegressGeneCounts.rds")
iso_transform <- readRDS("data/RegressIsoformCounts.rds")

################################################################################
# Soft thresholding                                                            #
################################################################################

# powers <- c(c(1:15), seq(from=16, to=30, by=2))
# sft_gene <- pickSoftThreshold(
#     t(gene_transform), powerVector = powers, verbose = 5
# )
# saveRDS(sft_gene, "data/genes/SoftThresholding.rds")
# sft_iso <- pickSoftThreshold(
#     t(iso_transform), powerVector = powers, verbose = 5
# )
# saveRDS(sft_iso, "data/isoforms/SoftThresholding.rds")

sft_gene <- readRDS("data/genes/SoftThresholding.rds")
sft_iso <- readRDS("data/isoforms/SoftThresholding.rds")
# 
# # Plotting
# # Gene
# scale_free_gene <- ggplot(
#     data = sft_gene$fitIndices,
#     mapping = aes(
#         x = Power, y = -sign(slope) * SFT.R.sq, label = Power
#     )
# ) +
#     geom_text(size = 5, colour = "red") +
#     geom_hline(yintercept = 0.8, colour = "red") +
#     labs(
#         x = "Soft Threshold Power", 
#         y = expression("Scale Free Topology Model Fit, signed R"^"2")
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 24)
#     )
# mean_conn_gene <- ggplot(
#     data = sft_gene$fitIndices,
#     mapping = aes(
#         x = Power, y = mean.k., label = Power
#     )
# ) +
#     geom_text(size = 5, colour = "red") +
#     labs(
#         x = "Soft Threshold Power", 
#         y = "Mean Connectivity"
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 24)
#     )
# sft_gene_plot <- cowplot::plot_grid(
#     scale_free_gene, mean_conn_gene, ncol = 2, nrow = 1, align = "hv"
# )
# ggsave(
#     filename = "data/genes/figures/SoftThresholding.pdf",
#     plot = sft_gene_plot,
#     device = "pdf", width = 16, height = 12
# )
# # Isoform
# scale_free_iso <- ggplot(
#     data = sft_iso$fitIndices,
#     mapping = aes(
#         x = Power, y = -sign(slope) * SFT.R.sq, label = Power
#     )
# ) +
#     geom_text(size = 5, colour = "red") +
#     geom_hline(yintercept = 0.8, colour = "red") +
#     labs(
#         x = "Soft Threshold Power", 
#         y = expression("Scale Free Topology Model Fit, signed R"^"2")
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 24)
#     )
# mean_conn_iso <- ggplot(
#     data = sft_iso$fitIndices,
#     mapping = aes(
#         x = Power, y = mean.k., label = Power
#     )
# ) +
#     geom_text(size = 5, colour = "red") +
#     labs(
#         x = "Soft Threshold Power", 
#         y = "Mean Connectivity"
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 24)
#     )
# sft_iso_plot <- cowplot::plot_grid(
#     scale_free_iso, mean_conn_iso, ncol = 2, nrow = 1, align = "hv"
# )
# ggsave(
#     filename = "data/isoforms/figures/SoftThresholding.pdf",
#     plot = sft_iso_plot,
#     device = "pdf", width = 16, height = 12
# )

################################################################################
# TOM construction                                                             #
################################################################################

# geneTOM <- 1 - TOMsimilarity(
#     adjacency(t(gene_transform), power = sft_gene$powerEstimate)
# )
# saveRDS(geneTOM, "data/genes/dissTOM.rds")
# isoformTOM <- 1 - TOMsimilarity(
#     adjacency(t(iso_transform), power = sft_gene$powerEstimate)
# )
# saveRDS(isoformTOM, "data/isoforms/dissTOM.rds")

################################################################################
# Network construction                                                         #
################################################################################

permute_params <- expand.grid(
    DS = c(1, 2, 3, 4),
    MM = seq(10, 50, by = 10)
) 
if (!is.na(as.numeric(i))) {
    permute_params[i, ] %>%
    apply(
        ., 1, function(param) {
            ds <- as.numeric(param[["DS"]])
            mm <- as.numeric(param[["MM"]])
            if(NETWORK_DATA == "gene") {
                geneNetwork <- blockwiseModules(
                    t(gene_transform),
                    maxBlockSize = 45000,
                    power = sft_gene$powerEstimate,
                    TOMType = "signed",
                    minModuleSize = mm,
                    reassignThreshold = 0,
                    mergeCutHeight = 0.25,
                    deepSplit =ds,
                    numericLabels = TRUE,
                    loadTOM = TRUE,
                    saveTOMs = TRUE,
                    saveTOMFileBase = "data/genes/geneTOM",
                    verbose = 5
                )
                saveRDS(
                    geneNetwork,
                    paste0(
                        "data/genes/Networks/Network_DS",
                        ds, "_MM", mm, ".rds"
                    )
                )
            } else if (NETWORK_DATA == "isoform") {
                isoformNetwork <- blockwiseModules(
                    t(iso_transform),
                    maxBlockSize = 45000,
                    power = sft_iso$powerEstimate,
                    TOMType = "signed",
                    minModuleSize = mm,
                    reassignThreshold = 0,
                    mergeCutHeight = 0.25,
                    deepSplit = ds,
                    numericLabels = TRUE,
                    loadTOM = TRUE,
                    saveTOMs = TRUE,
                    saveTOMFileBase = "data/isoforms/isoformTOM",
                    verbose = 5
                )
                saveRDS(
                    isoformNetwork,
                    paste0(
                        "data/isoforms/Networks/Network_DS",
                        ds, "_MM", mm, ".rds"
                    )
                )
            } else {
                geneNetwork <- blockwiseModules(
                    t(gene_transform),
                    maxBlockSize = 45000,
                    power = sft_gene$powerEstimate,
                    TOMType = "signed",
                    minModuleSize = mm,
                    reassignThreshold = 0,
                    mergeCutHeight = 0.25,
                    deepSplit =ds,
                    numericLabels = TRUE,
                    loadTOM = TRUE,
                    saveTOMs = TRUE,
                    saveTOMFileBase = "data/genes/geneTOM",
                    verbose = 5
                )
                saveRDS(
                    geneNetwork,
                    paste0(
                        "data/genes/Networks/Network_DS",
                        ds, "_MM", mm, ".rds"
                    )
                )
                isoformNetwork <- blockwiseModules(
                    t(iso_transform),
                    maxBlockSize = 45000,
                    power = sft_iso$powerEstimate,
                    TOMType = "signed",
                    minModuleSize = mm,
                    reassignThreshold = 0,
                    mergeCutHeight = 0.25,
                    deepSplit = ds,
                    numericLabels = TRUE,
                    loadTOM = TRUE,
                    saveTOMs = TRUE,
                    saveTOMFileBase = "data/isoforms/isoformTOM",
                    verbose = 5
                )
                saveRDS(
                    isoformNetwork,
                    paste0(
                        "data/isoforms/Networks/Network_DS",
                        ds, "_MM", mm, ".rds"
                    )
                )
            }
        }
    )
}