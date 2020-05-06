#!/usr/bin/env Rscript

library(tidyverse)
library(edgeR)
library(GGally)
library(WGCNA)
source("scripts/utility/transform_data.R")
enableWGCNAThreads()

dir.create("data/filtering_intermediates")
load("data/source/BrainSpan.RSEM_Quant.isoform.tpm.RData")
iso_tpm <- tpm
load("data/source/BrainSpan.RSEM_Quant.isoform.counts.RData")
iso_cts <- counts
load("data/source/RSEM_Quant.genes.counts.RData")
gene_cts <- counts[, grepl("BrainSpan", colnames(counts))]
load("data/source/RSEM_Quant.genes.tpm.RData")
gene_tpm <- tpm[, grepl("BrainSpan", colnames(counts))]
metadata <- read_tsv("data/source/brainSpan.phenotype.meta.final.tsv") %>%
    mutate(Sample = paste("BrainSpan", Braincode, Regioncode, sep = "_"))
anno <- read.csv("data/source/annotation.transcript.ensg75.txt", row.names = 1)

################################################################################
# TPM Filter --- TPM >= 0.1 in 25% of samples
################################################################################

retain_isoforms <- rownames(iso_tpm)[apply(
    iso_tpm, 1, function(iso) sum(iso >= 0.1) >= (0.25 * ncol(iso_tpm))
)]
retain_genes <- rownames(gene_tpm)[apply(
    gene_tpm, 1, function(gn) sum(gn >= 0.1) >= (0.25 * ncol(gene_tpm))
)]

anno_filter <- filter(
    anno,
    ensembl_gene_id %in% retain_genes & 
        ensembl_transcript_id %in% retain_isoforms
)

iso_tpm_filtered <- iso_tpm[as.character(anno_filter$ensembl_transcript_id), ]
iso_cts_filtered <- iso_cts[as.character(anno_filter$ensembl_transcript_id), ]
gene_tpm_filtered <- gene_tpm[as.character(anno_filter$ensembl_gene_id), ]
gene_cts_filtered <- gene_cts[as.character(anno_filter$ensembl_gene_id), ]

################################################################################
# Outlier analysis by sample connectivity
################################################################################

normadj <- (0.5 + 0.5 * bicor(log2(0.001 + iso_tpm_filtered)))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z_ku_isoform <- (ku-mean(ku))/sqrt(var(ku))
saveRDS(z_ku_isoform, "data/filtering_intermediates/z_ku_isoform.rds")
isoform_outliers <- z_ku_isoform[z_ku_isoform < -2]
normadj <- (0.5 + 0.5 * bicor(log2(0.001 + gene_tpm_filtered)))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z_ku_gene <- (ku-mean(ku))/sqrt(var(ku))
gene_outliers <- z_ku_gene[z_ku_gene < -2]
saveRDS(z_ku_gene, "data/filtering_intermediates/z_ku_gene.rds")
outliers <- unique(c(names(isoform_outliers), names(gene_outliers)))
writeLines(
    c(
        "Gene-level Outliers", 
        names(gene_outliers), 
        "Isoform-level Outliers", 
        names(isoform_outliers)
    ), 
    "data/filtering_intermediates/outliers"
)

################################################################################
# Subset and save
################################################################################

iso_tpm_filtered <- iso_tpm[
    as.character(anno_filter$ensembl_transcript_id), 
    !colnames(iso_tpm_filtered) %in% outliers
]
iso_cts_filtered <- iso_cts[
    as.character(anno_filter$ensembl_transcript_id), 
    !colnames(iso_cts_filtered) %in% outliers
]
gene_tpm_filtered <- gene_tpm[
    unique(as.character(anno_filter$ensembl_gene_id)), 
    !colnames(gene_tpm_filtered) %in% outliers
]
gene_cts_filtered <- gene_cts[
    unique(as.character(anno_filter$ensembl_gene_id)), 
    !colnames(gene_cts_filtered) %in% outliers
]


saveRDS(iso_tpm_filtered, "data/iso_tpm_filter.rds")
saveRDS(iso_cts_filtered, "data/iso_cts_filter.rds")
saveRDS(gene_tpm_filtered, "data/gene_tpm_filter.rds")
saveRDS(gene_cts_filtered, "data/gene_cts_filter.rds")

################################################################################
# Metadata format and filter
################################################################################

metadata <- metadata[match(colnames(iso_tpm_filtered), metadata$Sample), ]
write_tsv(metadata, "data/metadata.tsv")

# Plot PCA and RLE

metadata <- read_tsv("data/metadata.tsv")
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
rownames(mm) <- metadata$Sample
colnames(mm) <- make.names(colnames(mm))
transformed_iso_cts <- transform_data(
    cpm(calcNormFactors(DGEList(
        readRDS("data/iso_cts_filter.rds")
    ), method = "TMM")), 
    mm, 
    NULL
)
transformed_gn_cts <- transform_data(
    cpm(calcNormFactors(DGEList(
        readRDS("data/gene_cts_filter.rds")
    ), method = "TMM")), 
    mm, 
    NULL
)

# Isoform PCA
prc <- prcomp(t(transformed_iso_cts), scale = TRUE, center = TRUE)
pca <- cbind(prc$x, metadata) %>%
    mutate(
        PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
    )
pca_legend <- GGally::grab_legend(
    ggplot(
        data = pca, 
        mapping = aes(
            x = PC1, y = PC2, color = as.factor(Period), shape = PrePost)
        ) +
        geom_point(size = 5) +
        guides(
            color = guide_legend(title = "Period"), 
            shape = guide_legend(title = "")
        ) +
        theme_bw() +
        theme(
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 14),
            legend.box = "vertical",
            legend.position = "bottom"
        )
)
pca_annotation_fn <- function(data, mapping, ...) {
    return(
        ggplot(
            data = data,
            mapping = mapping
        ) +
            annotation_custom(
                grid::textGrob(
                    paste0(
                        rlang::quo_text(mapping$x), 
                        " (",
                        round(
                            summary(prc)$importance[
                                2, rlang::quo_text(mapping$x)
                                ] * 100, 
                            digits = 2
                        ),
                        "%)"
                    )
                ),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
            )
    )
}
gg_pairs <- ggpairs(
    data = pca,
    columns = c("PC1", "PC2", "PC3"),
    aes(shape = PrePost, colour = as.factor(Period)), 
    legend = pca_legend,
    diag = list(continuous = pca_annotation_fn), 
    upper = list(continuous = wrap("points", size = 1.5)),
    lower = list(continuous = "blank"),
    showStrips = FALSE
) +
    labs(title = "SV0") +
    theme_bw()
# Isoform RLE
lt_cts <- transform_data(
    voom(calcNormFactors(DGEList(
        readRDS("data/iso_cts_filter.rds")
    ), method = "TMM"))$E, mm, NULL
)
medians <- apply(lt_cts, 1, median)
cts_rle <- (lt_cts - medians) %>%
    reshape2::melt() %>%
    merge(., metadata, by.x = "Var2", by.y = "Sample")
rle_plot <- ggplot(
    data = cts_rle %>%
        mutate(Period = as.factor(Period)),
    aes(
        x = reorder(Var2, as.numeric(as.character(Period))), y = value
    )
) +
    facet_grid(. ~ Period, scales = "free_x") +
    geom_boxplot(
        aes(fill = Period), outlier.shape = NA,
        size = 0.1, colour = "black", fatten = 5
    ) +
    coord_cartesian(ylim = c(-5, 5)) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none"
    )
gg_mat <- ggmatrix(
    plots = list(
        gg_pairs[1, 1], gg_pairs[1, 2], gg_pairs[1, 3],
        pca_legend, gg_pairs[2, 2], gg_pairs[2, 3],
        rle_plot,
        gg_pairs[3, 2],
        gg_pairs[3, 3]
    ), nrow = 3, ncol = 3, title = "SV0"
)
saveRDS(gg_mat, paste0("data/isoforms/figures/PCA_RLE.rds"))
saveRDS(rle_plot, paste0("data/isoforms/figures/RLE.rds"))
saveRDS(gg_pairs, paste0("data/isoforms/figures/PCA.rds"))
ggsave(
    plot = gg_mat,
    filename = paste0("data/isoforms/figures/PCA_RLE.pdf"),
    width = 16, height = 9, device = "pdf"
)
ggsave(
    plot = gg_mat,
    filename = paste0("data/isoforms/figures/PCA_RLE.png"),
    width = 16, height = 9, device = "png", units = "in"
)

# Gene PCA
prc <- prcomp(t(transformed_gn_cts), scale = TRUE, center = TRUE)
pca <- cbind(prc$x, metadata) %>%
    mutate(
        PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
    )
pca_legend <- GGally::grab_legend(
    ggplot(
        data = pca, 
        mapping = aes(
            x = PC1, y = PC2, color = as.factor(Period), shape = PrePost)
    ) +
        geom_point(size = 5) +
        guides(
            color = guide_legend(title = "Period"), 
            shape = guide_legend(title = "")
        ) +
        theme_bw() +
        theme(
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 14),
            legend.box = "vertical",
            legend.position = "bottom"
        )
)
pca_annotation_fn <- function(data, mapping, ...) {
    return(
        ggplot(
            data = data,
            mapping = mapping
        ) +
            annotation_custom(
                grid::textGrob(
                    paste0(
                        rlang::quo_text(mapping$x), 
                        " (",
                        round(
                            summary(prc)$importance[
                                2, rlang::quo_text(mapping$x)
                                ] * 100, 
                            digits = 2
                        ),
                        "%)"
                    )
                ),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
            )
    )
}
gg_pairs <- ggpairs(
    data = pca,
    columns = c("PC1", "PC2", "PC3"),
    aes(shape = PrePost, colour = as.factor(Period)), 
    legend = pca_legend,
    diag = list(continuous = pca_annotation_fn), 
    upper = list(continuous = wrap("points", size = 1.5)),
    lower = list(continuous = "blank"),
    showStrips = FALSE
) +
    labs(title = "SV0") +
    theme_bw()
# Gene RLE
lt_cts <- transform_data(
    voom(calcNormFactors(DGEList(
        readRDS("data/gene_cts_filter.rds")
    ), method = "TMM"))$E, mm, NULL
)
medians <- apply(lt_cts, 1, median)
cts_rle <- (lt_cts - medians) %>%
    reshape2::melt() %>%
    merge(., metadata, by.x = "Var2", by.y = "Sample")
rle_plot <- ggplot(
    data = cts_rle %>%
        mutate(Period = as.factor(Period)),
    aes(
        x = reorder(Var2, as.numeric(as.character(Period))), y = value
    )
) +
    facet_grid(. ~ Period, scales = "free_x") +
    geom_boxplot(
        aes(fill = Period), outlier.shape = NA,
        size = 0.1, colour = "black", fatten = 5
    ) +
    coord_cartesian(ylim = c(-5, 5)) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none"
    )
gg_mat <- ggmatrix(
    plots = list(
        gg_pairs[1, 1], gg_pairs[1, 2], gg_pairs[1, 3],
        pca_legend, gg_pairs[2, 2], gg_pairs[2, 3],
        rle_plot,
        gg_pairs[3, 2],
        gg_pairs[3, 3]
    ), nrow = 3, ncol = 3, title = "SV0"
)
saveRDS(gg_mat, paste0("data/genes/figures/PCA_RLE.rds"))
saveRDS(rle_plot, paste0("data/genes/figures/RLE.rds"))
saveRDS(gg_pairs, paste0("data/genes/figures/PCA.rds"))
ggsave(
    plot = gg_mat,
    filename = paste0("data/genes/figures/PCA_RLE.pdf"),
    width = 16, height = 9, device = "pdf"
)
ggsave(
    plot = gg_mat,
    filename = paste0("data/genes/figures/PCA_RLE.png"),
    width = 16, height = 9, device = "png", units = "in"
)

# Individual figures

ggsave(
    filename = "data/genes/figures/PCA.pdf",
    plot = readRDS("data/genes/figures/PCA.rds")[1, 2],
    device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "data/genes/figures/RLE.pdf",
    plot = readRDS("data/genes/figures/RLE.rds"),
    device = "pdf", width = 16, height = 12
)
