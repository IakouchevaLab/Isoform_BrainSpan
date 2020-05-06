library(tidyverse)
library(edgeR)
library(GGally)
source("scripts/utility/transform_data.R")

# i <- 1
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

metadata <- read_tsv("data/metadata.tsv")
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
rownames(mm) <- metadata$Sample
colnames(mm) <- make.names(colnames(mm))

################################################################################
# ISOFORMS                                                                     #
################################################################################

iso_cts <- readRDS("data/iso_cts_filter.rds")
sv_iso <- readRDS(paste0("data/isoforms/sva_results/", i, ".rds"))
transformed_iso_cts <- transform_data(
    cpm(calcNormFactors(DGEList(iso_cts), method = "TMM")), 
    mm, 
    sv_iso$sv
)
# PCA
prc <- prcomp(t(transformed_iso_cts), scale = TRUE, center = TRUE)
pca <- cbind(prc$x, metadata) %>%
    mutate(
        PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
    )
pca_legend <- GGally::grab_legend(
    ggplot(data = pca, aes(x = PC1, y = PC2, color = as.factor(Period), shape = PrePost)) +
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
    labs(title = paste0("SV", i)) +
    theme_bw()
# RLE
lt_cts <- transform_data(
    voom(calcNormFactors(DGEList(iso_cts), method = "TMM"))$E, mm, sv_iso$sv
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
    ), nrow = 3, ncol = 3, title = paste0("SV", i)
)
dir.create("data/isoforms/figures/sva_checks")
saveRDS(gg_mat, paste0("data/isoforms/figures/sva_checks/", i, ".rds"))
saveRDS(rle, paste0("data/isoforms/figures/sva_checks/rle_", i, ".rds"))
saveRDS(gg_pairs, paste0("data/isoforms/figures/sva_checks/pca_", i, ".rds"))
pdf(
    paste0("data/isoforms/figures/sva_checks/pca_and_rle_", i, ".pdf"), 
    width = 16, height = 12
)
gg_pairs
rle
dev.off()
ggsave(
    plot = gg_mat,
    filename = paste0("data/isoforms/figures/sva_checks/", i, ".pdf"),
    width = 16, height = 9, device = "pdf"
)
ggsave(
    plot = gg_mat,
    filename = paste0("data/isoforms/figures/sva_checks/", i, ".png"),
    width = 16, height = 9, device = "png", units = "in"
)


################################################################################
# GENES                                                                        #
################################################################################

gene_cts <- readRDS("data/gene_cts_filter.rds")
sv_gn <- readRDS(paste0("data/genes/sva_results/", i, ".rds"))
transformed_gene_cts <- transform_data(
    cpm(calcNormFactors(DGEList(gene_cts), method = "TMM")), 
    mm, 
    sv_gn$sv
)
# PCA
prc <- prcomp(t(transformed_gene_cts), scale = TRUE, center = TRUE)
pca <- cbind(prc$x, metadata) %>%
    mutate(
        PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
    )
pca_legend <- GGally::grab_legend(
    ggplot(data = pca, aes(x = PC1, y = PC2, color = as.factor(Period), shape = PrePost)) +
        geom_point(size = 5) +
        guides(
            color = guide_legend(title = "(Starting) Period"), 
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
    labs(title = paste0("SV", i)) +
    theme_bw()
# RLE
lt_cts <- transform_data(
    voom(calcNormFactors(DGEList(gene_cts), method = "TMM"))$E, mm, sv_gn$sv
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
    ), nrow = 3, ncol = 3, title = paste0("SV", i)
)
dir.create("data/genes/figures/sva_checks")
saveRDS(gg_mat, paste0("data/genes/figures/sva_checks/", i, ".rds"))
saveRDS(rle, paste0("data/genes/figures/sva_checks/rle_", i, ".rds"))
saveRDS(gg_pairs, paste0("data/genes/figures/sva_checks/pca_", i, ".rds"))
pdf(
    paste0("data/genes/figures/sva_checks/pca_and_rle_", i, ".pdf"), 
    width = 16, height = 12
)
gg_pairs
rle
dev.off()
ggsave(
    plot = gg_mat,
    filename = paste0("data/genes/figures/sva_checks/", i, ".pdf"),
    width = 16, height = 9, device = "pdf"
)
ggsave(
    plot = gg_mat,
    filename = paste0("data/genes/figures/sva_checks/", i, ".png"),
    width = 16, height = 9, device = "png", units = "in"
)
