library(tidyverse)
library(edgeR)
library(doParallel)
source("scripts/utility/transform_data.R")

registerDoParallel(cores = detectCores() - 1)

metadata <- read_tsv("data/metadata.tsv")
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
rownames(mm) <- metadata$Sample
colnames(mm) <- make.names(colnames(mm))

iso_cts <- readRDS("data/iso_cts_filter.rds")
pca_tables_iso <- foreach(i = 1:20) %dopar% {
    sv_iso <- readRDS(paste0("data/isoforms/sva_results/", i, ".rds"))
    prc <-
    pca <- cbind(
        prcomp(t(transform_data(
            cpm(calcNormFactors(DGEList(iso_cts), method = "TMM")),
            mm,
            sv_iso$sv
        )), scale = TRUE, center = TRUE)$x,
        metadata
    ) %>%
        mutate(
            PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
        ) %>%
        mutate(SV = i) %>%
        dplyr::select(PC1, PC2, SV, Period, PrePost)
}
saveRDS(pca_tables_iso, "data/isoform_SV_PCAs.rds")

gene_cts <- readRDS("data/gene_cts_filter.rds")
pca_tables_gene <- foreach(i = 1:20) %dopar% {
    sv_iso <- readRDS(paste0("data/genes/sva_results/", i, ".rds"))
    prc <-
        pca <- cbind(
            prcomp(t(transform_data(
                cpm(calcNormFactors(DGEList(gene_cts), method = "TMM")),
                mm,
                sv_gene$sv
            )), scale = TRUE, center = TRUE)$x,
            metadata
        ) %>%
        mutate(
            PrePost = ifelse(Period <= 7, "Prenatal", "Postnatal")
        ) %>%
        mutate(SV = i) %>%
        dplyr::select(PC1, PC2, SV, Period, PrePost)
}
saveRDS(pca_tables_gene, "data/gene_SV_PCAs.rds")

pca_tables_iso <- readRDS("data/isoform_SV_PCAs.rds") %>%
    bind_rows() %>%
    mutate(PrePost = factor(PrePost, levels = c("Prenatal", "Postnatal")))
pca_tables_gene <- readRDS("data/gene_SV_PCAs.rds") %>%
    bind_rows() %>%
    mutate(PrePost = factor(PrePost, levels = c("Prenatal", "Postnatal")))

pca_iso <- ggplot(
    data = pca_tables_iso,
    mapping = aes(
        x = PC1, y = PC2, colour = as.factor(Period), shape = PrePost
    )
) +
    facet_wrap(~ SV) +
    geom_point() +
    guides(
        shape = guide_legend(title = ""),
        colour = guide_legend(title = "Period")
    ) +
    theme_minimal() +
    theme(
        text = element_text(size = 20)
    )
pca_gene <- ggplot(
    data = pca_tables_gene,
    mapping = aes(
        x = PC1, y = PC2, colour = as.factor(Period), shape = PrePost
    )
) +
    facet_wrap(~ SV) +
    geom_point() +
    guides(
        shape = guide_legend(title = ""),
        colour = guide_legend(title = "Period")
    ) +
    theme_minimal() +
    theme(
        text = element_text(size = 20)
    )
pca_iso
pca_gene
ggsave(
    filename = "data/SupplementaryFigures/SVA_ISOFORM.pdf",
    plot = pca_iso, device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "data/SupplementaryFigures/SVA_GENE.pdf",
    plot = pca_gene, device = "pdf", width = 16, height = 12
)
