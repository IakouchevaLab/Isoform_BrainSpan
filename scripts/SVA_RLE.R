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

# iso_cts <- readRDS("data/iso_cts_filter.rds")
# rle_iso <- list()
# for (i in 1:20) {
#     sv_iso <- readRDS(paste0("data/isoforms/sva_results/", i, ".rds"))
#     lt_cts <- transform_data(
#         voom(calcNormFactors(DGEList(iso_cts), method = "TMM"))$E, mm, sv_iso$sv
#     )
#     rle_iso[[i]] <- sweep(lt_cts, 1, apply(lt_cts, 1, median), "-") %>%
#         reshape2::melt() %>%
#         merge(., metadata, by.x = "Var2", by.y = "Sample") %>%
#         mutate(Period = as.factor(Period)) %>%
#         group_by(Var2) %>%
#         mutate(
#             Quantile75 = quantile(value, 0.75, na.rm = TRUE), 
#             Quantile25 = quantile(value, 0.25, na.rm = TRUE)
#         ) %>%
#         filter(value <= Quantile75 & value >= Quantile25) %>%
#         mutate(
#             Max = max(value),
#             Quantile75 = quantile(value, 0.75, na.rm = TRUE), 
#             Median = median(value),
#             Quantile25 = quantile(value, 0.25, na.rm = TRUE),
#             Min = min(value)
#         ) %>%
#         distinct(
#             Var2, Period, 
#             Max, Quantile75, Median, Quantile25, Min
#         )
# }
# gene_cts <- readRDS("data/gene_cts_filter.rds")
# rle_gene <- list()
# for (i in 1:20) {
#     sv_gene <- readRDS(paste0("data/genes/sva_results/", i, ".rds"))
#     lt_cts <- transform_data(
#         voom(calcNormFactors(DGEList(gene_cts), method = "TMM"))$E, mm, sv_gene$sv
#     )
#     rle_gene[[i]] <- sweep(lt_cts, 1, apply(lt_cts, 1, median), "-") %>%
#         reshape2::melt() %>%
#         merge(., metadata, by.x = "Var2", by.y = "Sample") %>%
#         mutate(Period = as.factor(Period)) %>%
#         group_by(Var2) %>%
#         mutate(
#             Quantile75 = quantile(value, 0.75, na.rm = TRUE), 
#             Quantile25 = quantile(value, 0.25, na.rm = TRUE)
#         ) %>%
#         filter(value <= Quantile75 & value >= Quantile25) %>%
#         mutate(
#             Max = max(value),
#             Quantile75 = quantile(value, 0.75, na.rm = TRUE), 
#             Median = median(value),
#             Quantile25 = quantile(value, 0.25, na.rm = TRUE),
#             Min = min(value)
#         ) %>%
#         distinct(
#             Var2, Period, 
#             Max, Quantile75, Median, Quantile25, Min
#         )
# }
# saveRDS(rle_iso, "data/isoform_SV_RLEs.rds")
# saveRDS(rle_gene, "data/gene_SV_RLEs.rds")


rle_iso <- readRDS("data/isoform_SV_RLEs.rds")
rle_iso_df <- mapply(
    function(rle, sv) {
        mutate(rle, SV = sv)
    }, rle_iso, seq_along(rle_iso),
    SIMPLIFY = FALSE
) %>% 
    bind_rows() %>% 
    mutate(Period = as.numeric(as.character(Period))) %>%
    left_join(
        metadata, by = c("Var2" = "Sample", "Period")
    ) %>%
    ungroup() %>%
    mutate(Sample = factor(Var2, levels = ordered_samples)) %>%
    mutate(
        xstart = as.numeric(Sample) - 0.5,
        xend = as.numeric(Sample) + 0.5
    )
rle_gene <- readRDS("data/gene_SV_RLEs.rds")
rle_gene_df <- mapply(
    function(rle, sv) {
        mutate(rle, SV = sv)
    }, rle_gene, seq_along(rle_gene),
    SIMPLIFY = FALSE
) %>% 
    bind_rows() %>% 
    mutate(Period = as.numeric(as.character(Period))) %>%
    left_join(
        metadata, by = c("Var2" = "Sample", "Period")
    ) %>%
    ungroup() %>%
    mutate(Sample = factor(Var2, levels = ordered_samples)) %>%
    mutate(
        xstart = as.numeric(Sample) - 0.5,
        xend = as.numeric(Sample) + 0.5
    )
isoform_rle_plot <- ggplot() +
    facet_wrap(~ SV)
    

ordered_samples <- metadata %>%
    arrange(Days, Regioncode) %>%
    pull(Sample)
    
gene_rle_plot <- ggplot() +
    facet_wrap(~ SV) +
    geom_segment(
        data = rle_gene_df,
        mapping = aes(
            x = as.numeric(Sample), xend = as.numeric(Sample),
            y = Min, yend = Max
        ),
        size = 0.1
    ) +
    geom_rect(
        data = rle_gene_df,
        mapping = aes(
            xmin = xstart, xmax = xend, ymin = Quantile25, ymax = Quantile75,
            fill = as.factor(Period)
        )
    ) +
    geom_segment(
        data = rle_gene_df,
        mapping = aes(
            x = xstart, y = Quantile75,
            xend = xend, yend = Quantile75
        )
    ) +
    geom_segment(
        data = rle_gene_df,
        mapping = aes(
            x = xstart, y = Quantile25,
            xend = xend, yend = Quantile25
        )
    ) +
    geom_segment(
        data = rle_gene_df,
        mapping = aes(
            x = xstart, y = Median,
            xend = xend, yend = Median
        )
    ) +
    labs(
        x = "Sample", y = "Relative Log Expression"
    ) +
    guides(
        fill = guide_legend(title = "Period")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
        text = element_text(size = 24),
        axis.text = element_blank()
    )

iso_rle_plot <- ggplot() +
    facet_wrap(~ SV) +
    geom_segment(
        data = rle_iso_df,
        mapping = aes(
            x = as.numeric(Sample), xend = as.numeric(Sample),
            y = Min, yend = Max
        ),
        size = 0.1
    ) +
    geom_rect(
        data = rle_iso_df,
        mapping = aes(
            xmin = xstart, xmax = xend, ymin = Quantile25, ymax = Quantile75,
            fill = as.factor(Period)
        )
    ) +
    geom_segment(
        data = rle_iso_df,
        mapping = aes(
            x = xstart, y = Quantile75,
            xend = xend, yend = Quantile75
        )
    ) +
    geom_segment(
        data = rle_iso_df,
        mapping = aes(
            x = xstart, y = Quantile25,
            xend = xend, yend = Quantile25
        )
    ) +
    geom_segment(
        data = rle_iso_df,
        mapping = aes(
            x = xstart, y = Median,
            xend = xend, yend = Median
        )
    ) +
    labs(
        x = "Sample", y = "Relative Log Expression"
    ) +
    guides(
        fill = guide_legend(title = "Period")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
        text = element_text(size = 24),
        axis.text = element_blank()
    )

ggsave(
    filename = "data/figures/IsoformSVARLE.pdf",
    plot = iso_rle_plot,
    device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "data/figures/GeneSVARLE.pdf",
    plot = gene_rle_plot,
    device = "pdf", width = 16, height = 12
)
