library(tidyverse)

pval_csvs <- list.files(
    path = "data/empirical_correlations_psychencode_pvalues/", 
    pattern = "*.csv",
    full.names = TRUE
) %>%
    setNames(nm = basename(.)) %>%
    lapply(function(f) read_csv(f) %>% mutate(filename = basename(f))) %>%
    lapply(function(x) mutate(x, filename = str_replace(str_replace(filename, "empirical_correlations_", ""), "_iso_p.value.csv", ""))) %>%
    bind_rows() %>%
    group_by(filename) %>%
    mutate(fdr = p.adjust(x, method = "fdr"))

pval_plt <- ggplot(
    data = pval_csvs,
    mapping = aes(x = "", y = x, group = filename)
) +
    facet_wrap(~filename) +
    geom_boxplot() +
    labs(y = "Pearson Correlation P-Value") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )
adj_pval_plt <- ggplot(
    data = pval_csvs,
    mapping = aes(x = "", y = fdr, group = filename)
) +
    facet_wrap(~filename) +
    geom_boxplot() +
    labs(y = "FDR-adjusted Pearson Correlation P-Value") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )
pdf("data/figures/rebuttal_2/PsychencodeIsoformDatasetCorrelationPValues.pdf", paper = "A4")
pval_plt
adj_pval_plt
dev.off()
