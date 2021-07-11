library(tidyverse)
library(ggsignif)

metadata <- read_tsv("data/metadata.tsv")
iso_tpm <- readRDS("data/iso_tpm_filter.rds")
gn_tpm <- readRDS("data/gene_tpm_filter.rds")
variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
) %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained", 
            "splice_donor_variant", "splice_acceptor_variant"
        )
    )

asd_lof_transcripts <- variants %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Transcript ID`)
asd_lof_genes <- variants %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Gene ID`)

average_period_expression_tx <- lapply(
    setNames(nm = sort(unique(metadata$Period))),
    function(p) {
        rowMeans(iso_tpm[, pull(filter(metadata, Period == p), Sample)])
    }
) %>%
    bind_cols() %>%
    mutate(ensembl_transcript_id = rownames(iso_tpm)) %>%
    pivot_longer(
        cols = as.character(sort(unique(metadata$Period))),
        names_to = "Period", values_to = "AveragePeriodicExpression"
    ) %>%
    mutate(Period = factor(as.character(Period), levels = as.character(sort(unique(metadata$Period))))) %>%
    mutate(Impacted = ifelse(ensembl_transcript_id %in% asd_lof_transcripts, TRUE, FALSE)) %>%
    mutate(FeatureType = "Isoform")
average_period_expression_gn <- lapply(
    setNames(nm = sort(unique(metadata$Period))),
    function(p) {
        rowMeans(gn_tpm[, pull(filter(metadata, Period == p), Sample)])
    }
) %>%
    bind_cols() %>%
    mutate(ensembl_gene_id = rownames(gn_tpm)) %>%
    pivot_longer(
        cols = as.character(sort(unique(metadata$Period))),
        names_to = "Period", values_to = "AveragePeriodicExpression"
    ) %>%
    mutate(Period = factor(as.character(Period), levels = as.character(sort(unique(metadata$Period))))) %>%
    mutate(Impacted = ifelse(ensembl_gene_id %in% asd_lof_genes, TRUE, FALSE)) %>%
    mutate(FeatureType = "Gene")

average_period_expression <- bind_rows(
    average_period_expression_gn, average_period_expression_tx
)

imp_nonimp_plt <- ggplot(
    data = average_period_expression %>%
        filter(FeatureType == "Isoform"),
    mapping = aes(
        x = Impacted, y = AveragePeriodicExpression,
        fill = Impacted
    )
) +
    facet_grid(FeatureType ~ Period) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
    geom_signif(
        comparisons = list(c("TRUE", "FALSE")),
        map_signif_level = TRUE, tip_length = 0, 
        y_position = log2(max(average_period_expression$AveragePeriodicExpression) * 1.1),
        test = 'wilcox.test'
    ) +
    scale_y_continuous(trans = "log2") +
    labs(
        title = "Expression of impacted vs. non-impacted transcripts",
        y = "Log2(AveragePeriodicExpression)"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "cm")
    )
imp_nonimp_plt
ggsave("data/figures/impact_nonimpact_expressions_isoforms.pdf", width = 16, height = 9)

for (p in unique(as.character(pull(average_period_expression, Period)))) {
    wilcox.test(
        x = average_period_expression %>%
            filter(FeatureType == "Isoform") %>%
            filter(Period == p) %>%
            filter(Impacted == TRUE) %>%
            pull(AveragePeriodicExpression),
        y = average_period_expression %>%
            filter(FeatureType == "Isoform") %>%
            filter(Period == p) %>%
            filter(Impacted == FALSE) %>%
            pull(AveragePeriodicExpression)
    ) %>%
        print()
}
p.adjust(rep(2.2e-16, times = 12), method = "fdr")
