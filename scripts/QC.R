library(tidyverse)

gene_tpm <- readRDS("data/RegressGeneCounts.rds")
samples <- str_remove(colnames(gene_tpm), "BrainSpan_")
brainspan_qc <- read_tsv("Synapse/PEC_DAC_RNAseq_QC/PEC_DAC_RNAseq_QC_matrix_v2.tsv") %>% 
    filter(Project == "BrainSpan") %>%
    separate(
        Sample_name, into = c("Individual", "RegionCode"), remove = FALSE, sep = '_'
    ) %>%
    arrange(Sample_name)
brainspan_fastqc <- read_tsv("Synapse/PEC_DAC_RNAseq_QC/PEC_DAC_RNAseq_FastQC.tsv") %>%
    filter(str_detect(Sample, "^HSB\\d+\\.\\w+")) %>%
    mutate(Sample_name = str_replace(Sample, "\\.", "_")) %>%
    separate(
        Sample_name, into = c("Individual", "RegionCode"), remove = FALSE, sep = '_'
    ) %>%
    arrange(Sample_name)

# FASTQC Deduplication
ggplot(
    data = brainspan_fastqc,
    mapping = aes(
        x = Individual, y = total_deduplicated_percentage
    )
) +
    geom_boxplot() +
    labs(
        x = "Sample", title = "Total (606) samples"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )
ggplot(
    data = brainspan_fastqc %>%
        filter(Sample_name %in% samples),
    mapping = aes(
        x = Individual, y = total_deduplicated_percentage
    )
) +
    geom_boxplot() +
    labs(
        x = "Sample", title = "Filtered (551) samples"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

# STAR uniquely mapped
ggplot(
    data = brainspan_qc,
    mapping = aes(
        x = Individual, y = STAR_uniquely_mapped_percent
    )
) +
    geom_boxplot() +
    labs(
        x = "Sample", title = "Total (606) samples"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0, units = "in"),
        # strip.text.x = element_text(angle = 60, hjust = 0),
        strip.background = element_blank()
    )
ggplot(
    data = brainspan_qc %>%
        filter(Sample_name %in% samples),
    mapping = aes(
        x = Individual, y = STAR_uniquely_mapped_percent
    )
) +
    geom_boxplot() +
    labs(
        x = "Sample", title = "Filtered (551) samples"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0, units = "in"),
        # strip.text.x = element_text(angle = 60, hjust = 0),
        strip.background = element_blank()
    )

# Picard 3'-5' Bias
ggplot(
    data = brainspan_qc %>%
        rename(`3' Bias` = Picard_CollectRnaSeqMetrics_Median_3prime_bias) %>%
        rename(`5' Bias` = Picard_CollectRnaSeqMetrics_Median_5prime_bias) %>%
        pivot_longer(
            cols = c(`3' Bias`, `5' Bias`),
            names_to = "BiasType", values_to = "Bias"
        ),
    mapping = aes(
        x = Individual, y = Bias, fill = BiasType
    )
) +
    geom_boxplot() +
    labs(
        x = "Sample", title = "Total (606) samples"
    ) +
    scale_y_continuous(limits = c(0, 3)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )
ggplot(
    data = brainspan_qc %>%
        rename(`3' Bias` = Picard_CollectRnaSeqMetrics_Median_3prime_bias) %>%
        rename(`5' Bias` = Picard_CollectRnaSeqMetrics_Median_5prime_bias) %>%
        pivot_longer(
            cols = c(`3' Bias`, `5' Bias`),
            names_to = "BiasType", values_to = "Bias"
        ) %>%
        filter(Sample_name %in% samples),
    mapping = aes(
        x = Individual, y = Bias, fill = BiasType
    )
) +
    geom_boxplot() +
    scale_y_continuous(limits = c(0, 3)) +
    labs(
        x = "Sample", title = "Filtered (551) samples"
    ) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

### Other datsets
qc <- read_tsv("Synapse/PEC_DAC_RNAseq_QC/PEC_DAC_RNAseq_QC_matrix_v2.tsv") %>% 
    arrange(Project, Sample_name) %>%
    bind_rows(
        read_tsv("Synapse/PEC_DAC_RNAseq_Freeze2_QC_matrix.tsv") %>%
            arrange(Project, Sample_name)
    )
# fastqc <- read_tsv("Synapse/PEC_DAC_RNAseq_QC/PEC_DAC_RNAseq_FastQC.tsv") %>%
#     rename(Sample_name = Sample) %>%
#     mutate(
#         Project = sapply(Sample_name, function(s) { pull(filter(qc, str_detect(s, Sample_name)), Project)[1]})
#     )
# ggplot(
#     data = fastqc,
#     mapping = aes(
#         x = Project, y = total_deduplicated_percentage
#     )
# ) +
#     geom_boxplot() +
#     labs(
#         x = "Project", title = "Deduplication rate"
#     ) +
#     scale_y_continuous(limits = c(0, 100)) +
#     theme(
#         text = element_text(size = 20),
#         axis.text.x = element_text(angle = 30, hjust = 1)
#     )
pdf("data/figures/qc_genes_unique_map_pc.pdf", width = 16, height = 9)
ggplot(
    data = qc %>%
        filter(
            Project %in% c("BrainSpan", "BrainGVEX", "CMC", "LIBD")
        ),
    mapping = aes(
        x = Project, y = STAR_uniquely_mapped_percent
    )
) +
    geom_boxplot() +
    labs(
        x = "Project", title = "Uniquely mapped reads"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0, units = "in"),
        strip.background = element_blank()
    )
dev.off()
pdf("data/figures/qc_genes_3prime_5prime.pdf", width = 16, height = 9)
ggplot(
    data = qc %>%
        rename(`3' Bias` = Picard_CollectRnaSeqMetrics_Median_3prime_bias) %>%
        rename(`5' Bias` = Picard_CollectRnaSeqMetrics_Median_5prime_bias) %>%
        pivot_longer(
            cols = c(`3' Bias`, `5' Bias`),
            names_to = "BiasType", values_to = "Bias"
        ) %>%
        filter(
            Project %in% c("BrainSpan", "BrainGVEX", "CMC", "LIBD")
        ),
    mapping = aes(
        x = Project, y = Bias, fill = BiasType
    )
) +
    geom_boxplot() +
    labs(
        x = "Project", title = "Read map bias"
    ) +
    scale_y_continuous(limits = c(0, 3)) +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )
dev.off()

gene_body_coverage <- qc %>%
    filter(
        Project %in% c("BrainSpan")
    ) %>%
    as.data.frame()
outlier_records <- subset(gene_body_coverage, gene_body_coverage$Picard_CollectRnaSeqMetrics_Median_3prime_bias %in% boxplot(gene_body_coverage$Picard_CollectRnaSeqMetrics_Median_3prime_bias)$out)
samples <- colnames(gene_tpm) %>% str_split("_") %>% sapply(function(x) paste(x[2], x[3], sep = "_"))
outlier_records$Picard_CollectRnaSeqMetrics_Median_3prime_bias[which(! outlier_records$Sample_name %in% samples)]
