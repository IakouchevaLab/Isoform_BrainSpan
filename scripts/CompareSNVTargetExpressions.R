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
        by = c("ensembl_transcript_id", "ensembl_gene_id")
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
# gtpm <- readRDS("data/gene_tpm_filter.rds")

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

################################################################################
# Relative expression of isoforms targeted by case mutations rel to control    #
# WILCOX TEST                                                                  #
################################################################################

case_features <- lof_mutations %>%
    filter(Affected_status == 2) %>%
    pull(ensembl_transcript_id) %>%
    intersect(rownames(itpm))
control_features <- lof_mutations %>%
    filter(Affected_status == 1) %>%
    pull(ensembl_transcript_id) %>%
    intersect(rownames(itpm))
notCase_features <- rownames(itpm)[! rownames(itpm) %in% case_features]
all_features <- list(
    Case = setdiff(case_features, c(control_features, notCase_features)),
    Control = setdiff(control_features, case_features),
    NonCase = notCase_features
)

compare_target_tpm <- lapply(
    setNames(nm = c(sort(unique(metadata$Period)), "Prenatal", "Postnatal")),
    function(period) {
        this_samples <- metadata %>%
            filter(
                if(period %in% c("Prenatal", "Postnatal")) {
                    if (period == "Prenatal") { Prenatal }
                    else { ! Prenatal }
                } else {
                    Period == as.numeric(period)
                }
            ) %>%
            pull(Sample)
        case_rowmeans <- data.frame(
            Feature = all_features[["Case"]],
            Meantpm = rowMeans(itpm[all_features[["Case"]], this_samples]),
            Data = "Case"
        )
        control_rowmeans <- data.frame(
            Feature = all_features[["Control"]],
            Meantpm = rowMeans(itpm[all_features[["Control"]], this_samples]),
            Data = "Control"
        )
        nonCase_rowmeans <- data.frame(
            Feature = all_features[["NonCase"]],
            Meantpm = rowMeans(itpm[all_features[["NonCase"]], this_samples]),
            Data = "NotCase"
        )
        case_control <- bind_rows(
            data.frame(
                Source = "Case",
                Against = "Control",
                ValueType = "Case",
                Mean = mean(case_rowmeans$Meantpm),
                SD = sd(case_rowmeans$Meantpm),
                Count = nrow(case_rowmeans)
            ),
            data.frame(
                Source = "Case",
                Against = "Control",
                ValueType = "Control",
                Mean = mean(control_rowmeans$Meantpm),
                SD = sd(control_rowmeans$Meantpm),
                Count = nrow(control_rowmeans)
            )
        )
        case_noncase <- bind_rows(
            data.frame(
                Source = "Case",
                Against = "NonCase",
                ValueType = "Case",
                Mean = mean(case_rowmeans$Meantpm),
                SD = sd(case_rowmeans$Meantpm),
                Count = nrow(case_rowmeans)
            ),
            data.frame(
                Source = "Case",
                Against = "NonCase",
                ValueType = "NonCase",
                Mean = mean(nonCase_rowmeans$Meantpm),
                SD = sd(nonCase_rowmeans$Meantpm),
                Count = nrow(nonCase_rowmeans)
            )
        )
        bind_rows(
            case_control,
            case_noncase
        ) %>%
            mutate(Period = period)
    }
) %>%
    bind_rows() %>%
    mutate(
        Period = factor(
            Period, levels = c(
                seq(2, 13), "Prenatal", "Postnatal"
            )
        )
    )

wilcox_test <- expand.grid(
    Source = "Case",
    Against = c("Control", "NonCase"),
    Period = unique(compare_target_tpm$Period)
) %>%
    bind_cols(
        .,
        apply(
            ., 1, function(param) {
                src <- param[["Source"]]
                tar <- param[["Against"]]
                period <- param[["Period"]]
                samples <- metadata %>%
                    filter(
                        if (period %in% c("Prenatal", "Postnatal")) {
                            if (period == "Prenatal") { Prenatal } 
                            else { ! Prenatal }
                        } else {
                            Period == as.numeric(period)
                        }
                    ) %>%
                    pull(Sample)
                src_features <- all_features[[src]]
                tar_features <- all_features[[tar]]
                src_means <- rowMeans(
                    itpm[intersect(rownames(itpm), src_features), samples]
                )
                tar_means <- rowMeans(
                    itpm[intersect(rownames(itpm), tar_features), samples]
                )
                wilcox.test(
                    x = src_means, y = tar_means, paired = FALSE,
                    alternative = "greater"
                ) %>%
                    broom::tidy()
            }
        ) %>%
            bind_rows()
    ) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    mutate(
        Star = sapply(
            adj.P.Val, 
            function(p) {
                if (p <= 0.05) {
                    if (p <= 0.01) {
                        if (p <= 0.001) {
                            return("*\n*\n*")
                        }
                        return("*\n*")
                    }
                    return("*")
                }
                return("")
            }
        )
    ) %>%
    mutate(
        Period = factor(
            Period, levels = c(
                seq(2, 13), "Prenatal", "Postnatal"
            )
        )
    ) %>%
    mutate(
        x = as.numeric(Period) - 0.33, 
        xend = as.numeric(Period) + 0.33
    ) %>%
    mutate(
        y_horiz = sapply(Period, function(per) {
            max(filter(compare_target_tpm, Period == per)$Mean) * 1.05
        })
    )

plot_data <- compare_target_tpm %>%
    left_join(
        wilcox_test %>%
            dplyr::select(Source, Against, Period, x, xend, y_horiz, Star),
        by = c("Source", "Against", "Period")
    ) %>%
    mutate(
        x_ver = mapply(
            function(period, vt) {
                if (vt == "Case") {
                    return(as.numeric(period) - 0.33)
                } else if (vt == "Control") {
                    return(as.numeric(period))
                } else {
                    return(as.numeric(period) + 0.33)
                }
            },
            Period, ValueType
        )
    ) %>%
    mutate(
        x_star = mapply(
            function(tar, per) {
                if (tar == "Control") {
                    return(as.numeric(per) - (0.33 / 2))
                } else {
                    return(as.numeric(per) + (0.33 / 2))
                }
            },
            Against, Period
        )
    )

tpmWilcox_compare_plot <- ggplot() +
    geom_bar(
        data = plot_data,
        mapping = aes(
            x = Period, y = Mean, fill = ValueType
        ),
        stat = "identity", position = "dodge", colour = "black"
    ) +
    geom_segment(
        data = plot_data,
        mapping = aes(
            x = x, xend = xend, y = y_horiz, yend = y_horiz
        )
    ) +
    geom_segment(
        data = plot_data,
        mapping = aes(
            x = x_ver, xend = x_ver, y = y_horiz, yend = Mean * 1.01
        )
    ) +
    geom_text(
        data = plot_data,
        mapping = aes(
            x = x_star, y = y_horiz * 1.01, label = Star
        ),
        fontface = "bold", size = 8,
        lineheight = 0.3, vjust = 0.25
    ) +
    labs(
        x = "Period",
        y = "Mean Isoform TPM"
    ) +
    scale_fill_manual(
        values = c(
            "Case" = rgb(210 / 255, 53 / 255, 43 / 255),
            "Control" = "green",
            "NonCase" = rgb(72 / 255, 126 / 255, 179 / 255)
        ),
        labels = c(
            "Case" = "ASD LoF Target",
            "Control" = "Control LoF Target",
            "NonCase" = "Untargeted by ASD LoF"
        )
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        legend.title = element_blank()
    )
saveRDS(
    tpmWilcox_compare_plot,
    "data/isoforms/figures/SNVTargets_tpmWilcoxCompare.rds"
)

# Region-specific
# Cortex vs else

cx <- c(
    "OFC", "DFC", "VFC", "MFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C",
    "FC", "PC", "TC", "OC" 
)
cx_samples <- metadata %>%
    filter(Regioncode %in% cx) %>%
    pull(Sample)
notCx_samples <- metadata %>%
    filter(! Regioncode %in% cx) %>%
    pull(Sample)

names(all_features)
mean_tpm_CXvElse_byPeriod <- expand.grid(
    feature_type = names(all_features),
    Period = c(sort(unique(metadata$Period)), "Prenatal", "Postnatal"),
    Cortex = c(TRUE, FALSE)
) %>%
    bind_cols(
        apply(., 1, function(param) {
            print(param)
            ft <- as.character(param[["feature_type"]])
            period <- as.character(param[["Period"]])
            cortex <- as.logical(param[["Cortex"]])
            region_samples <- if (cortex) {cx_samples} else {notCx_samples}
            period_samples <- metadata %>%
                filter(
                    if (period %in% c("Prenatal", "Postnatal")) {
                        Prenatal == (period == "Prenatal")
                    } else {
                        Period == as.numeric(period)
                    }
                ) %>%
                pull(Sample)
            samples <- intersect(region_samples, period_samples)
            expr <- itpm[all_features[[ft]], samples, drop = FALSE] %>%
                as.data.frame() %>%
                bind_cols(
                    data.frame(
                        MeanTPM = rowMeans(.),
                        ensembl_transcript_id = rownames(.)
                    )
                ) %>%
                dplyr::select(MeanTPM) %>%
                nest(MeanTPM, .key = "MeanTPM")
        }) %>%
            bind_rows()
    )

mean_tpm_CXvElse_byPeriod_wilcox <- lapply(
    setNames(nm = sort(unique(as.character(mean_tpm_CXvElse_byPeriod$Period)))),
    function(period) {
        wilcox.test(
            mean_tpm_CXvElse_byPeriod %>%
                filter(feature_type == "Case") %>%
                filter(Period == period) %>%
                filter(Cortex) %>%
                unnest(MeanTPM) %>%
                pull(MeanTPM) %>%
                as.numeric(),
            mean_tpm_CXvElse_byPeriod %>%
                filter(feature_type == "Case") %>%
                filter(Period == period) %>%
                filter(!Cortex) %>%
                unnest(MeanTPM) %>%
                pull(MeanTPM) %>%
                as.numeric(),
            paired = TRUE,
            alternative = "greater"
        ) %>%
            broom::tidy() %>%
            mutate(Period = period) %>%
            mutate(MeanCortex = mean(mean_tpm_CXvElse_byPeriod %>%
                                         filter(feature_type == "Case") %>%
                                         filter(Period == period) %>%
                                         filter(Cortex) %>%
                                         unnest(MeanTPM) %>%
                                         pull(MeanTPM) %>%
                                         as.numeric())) %>%
            mutate(MeanNotCortex = mean(mean_tpm_CXvElse_byPeriod %>%
                                         filter(feature_type == "Case") %>%
                                         filter(Period == period) %>%
                                         filter(!Cortex) %>%
                                         unnest(MeanTPM) %>%
                                         pull(MeanTPM) %>%
                                         as.numeric()))
    }
) %>% 
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    mutate(
        Star = sapply(adj.P.Val, function(p) {
            if (p <= 0.05) {
                if (p <= 0.01) {
                    if (p <= 0.001) {
                        return("*\n*\n*")
                    }
                    return("*\n*")
                }
                return("*")
            }
            return("")
        })
    ) %>%
    mutate(
        Period = factor(Period, levels = c(seq(2, 13), "Prenatal", "Postnatal"))
    ) %>%
    mutate(
        x = as.numeric(Period) - 0.25,
        xend = as.numeric(Period) + 0.25,
        y_horiz = mapply(max, MeanCortex, MeanNotCortex) * 1.05,
        y_star = mapply(max, MeanCortex, MeanNotCortex) * 1.1
    ) %>%
    gather(MeanCortex, MeanNotCortex, key = "ValueType", value = "Mean") %>%
    ungroup() %>%
    mutate(y_ver = Mean * 1.01) %>%
    mutate(x_ver = mapply(
        function(vt, x, xe) {
            if (vt == "MeanCortex") {return(x)} else {return(xe)}
        },
        ValueType, x, xend
    ))

mean_tpm_CXvElse_byPeriod_plot <- ggplot(
    data = mean_tpm_CXvElse_byPeriod_wilcox,
    mapping = aes(
        x = Period, y = Mean, fill = ValueType
    )
) +
    geom_bar(stat = "identity", position = "dodge", colour = "black") +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(
            x = x, xend = xend, y = y_horiz, yend = y_horiz
        )
    ) +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(
            x = x_ver, xend = x_ver, y = y_horiz, yend = y_ver
        )
    ) +
    geom_text(
        inherit.aes = FALSE,
        mapping = aes(
            x = Period, y = y_star, label = Star
        ),
        lineheight = 0.3, size = 8
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
saveRDS(
    mean_tpm_CXvElse_byPeriod_plot, 
    "data/isoforms/figures/SNVTargets_tpmWilcoxCompareCortex.rds"
)
ggsave(
    filename = "data/isoforms/figures/SNVTargets_tpmWilcoxCompareCortex.pdf",
    mean_tpm_CXvElse_byPeriod_plot, 
    device = "pdf", width = 16, height = 12
)
    