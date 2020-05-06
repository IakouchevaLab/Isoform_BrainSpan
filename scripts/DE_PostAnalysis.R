library(tidyverse)
library(biomaRt)
library(VennDiagram)

metadata <- read_tsv("data/metadata.tsv")
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
    left_join(anno_expand, by = "ensembl_gene_id")
tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")

# Compare DE iso against genes

de_gene_count <- tt_genes %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    distinct(ensembl_gene_id, Contrast) %>%
    count(Contrast) %>%
    # rename(DE_Genes = n) %>%
    mutate(Data = "Gene")
de_iso_count <- tt_iso %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    distinct(Contrast, ensembl_gene_id) %>%
    count(Contrast) %>%
    # rename(DE_Isoforms = n) %>%
    mutate(Data = "Isoform\n(Summarized to Genes)")
de_counts <- bind_rows(de_gene_count, de_iso_count)
de_frequency_combined_plot <- ggplot(
    data = de_counts,
    mapping = aes(
        x = Contrast, y = n
    )
) +
    geom_bar(
        mapping = aes(fill = Data, group = Data), 
        stat = "identity", position = "dodge"
    ) +
    labs(
        y = "DE Frequency"
    ) +
    scale_y_continuous(trans = "log1p", breaks = c(100, 1000, 10000, 20000)) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )

dte_gn <- tt_genes %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5))
dte_iso <- tt_iso %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5))

dte_overlaps <- data.frame(
    Contrast = sort(unique(tt_genes$Contrast))
) %>%
    mutate(overlap = sapply(Contrast, function(ctr) {
        length(intersect(
            dte_gn %>%
                filter(Contrast == ctr) %>%
                pull(ensembl_gene_id),
            dte_iso %>%
                filter(Contrast == ctr) %>%
                pull(ensembl_gene_id)
        ))
    }))

de_frequency_combined_plot <- de_frequency_combined_plot +
    geom_bar(
        data = dte_overlaps,
        mapping = aes(
            x = Contrast, y = overlap
        ),
        stat = "identity", fill = "black", colour = NA, alpha = 0.5
    )
saveRDS(
    de_frequency_combined_plot,
    "data/figures/DEFrequencyCombinedPlot.rds"
)







de_frequency_combined_plot_untrans <- ggplot(
    data = de_counts,
    mapping = aes(
        x = Contrast, y = n
    )
) +
    geom_bar(
        mapping = aes(fill = Data, group = Data), 
        stat = "identity", position = "dodge"
    ) +
    labs(
        y = "DE Frequency"
    ) +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
de_frequency_combined_plot_untrans <- de_frequency_combined_plot_untrans +
    geom_bar(
        data = dte_overlaps,
        mapping = aes(
            x = Contrast, y = overlap
        ),
        stat = "identity", fill = "black", colour = NA, alpha = 0.5
    )
saveRDS(
    de_frequency_combined_plot_untrans,
    "data/figures/DEFrequencyCombinedPlot_untrans.rds"
)




de_frequency_combined_plot_untrans_noPrePost <- ggplot(
    data = de_counts %>%
        filter(Contrast != "PrePost") %>%
        bind_rows(
            data.frame(
                Contrast = c("PrePost", "PrePost"),
                n = c(1, 1),
                Data = c("Gene", "Isoform\n(Summarized to Genes)")
            )
        ),
    mapping = aes(
        x = Contrast, y = n
    )
) +
    geom_bar(
        mapping = aes(fill = Data, group = Data), 
        stat = "identity", position = "dodge"
    ) +
    labs(
        y = "DE Frequency"
    ) +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
de_frequency_combined_plot_untrans_noPrePost <- de_frequency_combined_plot_untrans_noPrePost +
    geom_bar(
        data = dte_overlaps %>%
            filter(Contrast != "PrePost"),
        mapping = aes(
            x = Contrast, y = overlap
        ),
        stat = "identity", fill = "black", colour = NA, alpha = 0.5
    )
saveRDS(
    de_frequency_combined_plot_untrans_noPrePost,
    "data/figures/DEFrequencyCombinedPlot_untrans_noPrePost.rds"
)























# Specificity of DE for isoform or gene analysis against Satterstrom ASD

satterstromASD <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

lapply(
    setNames(
        nm = sort(unique(tt_genes$Contrast))
    ),
    function(ctr) {
        dte_gene <- dte_gn %>%
            filter(Contrast == ctr) %>%
            pull(external_gene_id)
        dte_iso_to_gene <- dte_iso %>%
            filter(Contrast == ctr) %>%
            pull(external_gene_id)
        unique_to_gene <- dte_gene[!dte_gene %in% dte_iso_to_gene]
        unique_to_iso <- dte_iso_to_gene[!dte_iso_to_gene %in% dte_gene]
        asd_gene <- length(intersect(unique_to_gene, satterstromASD))
        asd_iso <- length(intersect(unique_to_iso, satterstromASD))
        data.frame(
            Contrast = ctr,
            GeneASDOverlap = asd_gene,
            IsoformASDOverlap = asd_iso,
            GeneASDOverlapNorm = asd_gene / length(unique_to_gene),
            IsoformASDOverlapNorm = asd_iso / length(unique_to_iso)
        )
    }
) %>% 
    bind_rows() %>%
    write_tsv("data/DESpecificOverlapWithSatterstromASD.tsv")

# Test distribution of effect sizes

eff_size_test <- lapply(
    setNames(nm = c(sort(unique(tt_iso$Contrast)), "Overall")),
    function(contrast) {
        if (contrast != "Overall") {
            res <- t.test(
                tt_iso %>%
                    filter(adj.P.Val <= 0.05) %>%
                    filter(Contrast == contrast) %>%
                    pull(logFC) %>%
                    abs(),
                tt_genes %>%
                    filter(adj.P.Val <= 0.05) %>%
                    filter(Contrast == contrast) %>%
                    pull(logFC) %>%
                    abs(),
                # alternative = "greater"
            )
            mean_iso <- tt_iso %>%
                filter(adj.P.Val <= 0.05) %>%
                filter(Contrast == contrast) %>%
                pull(logFC) %>%
                abs() %>%
                mean()
            mean_gene <- tt_genes %>%
                filter(adj.P.Val <= 0.05) %>%
                filter(Contrast == contrast) %>%
                pull(logFC) %>%
                abs() %>%
                mean()
        } else {
            res <- t.test(
                tt_iso %>%
                    filter(adj.P.Val <= 0.05) %>%
                    pull(logFC) %>%
                    abs(),
                tt_genes %>%
                    filter(adj.P.Val <= 0.05) %>%
                    pull(logFC) %>%
                    abs(),
                # alternative = "greater"
            )
            mean_iso <- tt_iso %>%
                filter(adj.P.Val <= 0.05) %>%
                pull(logFC) %>%
                abs() %>%
                mean()
            mean_gene <- tt_genes %>%
                filter(adj.P.Val <= 0.05) %>%
                pull(logFC) %>%
                abs() %>%
                mean()
        }
        res %>%
            broom::tidy() %>%
            mutate(Contrast = contrast) %>%
            mutate(
                Mean_Iso = mean_iso,
                Mean_Gene = mean_gene
            )
    }
) %>%
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, "bonferroni")) %>%
    dplyr::select(
        Contrast, estimate, Mean_Iso, Mean_Gene, p.value, adj.P.Val
    )
t.test(
    tt_iso %>%
        filter(adj.P.Val <= 0.05) %>%
        pull(abs(logFC)), 
    tt_genes %>%
        filter(adj.P.Val <= 0.05) %>%
        pull(abs(logFC)),
    alternative = "greater"
)
tt_iso %>%
    filter(adj.P.Val <= 0.05) %>%
    pull(abs(logFC)) %>%
    mean()
tt_genes %>%
    filter(adj.P.Val <= 0.05) %>%
    pull(abs(logFC)) %>%
    mean()

lfc_dist_input <- bind_rows(
    tt_genes %>% mutate(Data = "Gene") %>% 
        mutate(biotype = gene_biotype) %>%
        dplyr::select(
            ensembl_gene_id, adj.P.Val, logFC, 
            Contrast, biotype, Data
        ) %>%
        rename(Feature = ensembl_gene_id) %>%
        distinct(),
    tt_iso %>% 
        mutate(Data = "Isoform") %>% 
        mutate(biotype = transcript_biotype) %>%
        dplyr::select(
            ensembl_transcript_id, adj.P.Val, logFC, 
            Contrast, biotype, Data
        ) %>%
        rename(Feature = ensembl_transcript_id) %>%
        distinct()
) %>% 
    filter(adj.P.Val <= 0.05) %>%
    filter(biotype == "protein_coding") %>%
    mutate(
        Star = ifelse(
            Contrast %in% pull(
                filter(eff_size_test, adj.P.Val <= 0.05), 
                Contrast
            ), "*", ""
        )
    ) %>%
    mutate(Contrast = paste0(Contrast, Star))

lfc_dist <- ggplot() +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Isoform"),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Gene"),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Isoform"),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Gene"),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    facet_grid(. ~ Contrast) +
    labs(
        x = expression("abs(log"["10"]*"Fold Change)"),
        y = expression(atop(
            "Significant (FDR"<="0.05)",
            "DE Features"
        )),
        title = "DE Effect Size Distribution",
        subtitle = expression("T-Test Bonferroni-Corrected *P"<="0.05")
    ) +
    guides(fill = guide_legend(title = "")) +
    scale_x_continuous(
        trans = "log1p", breaks  = c(0, 1, 5)
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 18),
        panel.spacing = unit(0, "lines"),
        legend.justification = c(1, 1),
        legend.position = c(0.15, 0.35),
        legend.title = element_blank()
    )
lfc_dist

ggsave(
    filename = "data/figures/LFC_DIST.pdf",
    plot = lfc_dist,
    device = "pdf", width = 16, height = 5
)

lfc_dist_overall <- ggplot() +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Isoform"),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Gene"),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Isoform"),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Gene"),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    labs(
        x = expression("abs(log"["2"]*"Fold Change)"),
        y = expression(atop(
            "Significant (FDR"<="0.05)",
            "DE Features"
        )),
        title = "DE Effect Size Distribution, Overall",
        subtitle = expression("T-Test Bonferroni-Corrected P"<="0.05")
    ) +
    guides(fill = guide_legend(title = "")) +
    scale_x_continuous(
        trans = "log1p", breaks  = c(0, 1, 5)
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 18),
        panel.spacing = unit(0, "lines"),
        legend.justification = c(1, 1),
        legend.position = c(0.95, 0.95),
        legend.title = element_blank()
    )
lfc_dist_overall

ggsave(
    filename = "data/figures/LFC_DIST_OVERALL.pdf",
    plot = lfc_dist_overall,
    device = "pdf", width = 16, height = 5
)

de_iso <- tt_iso %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
    mutate(Feature = ensembl_transcript_id)
de_gene <- tt_genes %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
    mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
    mutate(Feature = ensembl_gene_id)
de_iso_spec <- de_iso %>%
    filter(! selector %in% de_gene$selector)
de_gene_spec <- de_gene %>%
    filter(! selector %in% de_iso$selector)

lfc_dist_allAndSpec_input <- bind_rows(
    de_iso_spec %>%
        mutate(Data = "Isoform") %>%
        mutate(Class = "Specific") %>%
        rename(biotype = transcript_biotype) %>%
        dplyr::select(Feature, logFC, Data, Class, biotype, Contrast) %>%
        distinct(),
    de_gene_spec %>%
        mutate(Data = "Gene") %>%
        mutate(Class = "Specific") %>%
        rename(biotype = gene_biotype) %>%
        dplyr::select(Feature, logFC, Data, Class, biotype, Contrast) %>%
        distinct(),
    lfc_dist_input %>%
        mutate(Class = "All") %>%
        mutate(Contrast = gsub("\\*", "", Contrast))
) %>%
    filter(biotype == "protein_coding")
sig_test <- expand.grid(
    Contrast = sort(unique(lfc_dist_allAndSpec_input$Contrast)),
    Class = sort(unique(lfc_dist_allAndSpec_input$Class))
) %>%
    apply(
        ., 1, function(param) {
            iso <- lfc_dist_allAndSpec_input %>%
                filter(Data == "Isoform") %>%
                filter(
                    Contrast == param[["Contrast"]] & Class == param[["Class"]]
                )
            gene <- lfc_dist_allAndSpec_input %>%
                filter(Data == "Gene") %>%
                filter(
                    Contrast == param[["Contrast"]] & Class == param[["Class"]]
                )
            data.frame(
                Contrast = param[["Contrast"]],
                Class = param[["Class"]],
                p.value = tryCatch({
                    t.test(
                        abs(iso$logFC), abs(gene$logFC), alternative = "greater"
                    )$p.value
                }, error = function(e) {
                    return(1)
                })
            )
        }
    ) %>%
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni"))
sig_test <- bind_rows(
    sig_test %>% mutate(Data = "Gene"),
    sig_test %>% mutate(Data = "Isoform")
) %>%
    mutate(
        mean_eff = mapply(
            function(contrast, class_type, data_type) {
                mean(abs(pull(filter(
                    lfc_dist_allAndSpec_input, 
                    Contrast == contrast & Class == class_type & Data == data_type
                ), logFC)))
            },
            Contrast, Class, Data
        )
    )
lfc_dist_allAndSpec <- ggplot() +
    geom_vline(
        data = sig_test %>%
            filter(Class == "Specific") %>%
            mutate(
                Contrast = mapply(
                    function(ctr, cl) {
                        pval <- sig_test %>%
                            filter(Contrast == ctr & Class == cl) %>%
                            pull(adj.P.Val) %>%
                            max()
                        return(ifelse(pval <= 0.05, paste0(ctr, "*"), ctr))
                    },
                    Contrast, Class
                )
            ),
        mapping = aes(xintercept = mean_eff, colour = Data),
        show.legend = FALSE
    ) +
    # geom_rect(
    #     data = sig_test %>%
    #         filter(adj.P.Val <= 0.05),
    #     xmin = -1, xmax = max(lfc_dist_allAndSpec_input$logFC),
    #     ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1
    # ) +
    geom_histogram(
        data = lfc_dist_allAndSpec_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Isoform") %>%
            filter(Class == "Specific") %>%
            mutate(
                Contrast = mapply(
                    function(ctr, cl) {
                        pval <- sig_test %>%
                            filter(Contrast == ctr & Class == cl) %>%
                            pull(adj.P.Val) %>%
                            max()
                        return(ifelse(pval <= 0.05, paste0(ctr, "*"), ctr))
                    },
                    Contrast, Class
                )
            ),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_allAndSpec_input %>%
            filter(logFC > 0) %>%
            filter(Data == "Gene") %>%
            filter(Class == "Specific") %>%
            mutate(
                Contrast = mapply(
                    function(ctr, cl) {
                        pval <- sig_test %>%
                            filter(Contrast == ctr & Class == cl) %>%
                            pull(adj.P.Val) %>%
                            max()
                        return(ifelse(pval <= 0.05, paste0(ctr, "*"), ctr))
                    },
                    Contrast, Class
                )
            ),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_allAndSpec_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Isoform") %>%
            filter(Class == "Specific") %>%
            mutate(
                Contrast = mapply(
                    function(ctr, cl) {
                        pval <- sig_test %>%
                            filter(Contrast == ctr & Class == cl) %>%
                            pull(adj.P.Val) %>%
                            max()
                        return(ifelse(pval <= 0.05, paste0(ctr, "*"), ctr))
                    },
                    Contrast, Class
                )
            ),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    geom_histogram(
        data = lfc_dist_allAndSpec_input %>%
            filter(logFC < 0) %>%
            filter(Data == "Gene") %>%
            filter(Class == "Specific") %>%
            mutate(
                Contrast = mapply(
                    function(ctr, cl) {
                        pval <- sig_test %>%
                            filter(Contrast == ctr & Class == cl) %>%
                            pull(adj.P.Val) %>%
                            max()
                        return(ifelse(pval <= 0.05, paste0(ctr, "*"), ctr))
                    },
                    Contrast, Class
                )
            ),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        colour = "black", size = 0.1, alpha = 0.75
    ) +
    facet_grid(. ~ Contrast, scales = "free_y") +
    labs(
        x = expression("abs(log"["2"]*"Fold Change)"),
        y = expression(atop(
            "Significant (FDR"<="0.05)",
            "DE Features"
        )),
        title = "DE Effect Size Distribution",
        subtitle = expression("T-Test Bonferroni-Corrected *P"<="0.05")
    ) +
    guides(fill = guide_legend(title = "")) +
    scale_x_continuous(
        trans = "log1p", breaks  = c(0, 1, 5)
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 18),
        # plot.subtitle = element_text(colour = "red"),
        panel.spacing = unit(0, "lines"),
        # legend.justification = c(1, 1),
        # legend.position = c(0.15, 0.35),
        legend.title = element_blank()
    )
lfc_dist_allAndSpec
ggsave(
    filename = "data/figures/LFC_DIST_SPEC.pdf",
    plot = lfc_dist_allAndSpec,
    device = "pdf", width = 16, height = 5
)

################################################################################
# Final LFC Dist figure                                                        #
################################################################################

lfc_dist_combined <- cowplot::plot_grid(
    lfc_dist_overall +
        labs(
            title = "DE Effect Size Distribution",
            subtitle = expression("T-Test Bonferroni-Corrected *P"<="0.05")
        ) +
        theme(
            axis.title.x = element_blank(),
            legend.justification = c(1, 1),
            legend.position = c(0.95, 0.35)
        ),
    lfc_dist_allAndSpec + 
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            legend.position = "none"
        ),
    ncol = 1, axis = "lr", 
    labels = "AUTO", label_size = 40, label_fontface = "bold"
)
ggsave(
    filename = "data/figures/LFC_DIST_FINAL_COMBINED.pdf",
    plot = lfc_dist_combined,
    device = "pdf", width = 16, height = 9
)

################################################################################
# Regulation Direction Analysis                                                #
################################################################################

t.test(
    x = pull(filter(tt_iso, adj.P.Val <= 0.05), logFC),
    mu = 0
)
t.test(
    x = pull(filter(de_iso_spec, adj.P.Val <= 0.05), logFC),
    mu = 0
)
