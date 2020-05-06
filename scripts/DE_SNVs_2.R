library(tidyverse)
library(biomaRt)
library(ggrepel)

lof <- "acceptor|donor|stop_gain|start_lost|frameshift"

mutations <- read_tsv("data/SNVs/allSNVs_Masterfile_REACH_SSC.tsv")
mutations_expand <- mutations %>%
    separate_rows(Consequence, sep = ",")
lof_mutations <- mutations_expand %>%
    filter(str_detect(Consequence, lof))

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
                Location, Allele, Gene, Consequence, Phenotype_1ctrl_2case
            ),
        by = c("ensembl_gene_id" = "Gene")
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
                Location, Allele, Gene, Feature, 
                Consequence, Phenotype_1ctrl_2case
            ), 
        by = c("ensembl_transcript_id" = "Feature")
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

################################################################################
# Enrichment of isoform-specific lof targets vs gene-specific targets          #
################################################################################

de_iso_gene_lof_spec_enrichment <- apply(
    expand.grid(
        Contrast = sort(unique(tt_combined$Contrast)),
        Phenotype = c("Control", "Case"),
        Specific = c(TRUE)
    ), 
    1,
    function(param) {
        contrast <- as.character(param[["Contrast"]])
        phenotype <- ifelse(as.character(param[["Phenotype"]] == "Case"), 2, 1)
        spec <- as.logical(param[["Specific"]])
        de_iso <- tt_combined %>%
            filter(Data == "Isoform") %>%
            filter(Contrast == contrast) %>%
            filter(Significant) %>%
            filter(Specific == spec) %>%
            pull(ensembl_gene_id) %>%
            unique()
        de_genes <- tt_combined %>%
            filter(Data == "Gene") %>%
            filter(Contrast == contrast) %>%
            filter(Significant) %>%
            filter(Specific == spec) %>%
            pull(ensembl_gene_id) %>%
            unique()
        de_iso_targets <- tt_combined %>%
            filter(Data == "Isoform") %>%
            filter(Contrast == contrast) %>%
            filter(Significant) %>%
            filter(Specific == spec) %>%
            filter(!is.na(Consequence)) %>%
            filter(Phenotype_1ctrl_2case == phenotype) %>%
            pull(ensembl_gene_id) %>%
            unique()
        de_gene_targets <- tt_combined %>%
            filter(Data == "Gene") %>%
            filter(Contrast == contrast) %>%
            filter(Significant) %>%
            filter(Specific == spec) %>%
            filter(!is.na(Consequence)) %>%
            filter(Phenotype_1ctrl_2case == phenotype) %>%
            pull(ensembl_gene_id) %>%
            unique()
        fisher.test(
            matrix(c(
                length(de_iso_targets),
                length(de_gene_targets),
                length(de_iso[!de_iso %in% de_iso_targets]),
                length(de_genes[!de_genes %in% de_gene_targets])
            ), byrow = TRUE, ncol = 2, nrow = 2)
        ) %>%
            broom::tidy() %>%
            mutate(
                Contrast = contrast,
                Phenotype = as.character(param[["Phenotype"]]),
                Specific = as.logical(param[["Specific"]]),
                Overlap_Iso = length(de_iso_targets),
                Overlap_Genes = length(de_gene_targets)
            )
    }
) %>%
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni"))

barplot_input <- de_iso_gene_lof_spec_enrichment %>%
    dplyr::select(
        adj.P.Val, Contrast, Phenotype, Specific, starts_with("Overlap")
    ) %>%
    reshape2::melt(
        id.vars = c("Contrast", "Phenotype", "Specific", "adj.P.Val")
    ) %>%
    rename(Overlap_type = variable, Overlap = value) %>%
    filter(Phenotype == "Case", Specific) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.05, "*", "")) %>%
    mutate(Contrast = factor(Contrast)) %>%
    mutate(x = 0.75 + as.numeric(Contrast) - 1) %>%
    mutate(xend = 1.25 + as.numeric(Contrast) - 1) %>%
    mutate(
        y = mapply(
            function(contrast, star) {
                if (star == "") { return(NA) }
                else {
                    this_y<- filter(., Contrast == contrast) %>%
                        pull(Overlap) %>%
                        max()
                    this_y + 5
                }
            },
            Contrast, Star
        )
    ) %>%
    mutate(
        Overlap_type = factor(
            Overlap_type, levels = c("Overlap_Genes", "Overlap_Iso")
        )
    ) %>%
    mutate(
        Data = ifelse(str_detect(Overlap_type, "Gene"), "Gene", "Isoform")
    ) %>%
    mutate(
        Proportion = mapply(
            function(ctr, dat, overlap) {
                total <- tt_combined %>%
                    filter(Data == dat & Contrast == ctr) %>%
                    filter(Specific & Significant) %>%
                    nrow()
                overlap / total
            },
            Contrast, Data, Overlap
        )
    )

de_iso_gene_lof_spec_enrichment_barplot <- ggplot(
    data = barplot_input,
    mapping = aes(
        x = Contrast, y = Proportion, fill = Data
    )
) +
    geom_bar(
        stat = "identity", position = "dodge"
    ) +
    # geom_segment(
    #     inherit.aes = FALSE,
    #     mapping = aes(x = x, y = y, xend = xend, yend = y)
    # ) +
    geom_text(
        inherit.aes = FALSE,
        mapping = aes(
            x = Contrast, y = y + 1, label = Star
        ),
        size = 10
    ) +
    labs(
        title = "Specific DE genes targeted by Case LoF Mutations",
        subtitle = expression(
            "Fisher-exact test, Bonferroni Corrected, *P">="0.05"
        ),
        y = "Proportion DE Genes"
    ) +
    guides(
        fill = guide_legend(
            title = ""
        )
    ) +
    # scale_fill_discrete(
    #     breaks = c("Overlap_Genes", "Overlap_Iso"),
    #     labels = c(
    #         "Overlap_Genes" = "Gene",
    #         "Overlap_Iso" = "Isoform\n(Summarized to Genes)"
    #     )
    # ) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )
de_iso_gene_lof_spec_enrichment_barplot
saveRDS(
    de_iso_gene_lof_spec_enrichment_barplot,
    "data/figures/DE_Enr_LoFTargets_Bar.rds"
)
ggsave(
    filename = "data/figures/DE_Enr_LoFTargets_Bar.pdf",
    plot = de_iso_gene_lof_spec_enrichment_barplot,
    device = "pdf", width = 16, height = 12
)

################################################################################
# LoF Composition per Contrast                                                 #
################################################################################

lof_composition <- ggplot(
    data = tt_combined %>%
        filter(Phenotype_1ctrl_2case == 2) %>%
        filter(Significant & Specific) %>%
        dplyr::select(
            Location, Consequence, Data, Contrast
        ) %>%
        distinct() %>%
        count(Contrast, Consequence, Data) %>%
        mutate(
            total_n = mapply(
                function(ctr, dat) {
                    tt_combined %>%
                        filter(Phenotype_1ctrl_2case == 2) %>%
                        filter(Significant & Specific) %>%
                        filter(Data == dat & Contrast == ctr) %>%
                        nrow()
                },
                Contrast, Data
            )
        ) %>%
        mutate(rel_n = n / total_n) %>%
        mutate(
            Consequence = sapply(Consequence, function(con) {
                new_labels <- c(
                    "frameshift_variant" = "Frameshift",
                    "splice_acceptor_variant" = "Splice\nAcceptor",
                    "splice_donor_variant" = "Splice\nDonor",
                    "start_lost" = "Start Loss",
                    "stop_gained" = "Stop Gain"
                )
                new_labels[[con]]
            })
        ),
    mapping = aes(
        x = Contrast, y = rel_n, fill = Data, group = Data
    )
) +
    facet_grid(Consequence ~ .) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        y = "Relative Consequence Frequency"
    ) +
    scale_y_continuous(
        limits = c(0, 1.2), breaks = c(0, 0.5, 1), expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 20),
        panel.spacing = unit(0.25, "lines")
    )
lof_composition
saveRDS(
    lof_composition,
    "data/figures/LOF_DE_Breakdown.rds"
)

################################################################################
# Volcano Plots                                                                #
################################################################################

lof_volcano <- ggplot(
    data = tt_combined
) +
    facet_wrap(~ Contrast) +
    geom_point(
        data = tt_combined %>%
            filter(
                ! (Specific & Significant & 
                    Phenotype_1ctrl_2case == 2 & !is.na(Phenotype_1ctrl_2case))
            ),
        mapping = aes(x = logFC, y = -log10(adj.P.Val)),
        colour = "lightgrey", alpha = 0.1
    ) +
    geom_point(
        data = tt_combined %>%
            filter(
                Specific & Significant & 
                    Phenotype_1ctrl_2case == 2 & !is.na(Phenotype_1ctrl_2case)
            ),
        mapping = aes(
            x = logFC, y = -log10(adj.P.Val),
            colour = Data
        )
    ) +
    geom_vline(xintercept = c(log2(1.5), -log2(1.5)), colour = "black") +
    geom_hline(yintercept = -log10(0.05), colour = "black") +
    labs(
        title = "Data-specific LoF Targets",
        x = expression("log"["2"]*"(Fold Change)"),
        y = expression("-log"["10"]*"(FDR)")
    ) +
    guides(
        colour = guide_legend(title = ""),
        alpha = FALSE
    ) +
    scale_y_continuous(
        trans = "log1p", breaks = c(0, 25, 50, 100, 200)
    ) +
    scale_colour_manual(
        values = c(
            "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
            "Isoform" = rgb(84 / 255, 188 / 255, 194 / 255)
        )
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
saveRDS(
    lof_volcano,
    "data/figures/DE_LOF_Volcano.rds"
)
ggsave(
    filename = "data/figures/DE_LOF_Volcano.pdf",
    plot = lof_volcano,
    device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "data/figures/DE_LOF_Volcano.png",
    plot = lof_volcano,
    device = "png", width = 16, height = 12, units = "in"
)

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]


tt_combined %>%
    filter(
        Specific & Significant & 
            Phenotype_1ctrl_2case == 2 & !is.na(Phenotype_1ctrl_2case)
    ) %>%
    dplyr::select(
        Contrast, external_gene_id, Data
    ) %>%
    distinct()

lof_de_overlaps <- expand.grid(
    Contrast1 = sort(unique(tt_combined$Contrast)),
    Contrast2 = sort(unique(tt_combined$Contrast)),
    Data = c("Gene", "Isoform"), 
    stringsAsFactors = FALSE
) %>%
    mutate(
        Overlap = mapply(
            function(ctr1, ctr2, dat) {
                if(ctr1 == ctr2) { 
                    tt_combined %>%
                        filter(
                            Specific & Significant &
                                Phenotype_1ctrl_2case == 2 &
                                !is.na(Phenotype_1ctrl_2case)
                        ) %>%
                        dplyr::select(
                            Contrast, external_gene_id, Data
                        ) %>%
                        distinct() %>%
                        filter(Data == dat) %>%
                        filter(Contrast %in% c(ctr1, ctr2)) %>%
                        distinct() %>%
                        nrow()
                } else {
                    tt_combined %>%
                        filter(
                            Specific & Significant &
                                Phenotype_1ctrl_2case == 2 &
                                !is.na(Phenotype_1ctrl_2case)
                        ) %>%
                        dplyr::select(
                            Contrast, external_gene_id, Data
                        ) %>%
                        distinct() %>%
                        filter(Data == dat) %>%
                        filter(Contrast %in% c(ctr1, ctr2)) %>%
                        distinct() %>%
                        count(external_gene_id) %>%
                        filter(n > 1) %>%
                        nrow()
                }
            },
            Contrast1, Contrast2, Data
        )
    ) %>%
    mutate(
        Unique1 = mapply(
            function(ctr1, ctr2, dat) {
                if(ctr1 == ctr2) { return(0) }
                tt_combined %>%
                    filter(
                        Specific & Significant & 
                            Phenotype_1ctrl_2case == 2 & 
                            !is.na(Phenotype_1ctrl_2case)
                    ) %>%
                    dplyr::select(
                        Contrast, external_gene_id, Data
                    ) %>%
                    distinct() %>%
                    filter(Data == dat) %>%
                    filter(Contrast %in% c(ctr1, ctr2)) %>%
                    group_by(external_gene_id) %>%
                    filter(n() == 1) %>%
                    ungroup() %>%
                    filter(Contrast == ctr1) %>%
                    nrow()
            }, 
            Contrast1, Contrast2, Data
        )
    ) %>%
    mutate(
        Unique2 = mapply(
            function(ctr1, ctr2, dat) {
                if(ctr1 == ctr2) { return(0) }
                tt_combined %>%
                    filter(
                        Specific & Significant &
                            Phenotype_1ctrl_2case == 2 &
                            !is.na(Phenotype_1ctrl_2case)
                    ) %>%
                    dplyr::select(
                        Contrast, external_gene_id, Data
                    ) %>%
                    distinct() %>%
                    filter(Data == dat) %>%
                    filter(Contrast %in% c(ctr1, ctr2)) %>%
                    group_by(external_gene_id) %>%
                    filter(n() == 1) %>%
                    ungroup() %>%
                    filter(Contrast == ctr2) %>%
                    nrow()
            },
            Contrast1, Contrast2, Data
        )
    ) %>%
    mutate(
        Jaccard = Overlap / (Unique1 + Unique2 + Overlap)
    )

lof_de_overlap_plot <- ggplot(
    data = lof_de_overlaps %>%
        filter(Jaccard > 0 & Contrast1 != Contrast2),
    mapping = aes(
        x = Contrast1, y = Contrast2, size = Jaccard, colour = Data
    )
) +
    geom_point(
        data = lof_de_overlaps %>%
            filter(Jaccard > 0 & Contrast1 != Contrast2) %>%
            filter(Data == "Gene") %>%
            mutate(Contrast1_order = as.numeric(factor(Contrast1))) %>%
            mutate(Contrast2_order = as.numeric(factor(Contrast2))) %>%
            filter(Contrast1_order < Contrast2_order)
    ) +
    geom_point(
        data = lof_de_overlaps %>%
            filter(Jaccard > 0 & Contrast1 != Contrast2) %>%
            filter(Data == "Isoform") %>%
            mutate(Contrast1_order = as.numeric(factor(Contrast1))) %>%
            mutate(Contrast2_order = as.numeric(factor(Contrast2))) %>%
            filter(Contrast1_order > Contrast2_order)
    ) +
    scale_size_continuous(
        range = c(0, 15), limits = c(0, 1)
    ) +
    theme_bw() + 
    theme(
        text = element_text(size = 24),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 24, hjust = 1)
    )
lof_de_overlap_plot
saveRDS(
    lof_de_overlap_plot,
    "data/figures/LOF_DE_Overlaps.rds"
)
ggsave(
    filename = "data/figures/LOF_DE_Overlaps.pdf",
    plot = lof_de_overlap_plot,
    device = "pdf", width = 16, height = 12
)
