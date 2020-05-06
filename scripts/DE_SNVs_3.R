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

chisq_lof_target_test <- lapply(
    setNames(nm = sort(unique(tt_combined$Contrast))),
    function(contrast) {
        dei_spec <- tt_combined %>%
            filter(Contrast == contrast & Data == "Isoform") %>%
            filter(Significant & Specific) %>%
            pull(ensembl_transcript_id) %>%
            unique() %>%
            length()
        deg_spec <- tt_combined %>%
            filter(Contrast == contrast & Data == "Gene") %>%
            filter(Significant & Specific) %>%
            pull(ensembl_gene_id) %>%
            unique() %>%
            length()
        dei_spec_targ <- tt_combined %>%
            filter(Contrast == contrast & Data == "Isoform") %>%
            filter(Significant & Specific) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            pull(ensembl_transcript_id)
        deg_spec_targ <- tt_combined %>%
            filter(Contrast == contrast & Data == "Gene") %>%
            filter(Significant & Specific) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            pull(ensembl_gene_id)
        chisq.test(
            x = c(length(dei_spec_targ), length(deg_spec_targ)),
            p = c(
                dei_spec / (dei_spec + deg_spec),
                deg_spec / (dei_spec + deg_spec)
            ),
            simulate.p.value = TRUE
        ) %>%
            broom::tidy() %>%
            mutate(
                Contrast = contrast,
                DEITargets = length(dei_spec_targ),
                DEGTargets = length(deg_spec_targ),
                DEISpec = dei_spec,
                DEGSpec = deg_spec
            ) %>%
            mutate(
                DEITargets_Prop = DEITargets / DEISpec,
                DEGTargets_Prop = DEGTargets / DEGSpec,
            )
    }
) %>%
    bind_rows() %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH"))

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
# Fold change case vs control targets                                          #
################################################################################

ggplot(data = tt_combined) +
    facet_grid(. ~ Contrast) +
    # geom_histogram(
    #     data = tt_combined %>%
    #         filter(logFC > 0) %>%
    #         filter(Phenotype_1ctrl_2case == 1) %>%
    #         filter(Significant & Specific),
    #     mapping = aes(
    #         x = abs(logFC), y = ..count.., fill = Phenotype_1ctrl_2case
    #     ),
    #     alpha = 0.8
    # ) +
    # geom_histogram(
    #     data = tt_combined %>%
    #         filter(logFC < 0) %>%
    #         filter(Phenotype_1ctrl_2case == 1) %>%
    #         filter(Significant & Specific),
    #     mapping = aes(
    #         x = abs(logFC), y = -..count.., fill = Phenotype_1ctrl_2case
    #     ),
    #     alpha = 0.8
    # ) +
    geom_histogram(
        data = tt_combined %>%
            filter(logFC > 0) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            filter(Data == "Gene") %>%
            filter(Significant & Specific),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        alpha = 0.8
    ) +
    geom_histogram(
        data = tt_combined %>%
            filter(logFC < 0) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            filter(Data == "Gene") %>%
            filter(Significant & Specific),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        alpha = 0.8
    ) +
    geom_histogram(
        data = tt_combined %>%
            filter(logFC > 0) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            filter(Data == "Isoform") %>%
            filter(Significant & Specific),
        mapping = aes(
            x = abs(logFC), y = ..count.., fill = Data
        ),
        alpha = 0.8
    ) +
    geom_histogram(
        data = tt_combined %>%
            filter(logFC < 0) %>%
            filter(Phenotype_1ctrl_2case == 2) %>%
            filter(Data == "Isoform") %>%
            filter(Significant & Specific),
        mapping = aes(
            x = abs(logFC), y = -..count.., fill = Data
        ),
        alpha = 0.8
    )

################################################################################
# Are case mutations more likely to hit specifically DEI vs DEG rel. to ctr    #
################################################################################

lapply(
    setNames(nm = sort(unique(tt_combined$Contrast))),
    function(contrast) {
        fisher.test(
            x = matrix(c(
                
            ))
        )
    }
)