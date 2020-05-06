# Functional enrichment and curated list enrichments

library(tidyverse)
library(biomaRt)
library(gProfileR)
source("scripts/utility/plotting.R")

dir.create("data/genes/Enrichments")
dir.create("data/isoforms/Enrichments")

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
            attributes = c("ensembl_gene_id", "start_position", "end_position"),
            filters = "ensembl_gene_id",
            values = unique(.$ensembl_gene_id),
            mart = ensembl
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(gene_length = abs(start_position - end_position))

asd_lists <- lapply(
    setNames(
        nm = readxl::excel_sheets(
            "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx"
        )
    ),
    function(sheet) {
        if (sheet == "Bibliography") {
            return(NULL)
        }
        readxl::read_xlsx(
            "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
            sheet = sheet
        )[[1]]
    }
)
asd_lists <- asd_lists[!sapply(asd_lists, is.null)]

sv_gn <- 13
sv_tx <- 16

################################################################################
# Gene ASD List Enrichment                                                     #
################################################################################

tt_genes <- readRDS(
    paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_gene_id")
de_gene_list_enrichment <- expand.grid(
    Contrast = sort(unique(tt_genes$Contrast)),
    List = sort(unique(names(asd_lists)))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            lit_list <- as.character(param[["List"]])
            all_genes <- tt_genes %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_genes <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            list_genes <- asd_lists[[lit_list]]
            contingency <- matrix(c(
                length(
                    intersect(de_genes, list_genes)
                ),
                length(
                    intersect(de_genes, all_genes[!all_genes %in% list_genes])
                ),
                length(
                    intersect(all_genes[!all_genes %in% de_genes], list_genes)
                ),
                length(all_genes[!all_genes %in% c(de_genes, list_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    List = lit_list
                )
        }
    ) %>%
    bind_rows() 
# %>%
#     group_by(Contrast) %>%
#     mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
#     ungroup()

################################################################################
# Isoform ASD List Enrichment                                                  #
################################################################################

tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")
de_iso_list_enrichment <- expand.grid(
    Contrast = sort(unique(tt_iso$Contrast)),
    List = sort(unique(names(asd_lists)))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            lit_list <- as.character(param[["List"]])
            all_genes <- tt_iso %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_iso_to_gene <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            list_genes <- asd_lists[[lit_list]]
            contingency <- matrix(c(
                length(
                    intersect(de_iso_to_gene, list_genes)
                ),
                length(
                    intersect(
                        de_iso_to_gene, all_genes[!all_genes %in% list_genes]
                    )
                ),
                length(
                    intersect(
                        all_genes[!all_genes %in% de_iso_to_gene], list_genes
                    )
                ),
                length(all_genes[!all_genes %in% c(de_iso_to_gene, list_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    List = lit_list
                )
        }
    ) %>%
    bind_rows() 
# %>%
#     group_by(Contrast) %>%
#     mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
#     ungroup()

################################################################################
# Specific DE Enrichments                                                      #
################################################################################

de_gene_list_enrichment_specific <- expand.grid(
    Contrast = sort(unique(tt_genes$Contrast)),
    List = sort(unique(names(asd_lists)))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            lit_list <- as.character(param[["List"]])
            all_genes <- tt_genes %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_iso_gene_sym <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            de_genes <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                filter(! external_gene_id %in% de_iso_gene_sym) %>%
                pull(external_gene_id)
            list_genes <- asd_lists[[lit_list]]
            contingency <- matrix(c(
                length(
                    intersect(de_genes, list_genes)
                ),
                length(
                    intersect(de_genes, all_genes[!all_genes %in% list_genes])
                ),
                length(
                    intersect(all_genes[!all_genes %in% de_genes], list_genes)
                ),
                length(all_genes[!all_genes %in% c(de_genes, list_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    List = lit_list
                )
        }
    ) %>%
    bind_rows() 
de_iso_list_enrichment_specific <- expand.grid(
    Contrast = sort(unique(tt_iso$Contrast)),
    List = sort(unique(names(asd_lists)))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            lit_list <- as.character(param[["List"]])
            all_genes <- tt_iso %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_gene <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            de_iso_to_gene <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                distinct(external_gene_id) %>%
                filter(! external_gene_id %in% de_gene) %>%
                pull(external_gene_id)
            list_genes <- asd_lists[[lit_list]]
            contingency <- matrix(c(
                length(
                    intersect(de_iso_to_gene, list_genes)
                ),
                length(
                    intersect(
                        de_iso_to_gene, all_genes[!all_genes %in% list_genes]
                    )
                ),
                length(
                    intersect(
                        all_genes[!all_genes %in% de_iso_to_gene], list_genes
                    )
                ),
                length(all_genes[!all_genes %in% c(de_iso_to_gene, list_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    List = lit_list
                )
        }
    ) %>%
    bind_rows() 

################################################################################
# Plots                                                                        #
################################################################################

gene_list_enr <- de_gene_list_enrichment %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Gene")

isoform_list_enr <- de_iso_list_enrichment %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Isoform")

gene_list_enr_spec <- de_gene_list_enrichment_specific %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Gene, Specific")

isoform_list_enr_spec <- de_iso_list_enrichment_specific %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Isoform, Specific")

combined_list_enr <- bind_rows(
    de_gene_list_enrichment %>%
        mutate(Data = "Gene"),
    de_iso_list_enrichment %>%
        mutate(Data = "Isoform")
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

combined_list_enr_withSpec <- bind_rows(
    de_gene_list_enrichment %>%
        mutate(Data = "Gene, All"),
    de_iso_list_enrichment %>%
        mutate(Data = "Isoform, All"),
    de_gene_list_enrichment_specific %>%
        mutate(Data = "Gene, Specific"),
    de_iso_list_enrichment_specific %>%
        mutate(Data = "Isoform, Specific")
    
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

list_specOnly <- bind_rows(
    de_gene_list_enrichment_specific %>%
        mutate(Data = "Gene"),
    de_iso_list_enrichment_specific %>%
        mutate(Data = "Isoform")
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

list_enr_plot <- function(enr_data) {
    list_enr_plot_combined <- ggplot(
        data = enr_data %>%
            filter(
                ! List %in% c(
                    "VulnerableASD", "ASDSanders65",
                    "EichlerDNM_LGD_ASD_BD_SCZ", "EichlerDNM_LGD_MIS_ASD_BD_SCZ"
                )
            ),
        mapping = aes(
            x = Contrast, y = List, fill = -log10(adj.P.Val), label = Star
        )
    ) +
        facet_grid(Data ~ .) +
        geom_tile(
            colour = "lightgrey"
        ) +
        geom_text(fontface = "bold", vjust = 0.75, size = 10) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        # scale_fill_gradient2(
        #     low = "blue", mid = "white", high = "red", midpoint = 1
        # ) 
        scale_fill_gradient(low = "white", high = "red") +
        guides(
            fill = guide_colourbar(
                title = expression("-log"["10"]*"Adj. P-Value"),
                title.position = "top"
            )
        ) +
        theme_bw() +
        theme(
            text = element_text(size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.direction = "horizontal",
            legend.key.width = unit(3, units = "lines")
        )
}

enr_plots <- list(
    gene_only = list_enr_plot(enr_data = gene_list_enr),
    isoform_only = list_enr_plot(enr_data = isoform_list_enr),
    gene_spec = list_enr_plot(enr_data = gene_list_enr_spec),
    isoform_spec = list_enr_plot(enr_data = isoform_list_enr_spec),
    combined_all_only = list_enr_plot(enr_data = combined_list_enr),
    combined_all_andSpec = list_enr_plot(enr_data = combined_list_enr_withSpec),
    specOnly = list_enr_plot(enr_data = list_specOnly)
) %>%
    lapply(
        function(plt) {
            plt + theme(axis.text.x = element_text(angle = 30, hjust = 1))
        }
    )

saveRDS(
    enr_plots,
    "data/figures/DE_ASDListEnrichmentFisherCombinedPlot.rds"
)

pdf("data/figures/DE_ASDListEnrichmentsFisher.pdf", width = 16, height = 12)
invisible(lapply(enr_plots, print))
dev.off()

################################################################################
# Cell Types                                                                   #
################################################################################

cell_types <- list(
    readxl::read_xlsx(
        file.path(
            "data/source/CuratedLists",
            "ZhongNature2018_humanPreFrontalCortexCellTypes.xlsx"
        ),
        sheet = 1, skip = 4
    )[, c("cluster", "gene")] %>%
        filter(cluster == "NPCs"),
    readxl::read_xlsx(
        file.path(
            "data/source/CuratedLists",
            "WangScience2018CellTypesPsychencode.xlsx"
        ),
        sheet = 1, skip = 1
    ) %>%
        mutate(
            cluster = sapply(
                Cluster,
                function(clst) {
                    switch(
                        str_match(
                            clst,
                            "Ex|In|Astro|OPC|Oligo|Endo|Per|Microglia"
                        )[, 1],
                        "Ex" = "Excitatory neurons",
                        "In" = "Interneurons",
                        "Astro" = "Astrocytes",
                        "Oligo" = "Oligodendrocytes",
                        "OPC" = "OPC", 
                        "Endo" = "Endothelial",
                        "Per" = "Pericytes",
                        "Microglia" = "Microglia"
                    )
                }
            )
        ) %>%
        rename(gene = Gene) %>%
        dplyr::select(gene, cluster)
) %>%
    bind_rows() %>%
    distinct()

de_gene_ct_enrichment <- expand.grid(
    Contrast = sort(unique(tt_genes$Contrast)),
    CellType = sort(unique(cell_types$cluster))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            celltype <- as.character(param[["CellType"]])
            all_genes <- tt_genes %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_genes <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            ct_genes <- cell_types %>%
                filter(cluster == celltype) %>%
                pull(gene)
            contingency <- matrix(c(
                length(
                    intersect(de_genes, ct_genes)
                ),
                length(
                    intersect(de_genes, all_genes[!all_genes %in% ct_genes])
                ),
                length(
                    intersect(all_genes[!all_genes %in% de_genes], ct_genes)
                ),
                length(all_genes[!all_genes %in% c(de_genes, ct_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    CellType = celltype
                )
        }
    ) %>%
    bind_rows() 
de_iso_ct_enrichment <- expand.grid(
    Contrast = sort(unique(tt_iso$Contrast)),
    CellType = sort(unique(cell_types$cluster))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            celltype <- as.character(param[["CellType"]])
            all_genes <- tt_iso %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_iso_to_gene <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            ct_genes <- cell_types %>%
                filter(cluster == celltype) %>%
                pull(gene)
            contingency <- matrix(c(
                length(
                    intersect(de_iso_to_gene, ct_genes)
                ),
                length(
                    intersect(
                        de_iso_to_gene, all_genes[!all_genes %in% ct_genes]
                    )
                ),
                length(
                    intersect(
                        all_genes[!all_genes %in% de_iso_to_gene], ct_genes
                    )
                ),
                length(all_genes[!all_genes %in% c(de_iso_to_gene, ct_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    CellType = celltype
                )
        }
    ) %>%
    bind_rows()
de_gene_ct_enrichment_spec <- expand.grid(
    Contrast = sort(unique(tt_genes$Contrast)),
    CellType = sort(unique(cell_types$cluster))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            celltype <- as.character(param[["CellType"]])
            all_genes <- tt_genes %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_iso_gene_sym <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            de_genes <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                filter(! external_gene_id %in% de_iso_gene_sym) %>%
                pull(external_gene_id)
            ct_genes <- cell_types %>%
                filter(cluster == celltype) %>%
                pull(gene)
            contingency <- matrix(c(
                length(
                    intersect(de_genes, ct_genes)
                ),
                length(
                    intersect(de_genes, all_genes[!all_genes %in% ct_genes])
                ),
                length(
                    intersect(all_genes[!all_genes %in% de_genes], ct_genes)
                ),
                length(all_genes[!all_genes %in% c(de_genes, ct_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    CellType = celltype
                )
        }
    ) %>%
    bind_rows() 
de_iso_ct_enrichment_spec <- expand.grid(
    Contrast = sort(unique(tt_iso$Contrast)),
    CellType = sort(unique(cell_types$cluster))
) %>%
    apply(
        ., 1, function(param) {
            contrast <- as.character(param[["Contrast"]])
            celltype <- as.character(param[["CellType"]])
            all_genes <- tt_iso %>%
                distinct(external_gene_id) %>%
                pull(external_gene_id)
            de_gene <- tt_genes %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                pull(external_gene_id)
            de_iso_to_gene <- tt_iso %>%
                filter(Contrast == contrast) %>%
                filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                distinct(external_gene_id) %>%
                filter(! external_gene_id %in% de_gene) %>%
                pull(external_gene_id)
            ct_genes <- cell_types %>%
                filter(cluster == celltype) %>%
                pull(gene)
            contingency <- matrix(c(
                length(
                    intersect(de_iso_to_gene, ct_genes)
                ),
                length(
                    intersect(
                        de_iso_to_gene, all_genes[!all_genes %in% ct_genes]
                    )
                ),
                length(
                    intersect(
                        all_genes[!all_genes %in% de_iso_to_gene], ct_genes
                    )
                ),
                length(all_genes[!all_genes %in% c(de_iso_to_gene, ct_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(contingency) %>%
                broom::tidy() %>%
                rename(OR = estimate) %>%
                mutate(
                    Contrast = contrast,
                    CellType = celltype
                )
        }
    ) %>%
    bind_rows()

gene_ct_enr <- de_gene_ct_enrichment %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Gene")

isoform_ct_enr <- de_iso_ct_enrichment %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Isoform")

gene_ct_enr_spec <- de_gene_ct_enrichment_spec %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Gene, Specific")

isoform_ct_enr_spec <- de_iso_ct_enrichment_spec %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")) %>%
    mutate(Data = "Isoform, Specific")

combined_ct_enr <- bind_rows(
    de_gene_ct_enrichment %>%
        mutate(Data = "Gene"),
    de_iso_ct_enrichment %>%
        mutate(Data = "Isoform")
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

combined_ct_enr_withSpec <- bind_rows(
    de_gene_ct_enrichment %>%
        mutate(Data = "Gene, All"),
    de_iso_ct_enrichment %>%
        mutate(Data = "Isoform, All"),
    de_gene_ct_enrichment_spec %>%
        mutate(Data = "Gene, Specific"),
    de_iso_ct_enrichment_spec %>%
        mutate(Data = "Isoform, Specific")
    
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

ct_specOnly <- bind_rows(
    de_gene_ct_enrichment_spec %>%
        mutate(Data = "Gene"),
    de_iso_ct_enrichment_spec %>%
        mutate(Data = "Isoform")
) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.01, "*", ""))

ct_enr_plot <- function(enr_data) {
    list_enr_plot_combined <- ggplot(
        data = enr_data,
        mapping = aes(
            x = Contrast, y = CellType, fill = -log10(adj.P.Val), label = Star
        )
    ) +
        facet_grid(Data ~ .) +
        geom_tile(
            colour = "lightgrey"
        ) +
        geom_text(fontface = "bold", vjust = 0.75, size = 10) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        # scale_fill_gradient2(
        #     low = "blue", mid = "white", high = "red", midpoint = 1
        # ) +
        scale_fill_gradient(low = "white", high = "red") +
        guides(
            fill = guide_colourbar(
                title = expression("-log"["10"]*"Adj. P-Value"),
                title.position = "top"
            )
        ) +
        theme_bw() +
        theme(
            text = element_text(size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.direction = "horizontal",
            legend.key.width = unit(3, units = "lines")
        )
}

enr_plots <- list(
    gene_only = ct_enr_plot(enr_data = gene_ct_enr),
    isoform_only = ct_enr_plot(enr_data = isoform_ct_enr),
    gene_spec = ct_enr_plot(enr_data = gene_ct_enr_spec),
    isoform_spec = ct_enr_plot(enr_data = isoform_ct_enr_spec),
    combined_all_only = ct_enr_plot(enr_data = combined_ct_enr),
    combined_all_andSpec = ct_enr_plot(enr_data = combined_ct_enr_withSpec),
    specOnly = ct_enr_plot(enr_data = ct_specOnly)
) %>%
    lapply(
        function(plt) {
            plt + theme(axis.text.x = element_text(angle = 30, hjust = 1))
        }
    )

saveRDS(
    enr_plots,
    "data/figures/DE_CellTypeEnrichmentFisherCombinedPlot.rds"
)

pdf("data/figures/DE_CellTypeEnrichmentsFisher.pdf", width = 16, height = 12)
invisible(lapply(enr_plots, print))
dev.off()






ct_list_spec_combined <- ggplot(
    data = bind_rows(
        ct_specOnly %>%
            rename(List = CellType) %>%
            mutate(ListType = "Cell Type"),
        list_specOnly %>%
            mutate(ListType = "Curated Lists")
    ) %>%
        filter(
            ! List %in% c(
                "SFARISyndromicRisk12", "VulnerableASD", "ASDSanders65",
                "EichlerDNM_LGD_ASD_BD_SCZ", "EichlerDNM_LGD_MIS_ASD_BD_SCZ"
            )
        ),
    mapping = aes(
        y = Contrast, x = List, fill = -log10(adj.P.Val), label = Star
    )
) +
    facet_nested(. ~ ListType + Data, scales = "free_x") +
    geom_tile(
        colour = "lightgrey"
    ) +
    geom_text(fontface = "bold", vjust = 0.75, size = 10) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_gradient2(
    #     low = "blue", mid = "white", high = "red", midpoint = 1
    # ) +
    scale_fill_gradient(low = "white", high = "red") +
    guides(
        fill = guide_colourbar(
            title = expression("-log"["10"]*"Adj. P-Value"),
            title.position = "top"
        )
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "horizontal",
        legend.key.width = unit(3, units = "lines"),
        panel.spacing = unit(0, "lines")
    )
saveRDS(
    ct_list_spec_combined, 
    "data/figures/DE_FisherEnrich_CTandList_Combined.rds"
)

write_csv(
    bind_rows(
        ct_specOnly %>%
            rename(List = CellType) %>%
            mutate(ListType = "Cell Type"),
        list_specOnly %>%
            mutate(ListType = "Curated Lists")
    ) %>%
        filter(
            ! List %in% c(
                "SFARISyndromicRisk12", "VulnerableASD", "ASDSanders65",
                "EichlerDNM_LGD_ASD_BD_SCZ", "EichlerDNM_LGD_MIS_ASD_BD_SCZ"
            )
        ),
    "data/ListEnrichmentDE_Specifics.csv"
)
