library(tidyverse)
library(WGCNA)
library(ggdendro)
library(dendextend)
library(lmerTest)
library(cowplot)
source("scripts/utility/plotting.R")

metadata <- read_tsv("data/metadata.tsv")
feature_anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]
sv_gene <- readRDS("data/genes/sva_results/13.rds")$sv
sv_iso <- readRDS("data/isoforms/sva_results/16.rds")$sv
colnames(sv_gene) <- paste0("Gene", colnames(sv_gene))
colnames(sv_iso) <- paste0("Isoform", colnames(sv_iso))
metadata_sv <- inner_join(
    inner_join(
        metadata, 
        mutate(as.data.frame(sv_gene), Sample = rownames(sv_gene)), 
        by = "Sample"
    ),
    mutate(as.data.frame(sv_iso), Sample = rownames(sv_iso)), 
    by = "Sample"
)

gene_net <- readRDS("data/genes/Networks/Network_DS2_MM20.rds")
iso_net <- readRDS("data/isoforms/Networks/Network_DS2_MM20.rds")
module_assigns <- list(
    gene = read_tsv(
        "data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv"
    ) %>%
        left_join(
            feature_anno %>% 
                dplyr::select(-contains("transcript")) %>%
                distinct(),
            by = c("Feature" = "ensembl_gene_id")
        ),
    isoform = read_tsv(
        "data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv"
    ) %>%
        left_join(
            feature_anno, by = c("Feature" = "ensembl_transcript_id")
        )
)

module_maps <- bind_rows(
    module_assigns$gene %>%
        distinct(module_label, module_colour) %>%
        mutate(ME = paste0("gME", module_label)),
    module_assigns$isoform %>%
        distinct(module_label, module_colour) %>%
        mutate(ME = paste0("iME", module_label))
)

# ME Clustering

colnames(gene_net$MEs) <- paste0("g", colnames(gene_net$MEs))
colnames(iso_net$MEs) <- paste0("i", colnames(iso_net$MEs))

MEs <- inner_join(
    mutate(gene_net$MEs, Sample = rownames(gene_net$MEs)),
    mutate(iso_net$MEs, Sample = rownames(iso_net$MEs)),
    by = "Sample"
)

clustMEs <- hclust(dist(t(
    dplyr::select(MEs, -Sample, -gME0, -iME0)
)), method = "average")
dendMEs <- as.dendrogram(clustMEs)
dendDataMEs <- dendro_data(dendMEs)
segmentMEs <- segment(dendDataMEs)
colnames(segmentMEs) <- c("y", "x", "yend", "xend")
module_pos <- with(
    dendDataMEs$labels,
    data.frame(y_center = x, module = as.character(label), height = 1)
) %>%
    left_join(module_maps, by = c("module" = "ME")) %>%
    mutate(
        Data = sapply(
            module, function(m) {
                ifelse(substr(m, 1, 1) == "g", "Gene", "Isoform")
            }
        )
    )

dend_plot <- ggplot(
    data = segmentMEs
) +
    geom_segment(
        mapping = aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_point(
        data = module_pos,
        inherit.aes = FALSE,
        mapping = aes(y = y_center, x = 0, fill = module_colour, colour = Data),
        shape = 21, show.legend = FALSE, size = 5, stroke = 1
    ) +
    geom_tile(
        data = module_pos, 
        mapping = aes(y = y_center, x = 0),
        show.legend = FALSE, fill = NA, colour = NA
    ) +
    geom_text(
        data = module_pos,
        mapping = aes(
            y = y_center, x = -0.1, 
            label = gsub("ME", "M", module), colour = Data
        ),
        show.legend = FALSE, fontface = "bold", hjust = 0,
        size = 3
    ) +
    scale_y_discrete(
        breaks = module_pos$y_center,
        expand = c(0, 0)
    ) +
    scale_x_reverse(
        expand = expand_scale(mult = c(0, 0.05))
    ) +
    scale_fill_manual(
        values = setNames(nm = module_maps$module_colour)
    ) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )
saveRDS(
    dend_plot,
    "data/figures/GeneIsoformMEDendrogram.rds"
)

# Module Trait Association

module_assoc <- lapply(
    colnames(MEs)[grepl("ME", colnames(MEs))],
    function(me) {
        this_formula <- paste(
            me, 
            "~ 0 + Period + Regioncode + Sex + Ethnicity + Site"
        )
        if (substr(me, 1, 1) == "g") {
            this_formula <- paste(
                this_formula,
                "+",
                paste(
                    colnames(metadata_sv)[grepl(
                        "GeneSV", colnames(metadata_sv)
                    )],
                    collapse = " + "
                )
            )
        } else {
            this_formula <- paste(
                this_formula,
                "+",
                paste(
                    colnames(metadata_sv)[grepl(
                        "IsoformSV", colnames(metadata_sv)
                    )],
                    collapse = " + "
                )
            )
        }
        summary(lmer(
            formula = as.formula(
                paste(
                    this_formula, "+ (1|Braincode)" 
                )
            ),
            data = inner_join(MEs, metadata_sv, by = "Sample") %>% 
                mutate(Period = as.factor(Period))
        ))$coefficients %>%
            as.data.frame() %>%
            mutate(Factor = gsub("\\(Intercept\\)", "Period2", rownames(.))) %>%
            filter(str_detect(Factor, "Period")) %>%
            mutate(module = me)
    }
) %>%
    bind_rows() %>%
    bind_rows(
        lapply(
            colnames(MEs)[grepl("ME", colnames(MEs))],
            function(me) {
                this_formula <- paste(
                    me, 
                    "~ 0 + Prenatal + Regioncode + Sex + Ethnicity + Site"
                )
                if (substr(me, 1, 1) == "g") {
                    this_formula <- paste(
                        this_formula,
                        "+",
                        paste(
                            colnames(metadata_sv)[grepl(
                                "GeneSV", colnames(metadata_sv)
                            )],
                            collapse = " + "
                        )
                    )
                } else {
                    this_formula <- paste(
                        this_formula,
                        "+",
                        paste(
                            colnames(metadata_sv)[grepl(
                                "IsoformSV", colnames(metadata_sv)
                            )],
                            collapse = " + "
                        )
                    )
                }
                summary(lmer(
                    formula = as.formula(
                        paste(
                            this_formula, "+ (1|Braincode)" 
                        )
                    ),
                    data = inner_join(MEs, metadata_sv, by = "Sample") %>% 
                        mutate(Period = as.factor(Period))
                ))$coefficients %>%
                    as.data.frame() %>%
                    mutate(
                        Factor = gsub("\\(Intercept\\)", "Period2", rownames(.))
                    ) %>%
                    filter(str_detect(Factor, "natal")) %>%
                    mutate(
                        Factor = ifelse(
                            str_detect(Factor, "TRUE"), "Prenatal", "Postnatal"
                        )
                    ) %>%
                    mutate(module = me)
            }
        ) %>%
            bind_rows()
    ) %>%
    mutate(FDR = p.adjust(`Pr(>|t|)`, method = "fdr")) %>%
    left_join(module_pos, by = "module") %>%
    mutate(
        Factor = factor(
            Factor, 
            levels = c("Postnatal", "Prenatal", paste0("Period", seq(13, 2)))
        )
    )

module_assoc_plot <- ggplot(
    data = module_assoc %>%
        mutate(Star = ifelse(FDR <= 0.05, "*", "")) %>%
        mutate(
            Factor = factor(
                Factor,
                levels = c(
                    paste0("Period", seq(2, 13)), "Prenatal", "Postnatal"
                )
            )
        ) %>%
        mutate(Strip = "Module Association"),
    mapping = aes(
        x = Factor, 
        y = y_center,
        fill = Estimate, label = Star, colour = Data
    )
) +
    facet_grid(. ~ Strip) +
    geom_tile(colour = "lightgrey") +
    geom_text(colour = "black", size = 5, fontface = "bold") +
    guides(
        colour = FALSE, 
        fill = guide_colorbar(
            title = "Linear Regression Beta", title.position = "top"
        )
    ) +
    scale_x_discrete(
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = module_pos$y_center,
        expand = c(0, 0)
    ) +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_blank(),
        # legend.direction = "horizontal",
        legend.position = "top",
        legend.key.width = unit(5, "lines")
    )
saveRDS(
    module_assoc_plot,
    "data/figures/ModuleAssociationTiles.rds"
)

# Module-List Enrichments

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
cell_types <- list(
    readxl::read_xlsx(
        file.path(
            "data/source/CuratedLists",
            "ZhongNature2018_humanPreFrontalCortexCellTypes.xlsx"
        ),
        sheet = 1, skip = 4
    )[, c("cluster", "gene")],
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
cell_types_lists <- lapply(
    setNames(
        nm = sort(unique(cell_types$cluster))
    ),
    function(ct) {
        pull(filter(cell_types, cluster == ct), gene)
    }
)
all_lists <- c(asd_lists, cell_types_lists)

list_fisher_enrichment <- bind_rows(
    lapply(
        names(all_lists),
        function(l) {
            bind_rows(
                mutate(
                    distinct(
                        module_assigns$gene, module_label, module_colour
                    ),
                    Data = "Gene"
                ),
                mutate(
                    distinct(
                        module_assigns$isoform, module_label, module_colour
                    ),
                    Data = "Isoform"
                )
            ) %>%
                mutate(List = l)
        }
    )
) %>%
    mutate(
        List_Class = sapply(
            List,
            function(l) {
                if (l %in% names(asd_lists)) {
                    return("Curated Lists")
                } else if (l %in% names(cell_types_lists)) {
                    return("Cell Type")
                } else {
                    return("")
                }
            }
        )
    ) %>%
    apply(
        ., 1, function(param) {
            ml <- as.numeric(param[["module_label"]])
            mc <- as.character(param[["module_colour"]])
            data_type <- as.character(param[["Data"]])
            enr_list <- as.character(param[["List"]])
            list_type <- as.character(param[["List_Class"]])
            all_genes <- unique(c(
                module_assigns$gene$external_gene_id,
                module_assigns$isoform$external_gene_id
            ))
            if (data_type == "Gene") {
                module_genes <- module_assigns$gene
            } else {
                module_genes <- module_assigns$isoform
            }
            module_genes <- module_genes %>%
                filter(module_label == ml) %>%
                pull(external_gene_id)
            list_genes <- all_lists[[enr_list]]
            contingency <- matrix(c(
                length(intersect(
                    module_genes, list_genes
                )),
                length(intersect(
                    module_genes, all_genes[! all_genes %in% list_genes]
                )),
                length(intersect(
                    all_genes[! all_genes %in% module_genes], list_genes
                )),
                length(all_genes[! all_genes %in% c(module_genes, list_genes)])
            ), byrow = TRUE, ncol = 2, nrow = 2)
            fisher.test(x = contingency) %>%
                broom::tidy() %>%
                rename(OddsRatio = estimate) %>%
                mutate(
                    module_label = ml,
                    module_colour = mc,
                    Data = data_type,
                    List = enr_list,
                    List_Class = list_type
                )
        }
    ) %>%
    bind_rows() %>%
    left_join(
        module_pos,
        by = c("module_label", "module_colour", "Data")
    )

list_enrichment_plot <- ggplot(
    data = list_fisher_enrichment %>%
        filter(module_label != 0) %>%
        filter(
            ! List %in% c(
                "VulnerableASD", "ASDSanders65",
                "EichlerDNM_LGD_ASD_BD_SCZ", "EichlerDNM_LGD_MIS_ASD_BD_SCZ",
                "Pericytes"
            )
        ) %>%
        mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")),
    mapping = aes(
        y = y_center, x = List, 
        # fill = OddsRatio, 
        fill = -log10(adj.P.Val),
        label = Star
    )
) +
    # facet_nested(Data + List_Class ~ ., scales = "free_y") +
    facet_grid(. ~ List_Class, scales = "free_x") +
    geom_tile(colour = "lightgrey") +
    geom_text(colour = "black", size = 5, fontface = "bold") +
    guides(
        colour = FALSE, 
        fill = guide_colorbar(
            title = expression("-log"["10"]*"Adj. P-Value"),
            title.position = "top"
        )
    ) +
    scale_x_discrete(
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = module_pos$y_center,
        expand = c(0, 0)
    ) +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0, 
        # trans = "log1p", breaks = c(0, 1, 10, 50)
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.spacing = unit(0, units = "lines"),
        strip.background = element_rect(colour = "black"),
        axis.title = element_blank(),
        legend.key.height = unit(2, units = "lines"),
        # legend.direction = "horizontal",
        legend.position = "top",
        legend.key.width = unit(5, "lines")
    )
saveRDS(
    list_enrichment_plot,
    "data/figures/ModuleListEnrichmentFisher.rds"
)

# Gene-isoform module pair overlaps

analogous_modules <- data.frame(
    Gene = seq(1, 8),
    Isoform = c(1, 2, 4, 35, 11, 9, 23, 41)
)

analogous_modules_overlaps <- apply(
    analogous_modules, 1, function(x) {
        gml <- as.numeric(x[["Gene"]])
        iml <- as.numeric(x[["Isoform"]])
        g_genes <- module_assigns$gene %>%
            filter(module_label == gml) %>%
            pull(Feature)
        i_genes <- module_assigns$isoform %>%
            filter(module_label == iml) %>%
            pull(ensembl_gene_id)
        data.frame(
            gene_module_label = gml,
            isoform_module_label = iml,
            gene_module_size = length(g_genes),
            isoform_module_size = length(i_genes),
            gene_unique = length(g_genes[!g_genes %in% i_genes]),
            isoform_unique = length(i_genes[!i_genes %in% g_genes]),
            overlap = length(intersect(g_genes, i_genes)),
            jaccard_index = 1 - (
                length(intersect(g_genes, i_genes)) /
                    length(union(g_genes, i_genes))
            )
        )
    }
) %>%
    bind_rows()

jaccard_analysis <- apply(
    analogous_modules, 1, function(x) {
        gml <- as.numeric(x[["Gene"]])
        iml <- as.numeric(x[["Isoform"]])
        all_genes <- unique(c(
            module_assigns$gene$Feature, 
            module_assigns$isoform$ensembl_gene_id
        ))
        g_genes <- module_assigns$gene %>%
            filter(module_label == gml) %>%
            pull(Feature)
        i_genes <- module_assigns$isoform %>%
            filter(module_label == iml) %>%
            pull(ensembl_gene_id)
        g_mask <- setNames(
            vector(mode = "numeric", length = length(all_genes)),
            all_genes
        )
        g_mask[which(names(g_mask) %in% g_genes)] <- 1
        i_mask <- setNames(
            vector(mode = "numeric", length = length(all_genes)),
            all_genes
        )
        i_mask[which(names(i_mask) %in% g_genes)] <- 1
        jaccard::jaccard.test(
            g_mask, i_mask, method = "bootstrap", B = 10000
        )
    }
)
jaccard_compile <- analogous_modules %>%
    bind_cols(
        bind_rows(lapply(jaccard_analysis, function(j) {
            data.frame(
                jaccard = j$statistics,
                p_value = j$pvalue,
                expect = j$expectation
            )
        }))
    )

write_csv(module_assoc, "data/WGCNA/ModulePeriodAssociation.csv")
write_csv(
    list_fisher_enrichment %>%
        filter(module_label != 0) %>%
        filter(
            ! List %in% c(
                "VulnerableASD", "ASDSanders65",
                "EichlerDNM_LGD_ASD_BD_SCZ", "EichlerDNM_LGD_MIS_ASD_BD_SCZ",
                "Pericytes"
            )
        ) %>%
        mutate(adj.P.Val = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(Star = ifelse(adj.P.Val <= 0.01, "*", "")), 
    "data/WGCNA/ModuleListEnrichment.csv"
)
