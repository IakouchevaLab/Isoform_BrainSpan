library(tidyverse)
library(igraph)
library(ggrepel)
source("scripts/utility/graph.R")

# # Vertex tables
# annotations <- read.table(
#     "data/source/annotation.transcript.ensg75.txt", 
#     header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE
# )
# variants <- readxl::read_xlsx(
#     "data/SupplementaryTables/Supplementary Table 6.xlsx",
#     sheet = 2
# )
# control_targets <- c(
#     filter(variants, `Affected status` == 1)$`Ensembl Gene ID`,
#     filter(variants, `Affected status` == 1)$`Ensembl Transcript ID`
# )
# asd_targets <- c(
#     filter(variants, `Affected status` == 2)$`Ensembl Gene ID`,
#     filter(variants, `Affected status` == 2)$`Ensembl Transcript ID`
# )
# 
# SFARI1 <- c(
#     "ADNP", "ARID1B", "DYRK1A", "CHD2", "ASH1L", "ASXL3", "CHD8", "ANK2", 
#     "CUL3", "DSCAM", "MYT1L", "KMT2A", "KMT5B", "GRIN2B", "POGZ", "KATNAL2", 
#     "NAA15", "PTEN", "TRIP12", "RELN", "TBR1", "SCN2A", "SYNGAP1", "SETD5", 
#     "SHANK3"
# )
# SANDERS <- c(
#     "CHD8", "SCN2A", "ARID1B", "SYNGAP1", "DYRK1A", "CHD2", "ANK2", "KDM5B", 
#     "ADNP", "POGZ", "KMT5B", "TBR1", "GRIN2B", "DSCAM", "KMT2C", "TCF7L2", 
#     "TRIP12", "ASH1L", "CUL3", "KATNAL2", "GIGYF1", "TNRC6B", "WAC", "NCKAP1", 
#     "RANBP17", "KDM6B", "ILF2", "SPAST", "FOXP1", "AKAP9", "CMPK2", "DDX3X", 
#     "WDFY3", "PHF2", "BCL11A", "KMT2E", "CACNA2D3", "NRXN1", "SHANK2", "PTEN", 
#     "SHANK3", "SETD5", "DNMT3A", "MYT1L", "RAPGEF4", "PRKAR1B", "MFRP", 
#     "GABRB3", "P2RX5", "ETFB", "CTTNBP2", "INTS6", "USP45", "ERBIN", "TMEM39B", 
#     "TSPAN4", "MLANA", "SMURF1", "C16orf13", "BTRC", "CCSER1", "FAM98C", 
#     "SLC6A1", "ZNF559", "CAPN12", "GRIA1", "PCM1", "MYO5A", "UIMC1"
# )
# SATTERSTROM <- readxl::read_xlsx(
#     "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
#     sheet = "SatterstromASD"
# )[[1]]
# vertices <- bind_rows(
#     lapply(
#         setNames(
#             nm = readxl::excel_sheets(
#                 "data/SupplementaryTables/Supplementary Table 7.xlsx"
#             )[-1]
#         ),
#         function(sheet) {
#             readxl::read_xlsx(
#                 "data/SupplementaryTables/Supplementary Table 7.xlsx",
#                 sheet = sheet
#             )
#         }
#     ) %>%
#         bind_rows() %>%
#         mutate(name = `Ensembl Gene ID`) %>%
#         mutate(Data = "Gene"),
#     lapply(
#         setNames(
#             nm = readxl::excel_sheets(
#                 "data/SupplementaryTables/Supplementary Table 8.xlsx"
#             )[-1]
#         ),
#         function(sheet) {
#             readxl::read_xlsx(
#                 "data/SupplementaryTables/Supplementary Table 8.xlsx",
#                 sheet = sheet
#             )
#         }
#     ) %>%
#         bind_rows() %>%
#         mutate(name = `Ensembl Transcript ID`) %>%
#         mutate(Data = "Isoform")
# ) %>%
#     dplyr::select(name, everything()) %>%
#     mutate(
#         Target = case_when(
#             name %in% control_targets & name %in% asd_targets ~ "Both",
#             ! name %in% control_targets & name %in% asd_targets ~ "ASD",
#             name %in% control_targets & ! name %in% asd_targets ~ "Control",
#             ! name %in% control_targets & ! name %in% asd_targets ~ "Neither"
#         )
#     ) %>%
#     mutate(
#         SATTERSTROM = `Gene Symbol` %in% SATTERSTROM,
#         SANDERS = `Gene Symbol` %in% SANDERS,
#         SFARI1 = `Gene Symbol` %in% SFARI1
#     ) %>%
#     mutate(
#         ASD = mapply(
#             function(satterstrom, sanders, sfari1) {
#                 return(any(c(satterstrom, sanders, sfari1)))
#             }, SATTERSTROM, SANDERS, SFARI1
#         )
#     )
# gM1_vertices <- filter(vertices, Data == "Gene" & `Module Label` == 1)
# iM1_vertices <- filter(vertices, Data == "Isoform" & `Module Label` == 1)
# iM30_vertices <- filter(vertices, Data == "Isoform" & `Module Label` == 30)
# 
# # Edges
# gene_ppi <- read.table(
#     "data/source/CuratedLists/allVidal.lit17.mvp.hprdInvivo.hintBinary.txt",
#     header = TRUE, stringsAsFactors = FALSE
# )[, 1:2] %>%
#     setNames(nm = c("head", "tail"))
# gexpr <- readRDS("data/RegressGeneCounts.rds")
# iexpr <- readRDS("data/RegressIsoformCounts.rds")
# 
# gM1_edges <- gene_ppi %>%
#     filter(head %in% gM1_vertices$`Ensembl Gene ID`) %>%
#     filter(tail %in% gM1_vertices$`Ensembl Gene ID`) %>%
#     distinct() %>%
#     filter(head != tail) %>%
#     mutate(
#         weight = mapply(function(h, t) cor(gexpr[h, ], gexpr[t, ]), head, tail)
#     )
# iM1_edges <- gene_ppi %>%
#     filter(head %in% iM1_vertices$`Ensembl Gene ID`) %>%
#     filter(tail %in% iM1_vertices$`Ensembl Gene ID`) %>%
#     filter(head != tail) %>%
#     distinct() %>%
#     left_join(
#         dplyr::select(iM1_vertices, `Ensembl Gene ID`, `Ensembl Transcript ID`),
#         by = c("head" = "Ensembl Gene ID")
#     ) %>%
#     dplyr::select(-head) %>%
#     rename(head = `Ensembl Transcript ID`) %>%
#     left_join(
#         dplyr::select(iM1_vertices, `Ensembl Gene ID`, `Ensembl Transcript ID`),
#         by = c("tail" = "Ensembl Gene ID")
#     ) %>%
#     dplyr::select(-tail) %>%
#     rename(tail = `Ensembl Transcript ID`) %>%
#     filter(head != tail) %>%
#     mutate(
#         weight = mapply(function(h, t) cor(iexpr[h, ], iexpr[t, ]), head, tail)
#     )
# iM30_edges <- gene_ppi %>%
#     filter(head %in% iM30_vertices$`Ensembl Gene ID`) %>%
#     filter(tail %in% iM30_vertices$`Ensembl Gene ID`) %>%
#     filter(head != tail) %>%
#     distinct() %>%
#     left_join(
#         dplyr::select(iM30_vertices, `Ensembl Gene ID`, `Ensembl Transcript ID`),
#         by = c("head" = "Ensembl Gene ID")
#     ) %>%
#     dplyr::select(-head) %>%
#     rename(head = `Ensembl Transcript ID`) %>%
#     left_join(
#         dplyr::select(iM30_vertices, `Ensembl Gene ID`, `Ensembl Transcript ID`),
#         by = c("tail" = "Ensembl Gene ID")
#     ) %>%
#     dplyr::select(-tail) %>%
#     rename(tail = `Ensembl Transcript ID`) %>%
#     filter(head != tail) %>%
#     mutate(
#         weight = mapply(function(h, t) cor(iexpr[h, ], iexpr[t, ]), head, tail)
#     )
# 
# # Create graphs
# gM1_graph <- graph_from_data_frame(
#     d = gM1_edges, vertices = gM1_vertices, directed = FALSE
# )
# iM1_graph <- graph_from_data_frame(
#     d = iM1_edges, vertices = iM1_vertices, directed = FALSE
# )
# iM30_graph <- graph_from_data_frame(
#     d = iM30_edges, vertices = iM30_vertices, directed = FALSE
# )
# 
# # Write graphs to files
# writexl::write_xlsx(graph_to_tables(gM1_graph), "data/ppi/gM1_FullGraph.xlsx")
# writexl::write_xlsx(graph_to_tables(iM1_graph), "data/ppi/iM1_FullGraph.xlsx")
# writexl::write_xlsx(graph_to_tables(iM30_graph), "data/ppi/iM30_FullGraph.xlsx")
# saveRDS(gM1_graph, "data/ppi/gM1_FullGraph.rds")
# saveRDS(iM1_graph, "data/ppi/iM1_FullGraph.rds")
# saveRDS(iM30_graph, "data/ppi/iM30_FullGraph.rds")

gM1_graph <- readRDS("data/ppi/gM1_FullGraph.rds")
iM1_graph <- readRDS("data/ppi/iM1_FullGraph.rds")
iM30_graph <- readRDS("data/ppi/iM30_FullGraph.rds")

# ggplot() +
#     geom_density(
#         data = data.frame(
#             module = "gM1", weight = E(gM1_graph)$weight
#         ),
#         mapping = aes(
#             x = weight, fill = module, group = module
#         ),
#         alpha = 0.5
#     ) +
#     geom_density(
#         data = data.frame(
#             module = "iM1", weight = E(iM1_graph)$weight
#         ),
#         mapping = aes(
#             x = weight, fill = module, group = module
#         ),
#         alpha = 0.5
#     ) +
#     geom_density(
#         data = data.frame(
#             module = "iM30", weight = E(iM30_graph)$weight
#         ) %>%
#             filter(weight != 1), 
#         mapping = aes(
#             x = weight, fill = module, group = module
#         ),
#         alpha = 0.5
#     )

# Targets connected to at least one ASD
gM1_direct_q20 <- gM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= quantile(abs(E(.)$weight), 0.8))]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )
gM1_direct_q10 <- gM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= quantile(abs(E(.)$weight), 0.9))]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )
# plot_graph_byData(gM1_direct, "white")
iM1_direct_q20 <- iM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= quantile(abs(E(.)$weight), 0.8))]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )
plot_graph_byData(iM1_direct_q20, "white")
iM1_direct_q10 <- iM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= quantile(abs(E(.)$weight), 0.9))]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )

iM30_graph %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD | V(net)$Target %in% c("ASD", "Both"))) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    plot_graph_byData(colour = "white", plot_singletons = TRUE)

iM1_direct_PCC07 <- iM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= 0.7)]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )
gM1_direct_PCC07 <- gM1_graph %>%
    delete_edges(
        E(.)[which(abs(E(.)$weight) <= 0.7)]
    ) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD)) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[which(!sapply(
            V(.),
            function(v) {
                if (V(.)[[v]]$ASD) {
                    return(v)
                } else if (V(.)[[v]]$Target %in% c("ASD", "Both")) {
                    return(any(vertex_attr(., "ASD", neighbors(., v))))
                } else {
                    return(FALSE)
                }
            }
        ))]
    )

pdf("data/ppi/DirectConnect_Target_ASD_PCCcutoffs.pdf", width = 16, height = 12) 
print(plot_graph_byData(gM1_direct_q20, "white") + labs(
    title = "gM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "Top 20% PCC"
))
print(plot_graph_byData(gM1_direct_q10, "white") + labs(
    title = "gM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "Top 10% PCC"
))
print(plot_graph_byData(gM1_direct_PCC07, "white") + labs(
    title = "gM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "abs(PCC) >= 0.7"
) + 
    scale_size_continuous(breaks  = seq(0.7, 1, 0.05), range = c(0.01, 3)) +
    scale_alpha_continuous(breaks  = seq(0.7, 1, 0.05))
)
print(plot_graph_byData(iM1_direct_q20, "white") + labs(
    title = "iM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "Top 20% PCC"
))
print(plot_graph_byData(iM1_direct_q10, "white") + labs(
    title = "iM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "Top 10% PCC"
))
print(plot_graph_byData(iM1_direct_PCC07, "white") + labs(
    title = "iM1: Direct connections betw. ASD and ASD LoF Targets",
    subtitle = "abs(PCC) >= 0.7"
) + 
    scale_size_continuous(breaks  = seq(0.7, 1, 0.05), range = c(0.01, 3)) +
    scale_alpha_continuous(breaks  = seq(0.7, 1, 0.05))
)
# print(
#     plot_graph_byData_colorEdges(
#         gM1_direct_q20,
#         "white",
#         gene_isoform_overlapping_edges(gM1_direct_q20, iM1_graph)$gene_edges
#     ) +
#         labs(
#             title = "gM1: Direct connections betw. ASD and ASD LoF Targets",
#             subtitle = "Top 20% PCC, overlap with iM1 Top 10% PCC"
#         )
# )
# print(
#     plot_graph_byData_colorEdges(
#         iM1_direct_q10,
#         "white",
#         gene_isoform_overlapping_edges(gM1_graph, iM1_direct_q10)$isoform_edges
#     ) +
#         labs(
#             title = "iM1: Direct connections betw. ASD and ASD LoF Targets",
#             subtitle = "Top 10% PCC, overlap with gM1 Top 20% PCC"
#         )
# )
iM30_graph %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD | V(net)$Target %in% c("ASD", "Both"))) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact() %>%
    disjoint_union() %>% {
        plot_graph_byData(., colour = "white", plot_singletons = TRUE) +
            labs(
                title = "iM30: ASD and ASD LoF Targets"
            )
    } %>%
    print()
dev.off()    

# Finalized plots
# Grouped iM1
iM1_direct_q10_noSingle <- delete_vertices(
    iM1_direct_q10,
    V(iM1_direct_q10)[degree(iM1_direct_q10) == 0]
)
iM1_direct_q10_noSingle_groupByGene <- iM1_direct_q10_noSingle %>% {
    for (gene in unique(vertex_attr(., name = "Ensembl Gene ID"))) {
        nodes <- V(.)[V(.)$`Ensembl Gene ID` == gene]
        if (length(nodes) > 1) {
            . <- add_edges(
                ., 
                combn(V(.)[V(.)$`Ensembl Gene ID` == gene], 2), 
                attr = list(weight = 1.1)
            )
        }
    }
    return(.)
}
iM1_vertex_layout <- iM1_direct_q10_noSingle_groupByGene %>%
    layout_with_dh() %>%
    as.data.frame(stringsAsFactors = FALSE) %>% 
    setNames(nm = c("x", "y")) %>%
    mutate(name = names(V(iM1_direct_q10_noSingle_groupByGene))) %>%
    inner_join(
        as.data.frame(
            vertex_attr(iM1_direct_q10_noSingle_groupByGene), 
            stringsAsFactors = FALSE
        ), by = "name"
    ) %>%
    mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
iM1_edge_layout <- as_edgelist(iM1_direct_q10_noSingle) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(nm = c("head", "tail")) %>% 
    bind_cols(
        as.data.frame(edge_attr(iM1_direct_q10_noSingle))
    ) %>% {
        heads <- as.character(.$head)
        tails <- as.character(.$tail)
        head_coords <- iM1_vertex_layout[match(heads, iM1_vertex_layout$name), ] %>%
            dplyr::select(x, y) %>%
            setNames(nm = paste(names(.), "head", sep = "_"))
        tail_coords <- iM1_vertex_layout[match(tails, iM1_vertex_layout$name), ] %>%
            dplyr::select(x, y) %>%
            setNames(nm = paste(names(.), "tail", sep = "_"))
        bind_cols(., head_coords, tail_coords)
    }
datasets <- data.frame(
    Dataset = c("SATTERSTROM", "SANDERS", "SFARI1"),
    Fill = RColorBrewer::brewer.pal(n = 3, name = "Set1"),
    TextSize = c(3, 2, 1),
    stringsAsFactors = FALSE
)
iM1_plot <- ggplot() +
    geom_segment(
        data = iM1_edge_layout,
        mapping = aes(
            x = x_head, y = y_head, xend = x_tail, yend = y_tail,
            size = weight, alpha = weight
        ), 
        colour = "darkgrey", lineend = "round", linejoin = "round"
    ) +
    geom_point(
        data = datasets,
        mapping = aes(x = 0, y = 0, fill = Dataset), 
        size = 3, shape = 21, alpha = 0
    ) +
    scale_fill_manual(values = setNames(datasets$Fill, datasets$Dataset)) +
    geom_point(
        data = iM1_vertex_layout %>%
            filter(SATTERSTROM),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SATTERSTROM"], 
        size = datasets$TextSize[datasets$Dataset == "SATTERSTROM"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = iM1_vertex_layout %>%
            filter(SANDERS),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SANDERS"], 
        size = datasets$TextSize[datasets$Dataset == "SANDERS"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = iM1_vertex_layout %>%
            filter(SFARI1),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SFARI1"], 
        size = datasets$TextSize[datasets$Dataset == "SFARI1"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = iM1_vertex_layout %>%
            filter(!ASD),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = "white", 
        stroke = min(1 / 10),
        size = datasets$TextSize[datasets$Dataset == "SFARI1"]
    ) +
    geom_text_repel(
        data = iM1_vertex_layout,
        mapping = aes(
            x = x, y = y, 
            label = ifelse(
                Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
            )
        ), 
        size = 1.75, min.segment.length = 0.1,
        segment.size = min(1 / 10),
    ) +
    scale_size_continuous(
        breaks  = seq(0.7, 1, 0.05), range = c(0.1, 0.5)
    ) +
    scale_alpha_continuous(breaks  = seq(0.7, 1, 0.05)) +
    scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
    guides(
        size = guide_legend(title = "PCC"),
        alpha = guide_legend(title = "PCC"),
        shape = guide_legend(override.aes = list(fill = NA)),
        fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggnetwork::theme_blank() +
    theme(
        text = element_text(size = 8),
        legend.position = "none"
    )
ggsave(
    filename = "data/ppi/DirectConnect_iM1.pdf",
    plot = iM1_plot,
    device = "pdf", width = 13, height = 9, units = "cm", 
    useDingbats = FALSE
)





gM1_direct_q10_noSingle <- delete_vertices(
    gM1_direct_q10,
    V(gM1_direct_q10)[degree(gM1_direct_q10) == 0]
)
gM1_vertex_layout <- gM1_direct_q10_noSingle %>%
    layout_with_dh() %>%
    as.data.frame(stringsAsFactors = FALSE) %>% 
    setNames(nm = c("x", "y")) %>%
    mutate(name = names(V(gM1_direct_q10_noSingle))) %>%
    inner_join(
        as.data.frame(
            vertex_attr(gM1_direct_q10_noSingle), 
            stringsAsFactors = FALSE
        ), by = "name"
    ) %>%
    mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
gM1_edge_layout <- as_edgelist(gM1_direct_q10_noSingle) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(nm = c("head", "tail")) %>% 
    bind_cols(
        as.data.frame(edge_attr(gM1_direct_q10_noSingle))
    ) %>% {
        heads <- as.character(.$head)
        tails <- as.character(.$tail)
        head_coords <- gM1_vertex_layout[match(heads, gM1_vertex_layout$name), ] %>%
            dplyr::select(x, y) %>%
            setNames(nm = paste(names(.), "head", sep = "_"))
        tail_coords <- gM1_vertex_layout[match(tails, gM1_vertex_layout$name), ] %>%
            dplyr::select(x, y) %>%
            setNames(nm = paste(names(.), "tail", sep = "_"))
        bind_cols(., head_coords, tail_coords)
    }
datasets <- data.frame(
    Dataset = c("SATTERSTROM", "SANDERS", "SFARI1"),
    Fill = RColorBrewer::brewer.pal(n = 3, name = "Set1"),
    TextSize = c(3, 2, 1),
    stringsAsFactors = FALSE
)
gM1_plot <- ggplot() +
    geom_segment(
        data = gM1_edge_layout,
        mapping = aes(
            x = x_head, y = y_head, xend = x_tail, yend = y_tail,
            size = weight, alpha = weight
        ), 
        colour = "darkgrey", lineend = "round", linejoin = "round"
    ) +
    geom_point(
        data = datasets,
        mapping = aes(x = 0, y = 0, fill = Dataset), 
        size = 3, shape = 21, alpha = 0
    ) +
    scale_fill_manual(values = setNames(datasets$Fill, datasets$Dataset)) +
    geom_point(
        data = gM1_vertex_layout %>%
            filter(SATTERSTROM),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SATTERSTROM"], 
        size = datasets$TextSize[datasets$Dataset == "SATTERSTROM"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = gM1_vertex_layout %>%
            filter(SANDERS),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SANDERS"], 
        size = datasets$TextSize[datasets$Dataset == "SANDERS"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = gM1_vertex_layout %>%
            filter(SFARI1),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = datasets$Fill[datasets$Dataset == "SFARI1"], 
        size = datasets$TextSize[datasets$Dataset == "SFARI1"],
        stroke = min(1 / 10),
        show.legend = FALSE
    ) +
    geom_point(
        data = gM1_vertex_layout %>%
            filter(!ASD),
        mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
        fill = "white", 
        stroke = min(1 / 10),
        size = datasets$TextSize[datasets$Dataset == "SFARI1"]
    ) +
    geom_text_repel(
        data = gM1_vertex_layout,
        mapping = aes(
            x = x, y = y, 
            label = ifelse(
                Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
            )
        ), 
        size = 1.75, min.segment.length = 0.1,
        segment.size = min(1 / 10),
    ) +
    geom_segment(
        data = iM1_edge_layout, 
        mapping = aes(
            x = 0, y = 0, xend = 0, yend = 0, alpha = weight, size = weight
        )
    ) +
    scale_size_continuous(
        breaks  = seq(0.7, 1, 0.05), range = c(0.1, 0.5)
    ) +
    scale_alpha_continuous(breaks  = seq(0.7, 1, 0.05)) +
    scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
    guides(
        size = guide_legend(title = "PCC"),
        alpha = guide_legend(title = "PCC"),
        shape = guide_legend(override.aes = list(fill = NA)),
        fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggnetwork::theme_blank() +
    theme(
        text = element_text(size = 8),
        # legend.position = "none"
    )
ggsave(
    filename = "data/ppi/DirectConnect_gM1.pdf",
    plot = gM1_plot + theme(legend.position = "none"),
    device = "pdf", width = 6, height = 4, units = "cm", 
    useDingbats = FALSE
)
ggsave(
    filename = "data/ppi/DirectConnect_legend.pdf",
    plot = cowplot::get_legend(gM1_plot),
    device = "pdf", width = 12, height = 12, units = "cm", 
    useDingbats = FALSE
)

writexl::write_xlsx(
    graph_to_tables(gM1_direct_q10),
    "data/ppi/DirectConnect_gM1_q10.xlsx"
)
writexl::write_xlsx(
    graph_to_tables(iM1_direct_q10),
    "data/ppi/DirectConnect_iM1_q10.xlsx"
)

diffASDImpact <- iM1_direct_q10 %>%
    delete_vertices(
        V(.)[degree(.) < 1]
    ) %>% {
        graph_to_tables(.)$vertices
    } %>%
    group_by(., Gene.Symbol) %>%
    mutate(ctr = 1) %>%
    mutate(n = sum(ctr)) %>%
    filter(Target %in% c("ASD", "Both")) %>%
    mutate(n_target = sum(Target %in% c("ASD", "Both"))) %>%
    filter(n != n_target) %>%
    distinct(Gene.Symbol) %>%
    pull(Gene.Symbol) %>%
    as.character() %>%
    sort()

iM1_direct_q10 %>%
    as_edgelist() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    bind_cols(
        as.data.frame(edge_attr(iM1_direct_q10), stringsAsFactors = FALSE)
    ) %>% 
    setNames(c("head", "tail", "weight")) %>% {
        heads <- as.character(.$head)
        tails <- as.character(.$tail)
        vertex_attributes <- as.data.frame(
            vertex_attr(iM1_direct_q10), stringsAsFactors = FALSE
        )
        head_attr <- setNames(
            vertex_attributes,
            paste0("head_", colnames(vertex_attributes))
        )
        tail_attr <- setNames(
            vertex_attributes,
            paste0("tail_", colnames(vertex_attributes))
        )
        bind_cols(
            ., 
            head_attr[match(heads, head_attr$head_name), ], 
            tail_attr[match(tails, tail_attr$tail_name), ]
        )
    } %>%
    filter(
        head_Gene.Symbol %in% diffASDImpact | tail_Gene.Symbol %in% diffASDImpact
    ) %>% {
        heads <- dplyr::select(., starts_with("head")) %>%
            setNames(str_replace_all(colnames(.), "head", "tail"))
        tails <- dplyr::select(., starts_with("tail")) %>%
            setNames(str_replace_all(colnames(.), "tail", "head"))
        bind_rows(
            ., bind_cols(heads, tails, tibble(weight = .$weight))
        )
    } %>%
    filter(
        head_Gene.Symbol %in% diffASDImpact
    ) %>% 
    dplyr::select(starts_with("head"), weight) %>% {
        target_weights <- filter(., head_Target %in% c("ASD", "Both")) %>%
            pull(weight)
        nontarget_weights <- filter(., ! head_Target %in% c("ASD", "Both")) %>%
            pull(weight)
        wilcox_test <- wilcox.test(
            target_weights, nontarget_weights, 
            paired = FALSE, alternative = "g",
            conf.int = TRUE, conf.level = 0.95
        )
        return(wilcox_test)
    }
    
