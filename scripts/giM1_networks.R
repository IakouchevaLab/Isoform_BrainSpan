library(tidyverse)
library(igraph)
library(ggforce)

plot_graph <- function(this_graph, module_colour) {
    single_nodes <- this_graph %>%
        delete_vertices(
            V(.)[degree(.) > 0]
        )
    single_nodes_layout <- single_nodes %>%
        layout_on_grid() %>%
        as.data.frame() %>%
        setNames(nm = c("x", "y")) %>%
        cbind(
            single_nodes %>%
                vertex_attr() %>%
                as.data.frame() %>%
                filter(name %in% names(V(single_nodes)))
        ) %>%
        mutate(Target = Target %in% c("ASD", "Both"))
    connected_nodes <- this_graph %>%
        delete_vertices(
            V(.)[degree(.) == 0]
        )
    vertex_layout <-  connected_nodes%>%
        layout_with_fr() %>%
        as.data.frame() %>%
        setNames(nm = c("x", "y")) %>%
        cbind(
            connected_nodes %>%
                vertex_attr() %>%
                as.data.frame() %>%
                filter(name %in% names(V(connected_nodes)))
        ) %>%
        mutate(Target = Target %in% c("ASD", "Both"))
    edge_layout <- as.data.frame(as_edgelist(this_graph)) %>%
        setNames(nm = c("F1", "F2")) %>%
        cbind(edge_attr(this_graph)) %>%
        left_join(
            vertex_layout %>%
                dplyr::select(name, x, y),
            by = c("F1" = "name")
        ) %>%
        rename(x1 = x, y1 = y) %>%
        left_join(
            vertex_layout %>%
                dplyr::select(name, x, y),
            by = c("F2" = "name")
        ) %>%
        rename(x2 = x, y2 = y)
    connected_plot <- ggplot() +
        geom_segment(
            data = edge_layout,
            mapping = aes(
                x = x1, y = y1, xend = x2, yend = y2, size = weight
            ),
            colour = "lightgrey"
        ) +
        geom_point(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, shape = ASD, colour = Target
            ),
            size = 5, stroke = 1.5, fill = module_colour
        ) +
        ggrepel::geom_text_repel(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    vertex_layout$Data == "Gene",
                    as.character(Gene.Symbol), as.character(Transcript.Symbol)
                )
            ),
            size = 5
        ) +
        scale_shape_manual(
            values = c("TRUE" = 24, "FALSE" = 21),
            limits = c("TRUE", "FALSE"),
            breaks = c("TRUE", "FALSE")
        ) +
        scale_colour_manual(
            values = c("TRUE" = "red", "FALSE" = "black"),
            limits = c("TRUE", "FALSE")
        ) +
        scale_size_continuous(
            limits = c(
                floor(min(edge_layout$weight) * 10) / 10,
                ceiling(max(edge_layout$weight) * 10) / 10
            ), 
            range = c(0.5, 3)
        ) +
        ggnetwork::theme_blank()
    single_vertex_plot <- ggplot() +
        geom_point(
            data = single_nodes_layout,
            mapping = aes(
                x = x, y = y, shape = ASD, colour = Target
            ),
            size = 5, stroke = 1.5, fill = module_colour
        ) +
        ggrepel::geom_text_repel(
            data = single_nodes_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    single_nodes_layout$Data == "Gene",
                    as.character(Gene.Symbol), as.character(Transcript.Symbol)
                )
            ),
            size = 5
        ) +
        scale_shape_manual(
            values = c("TRUE" = 24, "FALSE" = 21),
            limits = c("TRUE", "FALSE"),
            breaks = c("TRUE", "FALSE")
        ) +
        scale_colour_manual(
            values = c("TRUE" = "red", "FALSE" = "black"),
            limits = c("TRUE", "FALSE")
        ) +
        scale_size_continuous(
            limits = c(
                floor(min(edge_layout$weight) * 10) / 10,
                ceiling(max(edge_layout$weight) * 10) / 10
            ), 
            range = c(0.5, 3)
        ) +
        ggnetwork::theme_blank()
    print(connected_plot)
}
module_colour <- "turquoise"
gM1_edges <- read_tsv("data/ppi/genePPI_gM1_edges.tsv")
gM1_vertices <- read_tsv("data/ppi/genePPI_gM1_vertices.tsv")
iM1_edges <- read_tsv("data/ppi/genePPI_iM1_edges.tsv")
iM1_vertices <- read_tsv("data/ppi/genePPI_iM1_vertices.tsv")

gM1_network <- graph_from_data_frame(
    gM1_edges, vertices = gM1_vertices, directed = FALSE
)
iM1_network <- graph_from_data_frame(
    iM1_edges, vertices = iM1_vertices, directed = FALSE
)

# Differential networks
gM1_diffNetwork <- delete_vertices(
    gM1_network,
    v = V(gM1_network)[(
        vertex_attr(gM1_network, name = "Ensembl Gene ID") %in% 
            vertex_attr(iM1_network, name = "Ensembl Gene ID")
    )]
) %>%
    delete_vertices(
        v = V(.)[degree(., V(.)) < 1]
    )
iM1_diffNetwork <- delete_vertices(
    iM1_network,
    v = V(iM1_network)[(
        vertex_attr(iM1_network, name = "Ensembl Gene ID") %in% 
            vertex_attr(gM1_network, name = "Ensembl Gene ID")
    )]
)

# Components with targeted ASD genes
iM1_ASDComponents <- decompose(iM1_diffNetwork) %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD & (V(net)$Target %in% c("ASD", "Both")))) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    plyr::compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[degree(., V(.)) < 1]
    ) %>%
    decompose()

# Functional enrichment analysis of iM1_ASD_DiffNet
iM1_ASD_DiffNet_Components_GO <- lapply(
    iM1_ASDComponents,
    function(iM1_ASD_component) {
        gprofiler2::gost(
            query = names(V(iM1_ASD_component))[order(
                V(iM1_ASD_component)$kME, decreasing = TRUE
            ) & V(iM1_ASD_component)$`Gene Symbol` != "SCN2A"],
            organism = "hsapiens",
            ordered_query = TRUE, 
            exclude_iea = TRUE,
            correction_method = "fdr",
            custom_bg = rownames(readRDS("data/RegressIsoformCounts.rds")),
            sources = c("GO:BP", "GO:MF")
        )$result
    }
)

# Annotate components of graph
component_functions <- list(
    c(
        "regulation of NMDA receptor activity", 
        "regulation of glutamate receptor signaling pathway",
        "synaptic signaling"
    ),
    c(
        "cytoskeletal anchoring at plasma membrane", 
        "transmission of nerve impulse",
        "maintenance of protein location"
    )
) %>%
    lapply(paste, collapse = "\n")


# Plot this graph
iM1_ASDComponents_merge <- mapply(
    function(component, func) {
        set_vertex_attr(component, name = "Function", value = func)
    },
    iM1_ASDComponents, component_functions, SIMPLIFY = FALSE
) %>%
    disjoint_union()
iM1_ASDComponents_grouped <- iM1_ASDComponents_merge
for (i in unique(V(iM1_ASDComponents_grouped)$`Ensembl Gene ID`)) {
    groupV <- which(V(iM1_ASDComponents_grouped)$`Ensembl Gene ID` == i)
    if (length(groupV) < 2) { next }
    iM1_ASDComponents_grouped <- add_edges(
        iM1_ASDComponents_grouped, combn(groupV, 2), attr = list(weight = 1.1)
    )
}
vertex_layout <- iM1_ASDComponents_grouped %>%
    layout_with_dh() %>%
    as.data.frame() %>%
    setNames(nm = c("x", "y")) %>%
    mutate(name = names(V(iM1_ASDComponents_merge))) %>%
    left_join(
        as.data.frame(vertex_attr(iM1_ASDComponents_merge)),
        by = "name"
    ) %>%
    mutate(Target = Target %in% c("ASD", "Both"))
colnames(vertex_layout) <- make.names(colnames(vertex_layout))
vertex_layout <- vertex_layout %>%
    mutate(counter = 1) %>%
    group_by(Ensembl.Gene.ID) %>%
    mutate(n = sum(counter)) %>%
    ungroup()
edge_layout <- as.data.frame(as_edgelist(iM1_ASDComponents_merge)) %>%
    setNames(nm = c("F1", "F2")) %>%
    cbind(as.data.frame(edge_attr(iM1_ASDComponents_merge))) %>%
    left_join(
        vertex_layout %>%
            dplyr::select(name, x, y),
        by = c("F1" = "name")
    ) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(
        vertex_layout %>%
            dplyr::select(name, x, y),
        by = c("F2" = "name")
    ) %>%
    rename(x2 = x, y2 = y)
iM1_ASDComponents_Grouped_Plot <- ggplot() +
    geom_mark_rect(
        data = vertex_layout %>%
            filter(n > 1),
        mapping = aes(
            x = x, y = y, description = Function
        ), 
        con.cap = 0, fill = "lightgrey",
        label.fontsize = 16, label.fill = NA
    ) +
    geom_segment(
        data = edge_layout,
        mapping = aes(
            x = x1, y = y1, xend = x2, yend = y2, 
            size = abs(weight), alpha = abs(weight)
        ),
        colour = "grey"
    ) +
    geom_point(
        data = vertex_layout,
        mapping = aes(
            x = x, y = y, shape = ASD, colour = Target
        ),
        size = 5, stroke = 1.1, fill = module_colour
    ) +
    geom_point(
        data = vertex_layout %>%
            filter(Target),
        mapping = aes(
            x = x, y = y, shape = ASD
        ),
        size = 5, stroke = 1.1, fill = "red", colour = "red"
    ) +
    geom_mark_ellipse(
        data = vertex_layout %>%
            filter(n > 1),
        mapping = aes(
            x = x, y = y, label = Gene.Symbol, group = Gene.Symbol
        ), 
        fill = "turquoise", con.cap = 0,
        label.fontsize = 14, label.hjust = 0.5, label.fill = NA
    ) +
    ggrepel::geom_text_repel(
        data = vertex_layout,
        mapping = aes(
            x = x, y = y, 
            label = ifelse(
                vertex_layout$Data == "Gene",
                as.character(Gene.Symbol), as.character(Transcript.Symbol)
            )
        ),
        size = 4
    ) +
    scale_shape_manual(
        values = c("TRUE" = 24, "FALSE" = 21),
        limits = c("TRUE", "FALSE"),
        breaks = c("TRUE", "FALSE")
    ) +
    scale_colour_manual(
        values = c("TRUE" = "red", "FALSE" = "black"),
        limits = c("TRUE", "FALSE")
    ) +
    scale_size_continuous(
        # limits = c(
        #     floor(min(edge_layout$weight) * 10) / 10,
        #     ceiling(max(edge_layout$weight) * 10) / 10
        # ), 
        range = c(0.05, 3)
    ) +
    scale_alpha_continuous(
        range = c(0.01, 1)
    ) +
    ggnetwork::theme_blank()

ggsave(
    filename = "data/ppi/iM1_DifferentialNetwork_SatterstromASDComponent.pdf",
    plot = iM1_ASDComponents_Grouped_Plot,
    device = "pdf", width = 12, height = 16,
    useDingbats = FALSE
)

# Components with targeted ASD genes (Nothing here)
gM1_ASDComponents <- decompose(gM1_diffNetwork) %>%
    lapply(
        function(net) {
            if (any(V(net)$ASD & (V(net)$Target %in% c("ASD", "Both")))) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    plyr::compact() %>%
    disjoint_union() %>%
    delete_vertices(
        V(.)[degree(., V(.)) < 1]
    ) %>%
    decompose()

# Overlapping networks
gM1_overlapping <- gM1_network %>%
    delete_vertices(
        V(.)[which(
            !V(.)$`Ensembl Gene ID` %in% V(iM1_network)$`Ensembl Gene ID`
        )]
    )
iM1_overlapping <- iM1_network %>%
    delete_vertices(
        V(.)[which(
            !V(.)$`Ensembl Gene ID` %in% V(gM1_network)$`Ensembl Gene ID`
        )]
    )


iM1Overlap_Components <- iM1_overlapping %>%
    delete_edges(E(.)[abs(E(.)$weight) <= 0.9]) %>%
    delete_vertices(V(.)[degree(.) < 1]) %>%
    decompose() %>%
    lapply(
        function(net) {
            if (any(V(net)$`ASD` & V(net)$`Target` %in% c("ASD", "Both"))) {
                return(net)
            } else {
                return(NULL)
            }
        }
    ) %>%
    compact()
plot_graph(iM1Overlap_Components[[2]], "turquoise")
