library(tidyverse)
library(igraph)

gM1_edges <- read_tsv("data/ppi/genePPI_gM1_edges.tsv")
gM1_vertices <- read_tsv("data/ppi/genePPI_gM1_vertices.tsv")
iM1_edges <- read_tsv("data/ppi/genePPI_iM1_edges.tsv")
iM1_vertices <- read_tsv("data/ppi/genePPI_iM1_vertices.tsv")
iM30_edges <- read_tsv("data/ppi/genePPI_iM30_edges.tsv")
iM30_vertices <- read_tsv("data/ppi/genePPI_iM30_vertices.tsv")

gM1_network <- graph_from_data_frame(gM1_edges, vertices = gM1_vertices)
iM1_network <- graph_from_data_frame(iM1_edges, vertices = iM1_vertices)
iM30_network <- graph_from_data_frame(iM30_edges, vertices = iM30_vertices)

# Remove vertices if they are (not targeted and only have one degree)
# Get degrees of vertices and add as attributes
gM1_network_filter <- gM1_network %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )
iM1_network_filter <- iM1_network %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )
iM30_network_filter <- iM30_network %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )

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

# plot_graph(gM1_network_filter)

network_list <- list(
    gM1 = gM1_network_filter,
    iM1 = iM1_network_filter,
    iM30 = iM30_network_filter
)

lapply(
    names(network_list),
    function(net) {
        write_tsv(
            as.data.frame(vertex_attr(network_list[[net]])),
            paste0("data/ppi/genePPI_", net, "_vertices_filter.tsv")
        )
        write_tsv(
            cbind(
                as_edgelist(network_list[[net]]),
                as.data.frame(edge_attr(network_list[[net]]))
            ),
            paste0("data/ppi/genePPI_", net, "_edges_filter.tsv")
        )
    }
)

# PCC Filter

gM1_network_PCCfilter <- gM1_network %>%
    delete_edges(
        E(.)[edge_attr(., name = "weight") <= 0.9]
    ) %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )
iM1_network_PCCfilter <- iM1_network %>%
    delete_edges(
        E(.)[edge_attr(., name = "weight") <= 0.9]
    ) %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )
iM30_network_PCCfilter <- iM30_network %>%
    delete_vertices(
        v = names(V(.))[! vertex_attr(., name = "Target") %in% c("ASD", "Both") & degree(.) < 2]
    )
network_list_pccFilter <- list(
    gM1 = gM1_network_PCCfilter,
    iM1 = iM1_network_PCCfilter,
    iM30 = iM30_network_PCCfilter
)
lapply(
    names(network_list_pccFilter),
    function(net) {
        write_tsv(
            as.data.frame(vertex_attr(network_list_pccFilter[[net]])),
            paste0("data/ppi/genePPI_", net, "_vertices_PCCfilter.tsv")
        )
        write_tsv(
            cbind(
                as_edgelist(network_list_pccFilter[[net]]),
                as.data.frame(edge_attr(network_list_pccFilter[[net]]))
            ),
            paste0("data/ppi/genePPI_", net, "_edges_PCCfilter.tsv")
        )
    }
)



diff_impact_genes <- c("ARID1B", "CHD8", "KMD5B", "KMT2A", "MED13L", "PCM1", "PHF12", "POGZ", "TCF4")

enst_diff <- as_tibble(vertex_attr(iM1_network)) %>%
    mutate(Differential = `Gene Symbol` %in% diff_impact_genes) %>%
    pull(`Ensembl Transcript ID`)
n <- ego(iM1_network, nodes = enst_diff) %>%
    sapply(., function(n) vertex_attr(iM1_network, name = "Ensembl Transcript ID", index = n)) %>%
    unlist()

net_v <- as_tibble(vertex_attr(iM1_network)) %>%
    mutate(Differential = `Gene Symbol` %in% diff_impact_genes) %>%
    mutate(DifferentialNeighbor = `Ensembl Transcript ID` %in% n) %>%
    mutate(ID = ifelse(Differential | DifferentialNeighbor, `Transcript Symbol`, ""))

write_csv(net_v, "data/slides_network_vertices.csv")
