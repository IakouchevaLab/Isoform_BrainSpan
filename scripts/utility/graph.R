library(ggrepel)

# Utility functions
graph_to_tables <- function(g) {
    vertex_table <- as.data.frame(vertex_attr(g))
    edge_table <- setNames(
        as.data.frame(as_edgelist(g)), nm = c("head", "tail")
    ) %>%
        bind_cols(as.data.frame(edge_attr(g)))
    return(list(vertices = vertex_table, edges = edge_table))
}

plot_graph <- function(g, colour, plot_singletons = FALSE) {
    # Create "dummy" edges among singleton nodes
    singletons <- V(g)[which(degree(g) < 1)]
    g_components <- g %>%
        delete_vertices(V(.)[which(V(.) %in% singletons)]) %>%
        set_edge_attr(name = "weight", value = abs(E(.)$weight))
    g_singletons <- g %>%
        delete_vertices(V(.)[which(! V(.) %in% singletons)])
    
    singleton_vertex_layout <- g_singletons %>%
        layout_on_grid() %>%
        as.data.frame(stringsAsFactors = FALSE) %>% 
        setNames(nm = c("x", "y")) %>%
        mutate(name = names(V(g_singletons))) %>%
        inner_join(
            as.data.frame(
                vertex_attr(g_singletons), stringsAsFactors = FALSE
            ), by = "name"
        ) %>%
        mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
    vertex_layout <- g_components %>%
        layout_with_dh() %>%
        as.data.frame(stringsAsFactors = FALSE) %>% 
        setNames(nm = c("x", "y")) %>%
        mutate(name = names(V(g_components))) %>%
        inner_join(
            as.data.frame(
                vertex_attr(g_components), stringsAsFactors = FALSE
            ), by = "name"
        ) %>%
        mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
    edge_layout <- as_edgelist(g_components) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        setNames(nm = c("head", "tail")) %>% 
        bind_cols(
            as.data.frame(edge_attr(g))
        ) %>% {
            heads <- as.character(.$head)
            tails <- as.character(.$tail)
            head_coords <- vertex_layout[match(heads, vertex_layout$name), ] %>%
                dplyr::select(x, y) %>%
                setNames(nm = paste(names(.), "head", sep = "_"))
            tail_coords <- vertex_layout[match(tails, vertex_layout$name), ] %>%
                dplyr::select(x, y) %>%
                setNames(nm = paste(names(.), "tail", sep = "_"))
            bind_cols(., head_coords, tail_coords)
        }
    components_plot <- ggplot() +
        geom_segment(
            data = edge_layout,
            mapping = aes(
                x = x_head, y = y_head, xend = x_tail, yend = y_tail,
                size = weight, alpha = weight
            ), 
            colour = "grey", lineend = "round", linejoin = "round"
        ) +
        geom_point(
            data = vertex_layout,
            mapping = aes(x = x, y = y, shape = ASDLoFTarget, fill = ASD), 
            size = 5
        ) +
        geom_text_repel(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
                )
            ), 
            size = 5
        ) +
        scale_size_continuous(range = c(0.01, 3)) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = colour)) +
        scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
        guides(
            size = guide_legend(title = "PCC"),
            alpha = guide_legend(title = "PCC")
        ) +
        ggnetwork::theme_blank()
    singletons_plot <- ggplot() +
        geom_point(
            data = singleton_vertex_layout,
            mapping = aes(x = x, y = y, shape = ASDLoFTarget, fill = ASD), 
            size = 5
        ) +
        geom_text_repel(
            data = singleton_vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
                )
            ), 
            size = 5
        ) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = colour)) +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        scale_alpha_continuous(
            range = c(0.15, 1),
            limits = c(
                floor(min(edge_layout$weight) * 10) / 10,
                ceiling(max(edge_layout$weight) * 10) / 10
            )
        ) +
        ggnetwork::theme_blank()
    if (plot_singletons) {
        return(
            cowplot::plot_grid(
                cowplot::plot_grid(
                    components_plot + theme(legend.position = "none"), 
                    singletons_plot + theme(legend.position = "none"),
                    ncol = 1, nrow = 2, rel_heights = c(1, 0.5)
                ),
                cowplot::get_legend(components_plot), 
                ncol = 2, nrow = 1, rel_widths = c(1, 0.2)
            )
        )
    } else {
        return(components_plot)
    }
}

plot_graph_byData <- function(
    g, 
    colour, 
    labelSize = 5,
    pointSize = c(1, 2, 3),
    plot_singletons = FALSE
) {
    # Create "dummy" edges among singleton nodes
    singletons <- V(g)[which(degree(g) < 1)]
    g_components <- g %>%
        delete_vertices(V(.)[which(V(.) %in% singletons)]) %>%
        set_edge_attr(name = "weight", value = abs(E(.)$weight))
    g_singletons <- g %>%
        delete_vertices(V(.)[which(! V(.) %in% singletons)])
    
    singleton_vertex_layout <- g_singletons %>%
        layout_on_grid() %>%
        as.data.frame(stringsAsFactors = FALSE) %>% 
        setNames(nm = c("x", "y")) %>%
        mutate(name = names(V(g_singletons))) %>%
        inner_join(
            as.data.frame(
                vertex_attr(g_singletons), stringsAsFactors = FALSE
            ), by = "name"
        ) %>%
        mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
    vertex_layout <- g_components %>%
        layout_with_dh() %>%
        as.data.frame(stringsAsFactors = FALSE) %>% 
        setNames(nm = c("x", "y")) %>%
        mutate(name = names(V(g_components))) %>%
        inner_join(
            as.data.frame(
                vertex_attr(g_components), stringsAsFactors = FALSE
            ), by = "name"
        ) %>%
        mutate(ASDLoFTarget = Target %in% c("ASD", "Both"))
    edge_layout <- as_edgelist(g_components) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        setNames(nm = c("head", "tail")) %>% 
        bind_cols(
            as.data.frame(edge_attr(g))
        ) %>% {
            heads <- as.character(.$head)
            tails <- as.character(.$tail)
            head_coords <- vertex_layout[match(heads, vertex_layout$name), ] %>%
                dplyr::select(x, y) %>%
                setNames(nm = paste(names(.), "head", sep = "_"))
            tail_coords <- vertex_layout[match(tails, vertex_layout$name), ] %>%
                dplyr::select(x, y) %>%
                setNames(nm = paste(names(.), "tail", sep = "_"))
            bind_cols(., head_coords, tail_coords)
        }
    datasets <- data.frame(
        Dataset = c("SATTERSTROM", "SANDERS", "SFARI1"),
        Fill = RColorBrewer::brewer.pal(n = 3, name = "Set1"),
        TextSize = sort(pointSize, decreasing = TRUE),
        stringsAsFactors = FALSE
    )
    components_plot <- ggplot() +
        geom_segment(
            data = edge_layout,
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
            data = vertex_layout %>%
                filter(SATTERSTROM),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SATTERSTROM"], 
            size = datasets$TextSize[datasets$Dataset == "SATTERSTROM"],
            stroke = min(pointSize / 10),
            show.legend = FALSE
        ) +
        geom_point(
            data = vertex_layout %>%
                filter(SANDERS),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SANDERS"], 
            size = datasets$TextSize[datasets$Dataset == "SANDERS"],
            stroke = min(pointSize / 10),
            show.legend = FALSE
        ) +
        geom_point(
            data = vertex_layout %>%
                filter(SFARI1),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SFARI1"], 
            size = datasets$TextSize[datasets$Dataset == "SFARI1"],
            stroke = min(pointSize / 10),
            show.legend = FALSE
        ) +
        geom_point(
            data = vertex_layout %>%
                filter(!ASD),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = colour, 
            stroke = min(pointSize / 10),
            size = datasets$TextSize[datasets$Dataset == "SFARI1"]
        ) +
        geom_text_repel(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
                )
            ), 
            size = labelSize,
            segment.size = min(pointSize / 10),
        ) +
        scale_size_continuous(range = c(0.01, 3)) +
        scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
        scale_alpha_continuous(range = c(0.15, 1)) +
        guides(
            size = guide_legend(title = "PCC"),
            alpha = guide_legend(title = "PCC"),
            shape = guide_legend(override.aes = list(fill = NA)),
            fill = guide_legend(override.aes = list(alpha = 1))
        ) +
        ggnetwork::theme_blank() +
        theme(
            text = element_text(size = 5)
        )
    singletons_plot <- ggplot() +
        geom_point(
            data = singleton_vertex_layout %>%
                filter(SATTERSTROM),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SATTERSTROM"], size = 3
        ) +
        geom_point(
            data = singleton_vertex_layout %>%
                filter(SANDERS),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SANDERS"], size = 2
        ) +
        geom_point(
            data = singleton_vertex_layout %>%
                filter(SFARI1),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = datasets$Fill[datasets$Dataset == "SFARI1"], size = 1
        ) +
        geom_point(
            data = singleton_vertex_layout %>%
                filter(!ASD),
            mapping = aes(x = x, y = y, shape = ASDLoFTarget), 
            fill = colour, size = 4
        ) +
        geom_text_repel(
            data = singleton_vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    Data == "Gene", `Gene.Symbol`, `Transcript.Symbol`
                )
            ), 
            size = 1
        ) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = colour)) +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21)) +
        scale_alpha_continuous(range = c(0.15, 1)) +
        ggnetwork::theme_blank()
    if (plot_singletons) {
        return(
            cowplot::plot_grid(
                cowplot::plot_grid(
                    components_plot + theme(legend.position = "none"), 
                    singletons_plot + theme(legend.position = "none"),
                    ncol = 1, nrow = 2, rel_heights = c(1, 0.5)
                ),
                cowplot::get_legend(components_plot), 
                ncol = 2, nrow = 1, rel_widths = c(1, 0.2)
            )
        )
    } else {
        return(components_plot)
    }
}

gene_isoform_overlapping_edges <- function(gene, isoform) {
    g_edges <- as_edgelist(gene) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        setNames(nm = c("head", "tail")) %>%
        mutate(edge_string = mapply(
            function(h, t) paste(sort(c(h, t)), collapse = "", sep = ""),
            head, tail
        ))
    i_edges <- as_edgelist(set_vertex_attr(
        isoform, name = "name", value = V(isoform)$`Ensembl Gene ID`
    )) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        setNames(nm = c("head", "tail")) %>%
        mutate(edge_string = mapply(
            function(h, t) paste(sort(c(h, t)), collapse = "", sep = ""),
            head, tail
        ))
    g_edges <- g_edges %>%
        mutate(overlap = edge_string %in% i_edges$edge_string)
    i_edges <- i_edges %>%
        mutate(overlap = edge_string %in% g_edges$edge_string)
    return(list(
        gene_edges = E(gene)[g_edges$overlap],
        isoform_edges = E(isoform)[i_edges$overlap]
    ))
}


# Filter uninteresting nodes without edges
filter_uninteresting <- function(g) {
    uninteresting_nodes <- V(g)[which(
        !V(g)$ASD & !V(g)$Target %in% c("ASD", "Both")
    )]
    uninteresting_nodes_degree <- degree(g, v = uninteresting_nodes)
    remove_nodes <- uninteresting_nodes[uninteresting_nodes_degree < 1]
    g %>%
        delete_vertices(remove_nodes) %>%
        decompose() %>%
        lapply(function(net) {
            if (any(V(net)$ASD | V(net)$Target %in% c("ASD", "Both"))) {
                return(net)
            } else {
                return(NULL)
            }
        }) %>%
        compact() %>%
        disjoint_union() %>%
        return()
}

# Filter non-ASD genes that connect to at least n ASD genes
filter_by_connection_ASD <- function(g, n) {
    non_asd_non_target <- V(g)[which(
        !V(g)$ASD & !V(g)$Target %in% c("ASD", "Both")
    )]
    num_connections <- sapply(
        non_asd_non_target,
        function(v) {
            v_neighbors <- neighbors(g, v = v, mode = "all")
            interesting_neighbors <- vertex_attr(
                graph = g, index = v_neighbors
            ) %>%
                as.data.frame() %>%
                filter(ASD)
            return(nrow(interesting_neighbors))
        }
    )
    remove_nodes <- non_asd_non_target[num_connections < n]
    return(delete_vertices(g, remove_nodes))
}

# Filter non-ASD genes that connect to at least n target genes
filter_by_connection_target <- function(g, n) {
    non_target_non_target <- V(g)[which(
        !V(g)$ASD & !V(g)$Target %in% c("ASD", "Both")
    )]
    num_connections <- sapply(
        non_target_non_target,
        function(v) {
            v_neighbors <- neighbors(g, v = v, mode = "all")
            interesting_neighbors <- vertex_attr(
                graph = g, index = v_neighbors
            ) %>%
                as.data.frame() %>%
                filter(Target %in% c("ASD", "Both"))
            return(nrow(interesting_neighbors))
        }
    )
    remove_nodes <- non_target_non_target[num_connections < n]
    return(delete_vertices(g, remove_nodes))
}

# Filter non-ASD genes that connect to at least n ASD and/or target genes
filter_by_connection_ASD_target <- function(g, n) {
    non_asd_non_target <- V(g)[which(
        !V(g)$ASD & !V(g)$Target %in% c("ASD", "Both")
    )]
    num_connections <- sapply(
        non_asd_non_target,
        function(v) {
            v_neighbors <- neighbors(g, v = v, mode = "all")
            interesting_neighbors <- vertex_attr(
                graph = g, index = v_neighbors
            ) %>%
                as.data.frame(stringsAsFactors = FALSE) %>%
                filter(Target %in% c("ASD", "Both") | ASD)
            return(nrow(interesting_neighbors))
        }
    )
    remove_nodes <- non_asd_non_target[num_connections < n]
    return(delete_vertices(g, remove_nodes))
}

# Filter non-ASD genes that connect to at least n SFARI1 genes
filter_by_connection_SFARI1 <- function(g, n) {
    v_uninteresting <- V(g)[which(
        !V(g)$ASD & !V(g)$Target %in% c("ASD", "Both")
    )]
    v_neighbors <- ego(g, order = 1, nodes = v_uninteresting, mode = "all") %>%
        lapply(
            function(v) {
                vertex_attr(g, index = v) %>%
                    as.data.frame(stringsAsFactors = FALSE) %>%
                    filter(SFARI1) %>%
                    nrow()
            }
        )
    remove_nodes <- v_uninteresting[v_neighbors < n]
    return(delete_vertices(g, remove_nodes))
}
