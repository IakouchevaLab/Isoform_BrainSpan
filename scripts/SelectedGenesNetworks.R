library(tidyverse)
library(igraph)

genes <- c(
    "ADNP", "ARID1B", "DYRK1A", "CHD2", "ASH1L", "ASXL3", "CHD8", "ANK2", 
    "CUL3", "DSCAM", "MYT1L", "KMT2A", "KMT5B", "GRIN2B", "POGZ", "KATNAL2", 
    "NAA15", "PTEN", "TRIP12", "RELN", "TBR1", "SCN2A", "SYNGAP1", "SETD5", 
    "SHANK3"
)

# Selected modules

gM1_edges <- read_tsv("data/ppi/genePPI_gM1_edges.tsv")
gM1_vertices <- read_tsv("data/ppi/genePPI_gM1_vertices.tsv") %>%
    mutate(SFARI1 = `Gene Symbol` %in% genes)
iM1_edges <- read_tsv("data/ppi/genePPI_iM1_edges.tsv")
iM1_vertices <- read_tsv("data/ppi/genePPI_iM1_vertices.tsv") %>%
    mutate(SFARI1 = `Gene Symbol` %in% genes)
iM30_edges <- read_tsv("data/ppi/genePPI_iM30_edges.tsv")
iM30_vertices <- read_tsv("data/ppi/genePPI_iM30_vertices.tsv") %>%
    mutate(SFARI1 = `Gene Symbol` %in% genes)

# How many 102ASD, SFARI1 genes per module
data.frame(
    Module = c("gM1", "iM1", "iM30"),
    SFARI1 = c(
        nrow(filter(gM1_vertices, SFARI1)), 
        nrow(distinct(filter(iM1_vertices, SFARI1), `Gene Symbol`)),
        nrow(distinct(filter(iM30_vertices, SFARI1), `Gene Symbol`))
    ),
    ASD = c(
        nrow(filter(gM1_vertices, ASD)), 
        nrow(distinct(filter(iM1_vertices, ASD), `Gene Symbol`)),
        nrow(distinct(filter(iM30_vertices, ASD), `Gene Symbol`))
    )
)
writexl::write_xlsx(
    list(
        "gM1" = filter(gM1_vertices, SFARI1 | ASD),
        "iM1" = filter(iM1_vertices, SFARI1 | ASD),
        "iM30" = filter(iM30_vertices, SFARI1 | ASD)
    ),
    "data/ppi/SFARI1_ASD_Vertices.xlsx"
)




gM1_graph <- graph_from_data_frame(
    d = gM1_edges, directed = FALSE, vertices = gM1_vertices
)
iM1_graph <- graph_from_data_frame(
    d = iM1_edges, directed = FALSE, vertices = iM1_vertices
)
iM30_graph <- graph_from_data_frame(
    d = iM30_edges, directed = FALSE, vertices = iM30_vertices
)

all_graphs <- list(
    gM1 = gM1_graph,
    iM1 = iM1_graph,
    iM30 = iM30_graph
)

# Select components containing selected genes

select_components_by_gene_symbol <- function(g, v) {
    decompose(g) %>%
        lapply(
            function(g) {
                symbols <- vertex_attr(g, name = "Gene Symbol")
                if (any(symbols %in% genes)) {
                    return(g)
                } else {
                    return(NULL)
                }
            }
        ) %>%
        plyr::compact() %>%
        disjoint_union()
}

all_graphs_selected_genes <- lapply(
    all_graphs, select_components_by_gene_symbol, v = genes
)

# Write to file function

write_graph_to_file <- function(g, gnm, write_to_dir, suffix) {
    write_tsv(
        as.data.frame(vertex_attr(g)),
        file.path(
            write_to_dir, paste0("genePPI_", gnm, "_vertices_", suffix, ".tsv")
        )
    )
    write_tsv(
        cbind(
            as_edgelist(g),
            as.data.frame(edge_attr(g))
        ),
        file.path(
            write_to_dir, paste0("genePPI_", gnm, "_edges_", suffix, ".tsv")
        )
    )
    return(0)
}
mapply(
    function(g, gnm) {
        write_graph_to_file(g, gnm, "data/ppi", "selected_genes")
    }, all_graphs_selected_genes, names(all_graphs_selected_genes)
)

# Filter neighbors (untargeted connected to at least 2 other targeted)
all_graphs_selected_genes_filter <- all_graphs_selected_genes %>%
    lapply(
        function(g) {
            untargeted_nodes <- names(V(g))[
                ! vertex_attr(g, name = "Target") %in% c("ASD", "Both")
            ]
            num_targeted_neighbors <- setNames(nm = untargeted_nodes) %>%
                sapply(
                    function(v) {
                        as.data.frame(vertex_attr(g)) %>%
                            filter(
                                name %in% names(unlist(
                                    ego(g, order = 1, nodes = v, mode = "all")
                                ))
                            ) %>%
                            filter(Target %in% c("ASD", "Both")) %>%
                            nrow()
                    }
                )
            keep_nodes <- c(
                untargeted_nodes[num_targeted_neighbors >= 2],
                names(V(g))[vertex_attr(g, "Target") %in% c("ASD", "Both")]
            )
            delete_vertices(g, V(g)[! names(V(g)) %in% keep_nodes])
        }
    ) %>%
    lapply(
        select_components_by_gene_symbol, v = genes
    )
mapply(
    function(g, gnm) {
        write_graph_to_file(g, gnm, "data/ppi", "selected_genes_filter")
    }, all_graphs_selected_genes_filter, names(all_graphs_selected_genes_filter)
)

# Distinct graphs between gM1 and iM1

gM1_selected_filter_distinct <- all_graphs_selected_genes_filter$gM1 %>%
    delete_vertices(
        V(.)[
            vertex_attr(., name = "Ensembl Gene ID") %in% vertex_attr(
                all_graphs$iM1, name = "Ensembl Gene ID"
            )]
    )
iM1_selected_filter_distinct <- all_graphs_selected_genes_filter$iM1 %>%
    delete_vertices(
        V(.)[
            vertex_attr(., name = "Ensembl Gene ID") %in% vertex_attr(
                all_graphs$gM1, name = "Ensembl Gene ID"
            )]
    )
distinct_igM1 <- lapply(
    list(
        gM1 = gM1_selected_filter_distinct, iM1 = iM1_selected_filter_distinct
    ),
    function(g) {
        untargeted_nodes <- names(V(g))[
            ! vertex_attr(g, name = "Target") %in% c("ASD", "Both")
            ]
        num_targeted_neighbors <- setNames(nm = untargeted_nodes) %>%
            sapply(
                function(v) {
                    as.data.frame(vertex_attr(g)) %>%
                        filter(
                            name %in% names(unlist(
                                ego(g, order = 1, nodes = v, mode = "all")
                            ))
                        ) %>%
                        filter(Target %in% c("ASD", "Both")) %>%
                        nrow()
                }
            )
        keep_nodes <- c(
            untargeted_nodes[num_targeted_neighbors >= 2],
            names(V(g))[vertex_attr(g, "Target") %in% c("ASD", "Both")]
        )
        delete_vertices(g, V(g)[! names(V(g)) %in% keep_nodes])
    }
) %>%
    lapply(
        select_components_by_gene_symbol, v = genes
    )
mapply(
    function(g, gnm) {
        write_graph_to_file(g, gnm, "data/ppi", "selected_filter_distinct")
    }, distinct_igM1, names(distinct_igM1)
)
