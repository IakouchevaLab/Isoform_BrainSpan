library(tidyverse)
library(igraph)
library(gprofiler2)

gM1_edges <- readxl::read_xlsx(
    "data/ppi/IsoformPPIAnnotated.xlsx", sheet = "gM1"
)
iM1_edges <- readxl::read_xlsx(
    "data/ppi/IsoformPPIAnnotated.xlsx", sheet = "iM1"
)
iM30_edges <- readxl::read_xlsx(
    "data/ppi/IsoformPPIAnnotated.xlsx", sheet = "iM30"
)
vertices <- readxl::read_xlsx(
    "data/ppi/IsoformPPIAnnotated.xlsx", sheet = "Vertices",
    col_types = "text"
)

gM1_graph <- graph_from_data_frame(
    setNames(gM1_edges, c("F1", "F2", "weight")), 
    vertices = vertices %>%
        filter(name %in% unlist(c(gM1_edges[, 1], gM1_edges[, 2])))
)
iM1_graph <- graph_from_data_frame(
    setNames(iM1_edges, c("F1", "F2", "weight")), 
    vertices = vertices %>%
        filter(name %in% unlist(c(iM1_edges[, 1], iM1_edges[, 2])))
)
iM30_graph <- graph_from_data_frame(
    setNames(iM30_edges, c("F1", "F2", "weight")), 
    vertices = vertices %>%
        filter(name %in% unlist(c(iM30_edges[, 1], iM30_edges[, 2])))
)

gM1_components <- decompose(gM1_graph)
iM1_components <- decompose(iM1_graph)
iM30_components <- decompose(iM30_graph)    # iM30 has nothing

gM1_components_enr <- gM1_components %>%
    setNames(
        nm = sapply(., function(component) {
            vertex_attr(component, name = "Gene Symbol")[[1]]
        })
    ) %>%
    lapply(
        function(component) {
            vertex_order <- order(
                as.numeric(vertex_attr(component, name = "kME")), 
                decreasing = TRUE
            )
            gost(
                query = names(V(component))[vertex_order],
                organism = "hsapiens",
                ordered_query = TRUE,
                exclude_iea = TRUE, 
                correction_method = "fdr",
                custom_bg = vertices$name[
                    vertices$`Module Label` == 1 & vertices$`Data` == "Gene"
                ]
            )$result
        }
    )

iM1_components_enr <- iM1_components %>%
    setNames(
        nm = sapply(., function(component) {
            vertex_attr(component, name = "Transcript Symbol")[[1]]
        })
    ) %>%
    lapply(
        function(component) {
            vertex_order <- order(
                as.numeric(vertex_attr(component, name = "kME")), 
                decreasing = TRUE
            )
            gost(
                query = names(V(component))[vertex_order],
                organism = "hsapiens",
                ordered_query = TRUE,
                exclude_iea = TRUE, 
                correction_method = "fdr",
                custom_bg = vertices$name[
                    vertices$`Module Label` == 1 & vertices$`Data` == "Isoform"
                    ]
            )$result
        }
    )
