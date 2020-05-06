library(tidyverse)
library(igraph)
library(WGCNA)
library(ggrepel)

################################################################################
# Globals                                                                      #
################################################################################

n_hubs <- 20

################################################################################
# Data                                                                         #
################################################################################

transcript_annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", 
    header = TRUE, sep = ",", row.names = 1
)

gexpr <- readRDS("data/RegressGeneCounts.rds")
gsft <- readRDS("data/genes/SoftThresholding.rds")$powerEstimate
iexpr <- readRDS("data/RegressIsoformCounts.rds")
isft <- readRDS("data/isoforms/SoftThresholding.rds")$powerEstimate

module_assignments <- bind_rows(
    read_tsv("data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        mutate(Data = "Gene") %>%
        distinct(),
    read_tsv("data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        mutate(Data = "Isoform") %>%
        distinct()
) %>%
    left_join(
        bind_rows(
            transcript_annotations %>% 
                dplyr::select(ensembl_gene_id, external_gene_id) %>%
                rename(Feature = ensembl_gene_id, Symbol = external_gene_id),
            transcript_annotations %>% 
                dplyr::select(
                    ensembl_transcript_id, external_transcript_id
                ) %>%
                rename(
                    Feature = ensembl_transcript_id, 
                    Symbol = external_transcript_id
                )
        ) %>%
            distinct(), 
        by = "Feature"
    )

################################################################################
# Analysis                                                                     #
################################################################################

# Get hubs per module, per dataset

hubs <- lapply(
    setNames(nm = sort(unique(module_assignments$Data))),
    function(data_type) {
        modules <- module_assignments %>%
            filter(Data == data_type) %>%
            pull(module_label) %>%
            as.numeric() %>% sort() %>% unique()
        lapply(
            setNames(nm = modules),
            function(module) {
                kME_col <- paste0("kME", module)
                this_hubs <- module_assignments %>%
                    filter(Data == data_type & module_label == module) %>%
                    top_n(n_hubs, !!sym(kME_col)) %>%
                    pull(Symbol)
                expand.grid(
                    source = this_hubs,
                    target = this_hubs
                ) %>%
                    mutate(Data = data_type, module_label = module) %>%
                    filter(source != target) %>%
                    mutate(
                        ensembl_source = module_assignments$Feature[
                            match(source, module_assignments$Symbol)
                        ],
                        ensembl_target = module_assignments$Feature[
                            match(target, module_assignments$Symbol)
                        ]
                    )
            }
        ) %>%
            bind_rows()
    }
) %>%
    bind_rows()

symbol_ensembl_map <- c(
    setNames(hubs$source, hubs$ensembl_source), 
    setNames(hubs$target, hubs$ensembl_target)
)

networks_all <- lapply(
    setNames(nm = sort(unique(hubs$Data))),
    function(data_type) {
        modules <- sort(unique(filter(hubs, Data == data_type)$module_label))
        lapply(
            setNames(nm = modules),
            function(module) {
                if (data_type == "Gene") {
                    expr <- gexpr
                    sft <- gsft
                } else {
                    expr <- iexpr
                    sft <- isft
                }
                rownames(expr) <- symbol_ensembl_map[rownames(expr)]
                this_hubs <- hubs %>%
                    filter(Data == data_type & module_label == module) %>%
                    dplyr::select(source, target) %>%
                    as.matrix() %>% as.character() %>% unique()
                adj_mat <- adjacency(
                    t(expr[this_hubs, ]), type = "signed", power = sft
                )
                adj_mat[adj_mat < quantile(adj_mat, 0.1)] = 0
                graph.adjacency(
                    as.matrix(adj_mat), mode = "undirected", weighted=T, diag=F
                )
            }
        )
    }
)

saveRDS(
    networks_all, 
    file = paste0("data/WGCNAHubNetworks_n", n_hubs, ".rds")
)

network_plots <- lapply(
    setNames(nm = names(networks_all)),
    function(data_type) {
        lapply(
            setNames(nm = names(networks_all[[data_type]])),
            function(module) {
                module_colour <- module_assignments %>%
                    filter(Data == data_type & module_label == module) %>%
                    pull(module_colour) %>%
                    unique()
                this_graph <- networks_all[[data_type]][[module]]
                plotcoord <- data.frame(layout_with_fr(this_graph))
                colnames(plotcoord) <- c("X1","X2")
                edgelist <- get.edgelist(this_graph, names = F)
                edges <- data.frame(
                    plotcoord[edgelist[, 1], ], plotcoord[edgelist[, 2], ]
                )
                colnames(edges) <- c("X1","Y1","X2","Y2")
                plotcoord = cbind(
                    plotcoord, 
                    data.frame(symbol = as_ids(V(this_graph)))
                )
                ggplot() + 
                    geom_segment(
                        data = edges,
                        mapping = aes(x = X1, y = Y1, xend = X2, yend = Y2),
                        size = 0.5, 
                        colour="grey"
                    ) +
                    geom_point(
                        data = plotcoord,
                        mapping = aes(x = X1, y = X2), 
                        colour = "black", fill = module_colour, 
                        size = 5, shape = 21
                    ) + 
                    geom_text_repel(
                        data = plotcoord,
                        mapping = aes(x = X1, y = X2 + .05, label = symbol),
                        fontface="bold", size = 8,
                        vjust = "inward", hjust = "inward"
                    ) +
                    labs(
                        title = paste(
                            data_type, module, module_colour
                        )
                    ) +
                    scale_x_continuous(expand = c(0.5, 0)) +
                    scale_y_continuous(expand = c(0.5, 0)) +
                    theme_classic() +
                    theme(
                        text = element_text(size = 30),
                        axis.line = element_blank(), 
                        axis.title = element_blank(),
                        axis.text = element_blank(), 
                        axis.ticks = element_blank()
                    )
            }
        )
    }
)

saveRDS(
    network_plots,
    "data/figures/WGCNA_networks.rds"
)

pdf("data/figures/WGCNA_networks.pdf", width = 16, height = 12)
invisible(lapply(network_plots, print))
dev.off()

pdf("data/genes/figures/WGCNA_networks.pdf", width = 16, height = 12)
invisible(lapply(network_plots[["Gene"]], print))
dev.off()

pdf("data/isoforms/figures/WGCNA_networks.pdf", width = 16, height = 12)
invisible(lapply(network_plots[["Isoform"]], print))
dev.off()
