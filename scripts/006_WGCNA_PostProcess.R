library(WGCNA)
library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]

gene_expr <- readRDS("data/RegressGeneCounts.rds")
isoform_expr <- readRDS("data/RegressIsoformCounts.rds")

process_network <- function(net, expr, data_type = c("gene", "isoform")) {
    module_assignments <- data.frame(
        Feature = names(net$colors)
    ) %>%
        mutate(
            module_label = net$colors
        ) %>%
        mutate(
            module_colour = WGCNA::labels2colors(
                module_label, zeroIsGrey = TRUE
            )
        ) %>%
        left_join(
            signedKME(t(expr), net$MEs) %>%
                as.data.frame() %>%
                mutate(Feature = rownames(.)),
            by = "Feature"
        )
}

# for (
#     gene_fn in list.files(
#         "data/genes/Networks", pattern = ".rds", full.names = TRUE
#     )
# ) {
#     write_tsv(
#         process_network(
#             net = readRDS(gene_fn),
#             expr = gene_expr,
#             data_type = "gene"
#         ),
#         file.path(
#             "data/genes/Networks",
#             gsub(".rds", "_ModuleAssign.tsv", basename(gene_fn))
#         )
#     )
# }
# for (
#     iso_fn in list.files(
#         "data/isoforms/Networks", pattern = ".rds", full.names = TRUE
#     )
# ) {
#     write_tsv(
#         process_network(
#             net = readRDS(iso_fn),
#             expr = isoform_expr,
#             data_type = "isoform"
#         ),
#         file.path(
#             "data/isoforms/Networks",
#             gsub(".rds", "_ModuleAssign.tsv", basename(iso_fn))
#         )
#     )
# }

gene_assign_fns <- list.files(
    "data/genes/Networks", pattern = ".tsv", full.names = TRUE
)
iso_assign_fns <- list.files(
    "data/isoforms/Networks", pattern = ".tsv", full.names = TRUE
)

################################################################################
# Hubs                                                                         #
################################################################################

# gene_hubs <- lapply(
#     setNames(
#         nm = gene_assign_fns
#     ),
#     function(fn) {
#         ma <- read_tsv(fn)
#         chooseTopHubInEachModule(
#             datExpr = t(gene_expr),
#             colorh = ma$module_label[match(rownames(gene_expr), ma$Feature)],
#             power = readRDS("data/genes/SoftThresholding.rds")$powerEstimate,
#             type = "signed"
#         )
#     }
# )
# iso_hubs <- lapply(
#     setNames(
#         nm = iso_assign_fns
#     ),
#     function(fn) {
#         ma <- read_tsv(fn)
#         chooseTopHubInEachModule(
#             datExpr = t(isoform_expr),
#             colorh = ma$module_label[match(rownames(isoform_expr), ma$Feature)],
#             power = readRDS("data/isoforms/SoftThresholding.rds")$powerEstimate,
#             type = "signed"
#         )
#     }
# )

################################################################################
# Functional Enrichments                                                       #
################################################################################

dir.create("data/genes/Networks/Enrichments")
dir.create("data/isoforms/Networks/Enrichments")
foreach(i = 1:length(gene_assign_fns)) %dopar% {
    module_assignments <- read_tsv(gene_assign_fns[[i]])
    module_labels <- sort(unique(module_assignments$module_label))
    enr <- lapply(
        setNames(nm = module_labels),
        function(ml) {
            features <- module_assignments %>%
                filter(module_label == ml) %>%
                arrange(desc(!!sym(paste0("kME", ml)))) %>%
                pull(Feature)
            gProfileR::gprofiler(
                query = features,
                organism = "hsapiens", 
                ordered_query = TRUE, 
                exclude_iea = TRUE,
                correction_method = "fdr", 
                hier_filtering = "none", 
                custom_bg = module_assignments$Feature,
                src_filter = c("GO:BP", "GO:MF")
            )
        }
    )
    saveRDS(
        enr,
        file.path(
            "data/genes/Networks/Enrichments",
            gsub(
                "ModuleAssign.tsv", "GO.rds", basename(gene_assign_fns[[i]])
            )
        )
    )
    writexl::write_xlsx(
        enr %>%
            lapply(
                .,
                function(e) {
                    dplyr::select(e, -intersection)
                }
            ),
        file.path(
            "data/genes/Networks/Enrichments",
            gsub(
                "ModuleAssign.tsv", "GO.xlsx", basename(gene_assign_fns[[i]])
            )
        )
    )
}
foreach(i = 1:length(iso_assign_fns)) %dopar% {
    module_assignments <- read_tsv(iso_assign_fns[[i]])
    module_labels <- sort(unique(module_assignments$module_label))
    enr <- lapply(
        setNames(nm = module_labels),
        function(ml) {
            features <- module_assignments %>%
                filter(module_label == ml) %>%
                arrange(desc(!!sym(paste0("kME", ml)))) %>%
                pull(Feature)
            gProfileR::gprofiler(
                query = features,
                organism = "hsapiens", 
                ordered_query = TRUE, 
                exclude_iea = TRUE,
                correction_method = "fdr", 
                hier_filtering = "none", 
                custom_bg = module_assignments$Feature,
                src_filter = c("GO:BP", "GO:MF")
            )
        }
    )
    saveRDS(
        enr,
        file.path(
            "data/isoforms/Networks/Enrichments",
            gsub(
                "ModuleAssign.tsv", "GO.rds", basename(gene_assign_fns[[i]])
            )
        )
    )
    writexl::write_xlsx(
        enr %>%
            lapply(
                .,
                function(e) {
                    dplyr::select(e, -intersection)
                }
            ),
        file.path(
            "data/isoforms/Networks/Enrichments",
            gsub(
                "ModuleAssign.tsv", "GO.xlsx", basename(iso_assign_fns[[i]])
            )
        )
    )
}

# lapply(
#     list.files(
#         "data/genes/Networks", pattern = ".tsv", full.names = TRUE
#     ),
#     function(fn) {
#         data.frame(
#             filename = basename(fn),
#             module_count = length(unique(read_tsv(fn)$module_colour))
#         )
#     }
# ) %>% 
#     bind_rows() %>%
#     View()
# lapply(
#     list.files(
#         "data/isoforms/Networks", pattern = ".tsv", full.names = TRUE
#     ),
#     function(fn) {
#         data.frame(
#             filename = basename(fn),
#             module_count = length(unique(read_tsv(fn)$module_colour))
#         )
#     }
# ) %>% 
#     bind_rows() %>%
#     View()
