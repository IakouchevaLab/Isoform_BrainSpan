library(tidyverse)
library(readxl)

slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
n_perm <- 10

isoform_neighbors <- readRDS("data/filtered_isoform_neighbors")
gene_neighbors <- readRDS("data/filtered_gene_neighbors")

deg_sheets <- excel_sheets("data/SupplementaryTables/Supplementary Table 2.xlsx")[2:13]
dei_sheets <- excel_sheets("data/SupplementaryTables/Supplementary Table 3.xlsx")[2:13]
gm_sheets <- excel_sheets("data/SupplementaryTables/Supplementary Table 8.xlsx")[2:10]
im_sheets <- excel_sheets("data/SupplementaryTables/Supplementary Table 9.xlsx")[2:57]
list_sheets <- excel_sheets("data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx")[2:12]
cell_types <- readxl::read_xlsx(
    file.path(
        "data/source/CuratedLists",
        "ZhongNature2018_humanPreFrontalCortexCellTypes.xlsx"
    ),
    sheet = 1, skip = 4
)[, c("cluster", "gene")]

asd_lists <- lapply(
    setNames(nm=list_sheets),
    function(s) {
        read_excel("data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx", sheet=s)[[1]]
    }
)

# DE Enrichments

deg_permutation_tests <- list()
# for (deg_sheet in deg_sheets) {
#     message(deg_sheet)
#     deg <- read_excel("data/SupplementaryTables/Supplementary Table 2.xlsx",
#                       sheet = deg_sheet) %>%
#         filter(FDR <= 0.05) %>%
#         filter(abs(logFC) >= log2(1.5))
#     this_deg_empirical_asd <- sapply(
#         setNames(nm=names(asd_lists)),
#         function(asd_list) {
#             length(intersect(pull(deg, `Gene Symbol`), asd_lists[[asd_list]]))
#         }
#     )
#     this_deg_empirical_ct <- sapply(
#         setNames(nm=sort(unique(cell_types$cluster))),
#         function(cell_type) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             length(intersect(pull(deg, `Gene Symbol`), cell_type_list))
#         }
#     )
#     this_deg_permutations_asd <- list()
#     this_deg_permutations_ct <- list()
#     for (i in seq(1, n_perm)) {
#         message(i)
#         random_genes <- sapply(pull(deg, `Ensembl Gene ID`), function(x) {
#             sample(gene_neighbors[[x]], 1)
#         })
#         for (asd_list in names(asd_lists)) {
#             this_deg_permutations_asd[[asd_list]] <- c(this_deg_permutations_asd[[asd_list]], length(intersect(pull(filter(deg, `Ensembl Gene ID` %in% random_genes), `Gene Symbol`), asd_lists[[asd_list]])))
#         }
#         for (cell_type in sort(unique(cell_types$cluster))) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             this_deg_permutations_ct[[cell_type]] <- c(this_deg_permutations_ct[[cell_type]], length(intersect(pull(filter(deg, `Ensembl Gene ID` %in% random_genes), `Gene Symbol`), cell_type_list)))
#         }
#     }
#     deg_permutation_tests[[deg_sheet]] <- list(
#         asd_empirical = this_deg_empirical_asd,
#         asd_permutation = this_deg_permutations_asd,
#         ct_empirical = this_deg_empirical_ct,
#         ct_permutation = this_deg_permutations_ct
#     )
# }
# saveRDS(deg_permutation_tests,
#         paste0("data/list_permutation_results/deg_permutations", slurm_task, ".rds"))

deg_initial <- readRDS("data/list_permutation_results/deg_permutations.rds")
deg_initial_tbl <- lapply(names(deg_initial), function(p) {
    asd_empirical <- deg_initial[[p]][["asd_empirical"]]
    ct_empirical <- deg_initial[[p]][["ct_empirical"]]
    bind_rows(
        tibble(
            contrast = p,
            list = names(asd_empirical),
            empirical = asd_empirical
        ),
        tibble(
            contrast = p,
            list = names(ct_empirical),
            empirical = ct_empirical
        )
    )
}) %>%
    bind_rows()
deg_permutations <- list()
for (dat in list.files("data/list_permutation_results/", pattern="deg_permutations\\d+.rds", full.names=TRUE)) {
        l <- readRDS(dat)
        for (my_contrast in names(l)) {
            asd_permutations <- l[[my_contrast]][["asd_permutation"]]
            ct_permutations <- l[[my_contrast]][["ct_permutation"]]
            for (asd_list in names(asd_permutations)) {
                deg_permutations[[my_contrast]][[asd_list]] <- c(
                    unlist(deg_permutations[[my_contrast]][[asd_list]]),
                    unlist(asd_permutations[[asd_list]], use.names=FALSE)
                )
            }
            for (ct_list in names(ct_permutations)) {
                deg_permutations[[my_contrast]][[ct_list]] <- c(
                    unlist(deg_permutations[[my_contrast]][[ct_list]]),
                    unlist(ct_permutations[[ct_list]], use.names=FALSE)
                )
            }
        }
}
deg_pval <- deg_initial_tbl %>%
    mutate(
        pval = apply(., 1, function(x) {
            my_contrast <- x[["contrast"]]
            my_list <- x[["list"]]
            my_empirical <- as.numeric(x[["empirical"]])
            1 - (sum(my_empirical > deg_permutations[[my_contrast]][[my_list]]) / (length(deg_permutations[[my_contrast]][[my_list]]) + 1))
        })
    ) %>%
    mutate(
        padj = p.adjust(pval, method = "fdr")
    )

# dei_permutation_tests <- list()
# for (dei_sheet in dei_sheets) {
#     message(dei_sheet)
#     dei <- read_excel("data/SupplementaryTables/Supplementary Table 3.xlsx",
#                       sheet = dei_sheet) %>%
#         filter(FDR <= 0.05) %>%
#         filter(abs(logFC) >= log2(1.5))
#     this_dei_empirical_asd <- sapply(
#         setNames(nm=names(asd_lists)),
#         function(asd_list) {
#             length(intersect(pull(dei, `Gene Symbol`), asd_lists[[asd_list]]))
#         }
#     )
#     this_dei_empirical_ct <- sapply(
#         setNames(nm=sort(unique(cell_types$cluster))),
#         function(cell_type) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             length(intersect(pull(dei, `Gene Symbol`), cell_type_list))
#         }
#     )
#     this_dei_permutations_asd <- list()
#     this_dei_permutations_ct <- list()
#     for (i in seq(1, n_perm)) {
#         message(i)
#         random_isoforms <- sapply(pull(dei, `Ensembl Transcript ID`), function(x) {
#             sample(isoform_neighbors[[x]], 1)
#         })
#         for (asd_list in names(asd_lists)) {
#             this_dei_permutations_asd[[asd_list]] <- c(this_dei_permutations_asd[[asd_list]], length(intersect(pull(filter(dei, `Ensembl Transcript ID` %in% random_isoforms), `Gene Symbol`), asd_lists[[asd_list]])))
#         }
#         for (cell_type in sort(unique(cell_types$cluster))) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             this_dei_permutations_ct[[cell_type]] <- c(this_dei_permutations_ct[[cell_type]], length(intersect(pull(filter(dei, `Ensembl Transcript ID` %in% random_isoforms), `Gene Symbol`), cell_type_list)))
#         }
#     }
#     dei_permutation_tests[[dei_sheet]] <- list(
#         asd_empirical = this_dei_empirical_asd,
#         asd_permutation = this_dei_permutations_asd,
#         ct_empirical = this_dei_empirical_ct,
#         ct_permutation = this_dei_permutations_ct
#     )
# }
# saveRDS(dei_permutation_tests,
#         paste0("data/list_permutation_results/dei_permutations", slurm_task, ".rds"))

dei_initial <- readRDS("data/list_permutation_results/dei_permutations.rds")
dei_initial_tbl <- lapply(names(dei_initial), function(p) {
    asd_empirical <- dei_initial[[p]][["asd_empirical"]]
    ct_empirical <- dei_initial[[p]][["ct_empirical"]]
    bind_rows(
        tibble(
            contrast = p,
            list = names(asd_empirical),
            empirical = asd_empirical
        ),
        tibble(
            contrast = p,
            list = names(ct_empirical),
            empirical = ct_empirical
        )
    )
}) %>%
    bind_rows()
dei_permutations <- list()
for (dat in list.files("data/list_permutation_results/", pattern="dei_permutations\\d+.rds", full.names=TRUE)) {
    l <- readRDS(dat)
    for (my_contrast in names(l)) {
        asd_permutations <- l[[my_contrast]][["asd_permutation"]]
        ct_permutations <- l[[my_contrast]][["ct_permutation"]]
        for (asd_list in names(asd_permutations)) {
            dei_permutations[[my_contrast]][[asd_list]] <- c(
                unlist(dei_permutations[[my_contrast]][[asd_list]]),
                unlist(asd_permutations[[asd_list]], use.names=FALSE)
            )
        }
        for (ct_list in names(ct_permutations)) {
            dei_permutations[[my_contrast]][[ct_list]] <- c(
                unlist(dei_permutations[[my_contrast]][[ct_list]]),
                unlist(ct_permutations[[ct_list]], use.names=FALSE)
            )
        }
    }
}
dei_pval <- dei_initial_tbl %>%
    mutate(
        pval = apply(., 1, function(x) {
            my_contrast <- x[["contrast"]]
            my_list <- x[["list"]]
            my_empirical <- as.numeric(x[["empirical"]])
            1 - (sum(my_empirical > dei_permutations[[my_contrast]][[my_list]]) / (length(dei_permutations[[my_contrast]][[my_list]]) + 1))
        })
    ) %>%
    mutate(
        padj = p.adjust(pval, method = "fdr")
    )

# Module enrichments

# gm_permutation_tests <- list()
# for (gm_sheet in gm_sheets) {
#     message(gm_sheet)
#     gm <- read_excel("data/SupplementaryTables/Supplementary Table 8.xlsx",
#                       sheet = gm_sheet)
#     this_gm_empirical_asd <- sapply(
#         setNames(nm=names(asd_lists)),
#         function(asd_list) {
#             length(intersect(pull(gm, `Gene Symbol`), asd_lists[[asd_list]]))
#         }
#     )
#     this_gm_empirical_ct <- sapply(
#         setNames(nm=sort(unique(cell_types$cluster))),
#         function(cell_type) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             length(intersect(pull(gm, `Gene Symbol`), cell_type_list))
#         }
#     )
#     this_gm_permutations_asd <- list()
#     this_gm_permutations_ct <- list()
#     for (i in seq(1, n_perm)) {
#         message(i)
#         random_genes <- sapply(pull(gm, `Ensembl Gene ID`), function(x) {
#             sample(gene_neighbors[[x]], 1)
#         })
#         for (asd_list in names(asd_lists)) {
#             this_gm_permutations_asd[[asd_list]] <- c(this_gm_permutations_asd[[asd_list]], length(intersect(pull(filter(gm, `Ensembl Gene ID` %in% random_genes), `Gene Symbol`), asd_lists[[asd_list]])))
#         }
#         for (cell_type in sort(unique(cell_types$cluster))) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             this_gm_permutations_ct[[cell_type]] <- c(this_gm_permutations_ct[[cell_type]], length(intersect(pull(filter(gm, `Ensembl Gene ID` %in% random_genes), `Gene Symbol`), cell_type_list)))
#         }
#     }
#     gm_permutation_tests[[gm_sheet]] <- list(
#         asd_empirical = this_gm_empirical_asd,
#         asd_permutation = this_gm_permutations_asd,
#         ct_empirical = this_gm_empirical_ct,
#         ct_permutation = this_gm_permutations_ct
#     )
# }
# saveRDS(gm_permutation_tests,
#         paste0("data/list_permutation_results/gm_permutations", slurm_task, ".rds"))

gm_initial <- readRDS("data/list_permutation_results/gm_permutations.rds")
gm_initial_tbl <- lapply(names(gm_initial), function(p) {
    asd_empirical <- gm_initial[[p]][["asd_empirical"]]
    ct_empirical <- gm_initial[[p]][["ct_empirical"]]
    bind_rows(
        tibble(
            contrast = p,
            list = names(asd_empirical),
            empirical = asd_empirical
        ),
        tibble(
            contrast = p,
            list = names(ct_empirical),
            empirical = ct_empirical
        )
    )
}) %>%
    bind_rows()
gm_permutations <- list()
for (dat in list.files("data/list_permutation_results/", pattern="gm_permutations\\d+.rds", full.names=TRUE)) {
    l <- readRDS(dat)
    for (my_contrast in names(l)) {
        asd_permutations <- l[[my_contrast]][["asd_permutation"]]
        ct_permutations <- l[[my_contrast]][["ct_permutation"]]
        for (asd_list in names(asd_permutations)) {
            gm_permutations[[my_contrast]][[asd_list]] <- c(
                unlist(gm_permutations[[my_contrast]][[asd_list]]),
                unlist(asd_permutations[[asd_list]], use.names=FALSE)
            )
        }
        for (ct_list in names(ct_permutations)) {
            gm_permutations[[my_contrast]][[ct_list]] <- c(
                unlist(gm_permutations[[my_contrast]][[ct_list]]),
                unlist(ct_permutations[[ct_list]], use.names=FALSE)
            )
        }
    }
}
gm_pval <- gm_initial_tbl %>%
    mutate(
        pval = apply(., 1, function(x) {
            my_contrast <- x[["contrast"]]
            my_list <- x[["list"]]
            my_empirical <- as.numeric(x[["empirical"]])
            1 - (sum(my_empirical > gm_permutations[[my_contrast]][[my_list]]) / (length(gm_permutations[[my_contrast]][[my_list]]) + 1))
        })
    ) %>%
    mutate(
        padj = p.adjust(pval, method = "fdr")
    )

# im_permutation_tests <- list()
# for (im_sheet in im_sheets) {
#     message(im_sheet)
#     im <- read_excel("data/SupplementaryTables/Supplementary Table 9.xlsx",
#                      sheet = im_sheet)
#     this_im_empirical_asd <- sapply(
#         setNames(nm=names(asd_lists)),
#         function(asd_list) {
#             length(intersect(pull(im, `Gene Symbol`), asd_lists[[asd_list]]))
#         }
#     )
#     this_im_empirical_ct <- sapply(
#         setNames(nm=sort(unique(cell_types$cluster))),
#         function(cell_type) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             length(intersect(pull(im, `Gene Symbol`), cell_type_list))
#         }
#     )
#     this_im_permutations_asd <- list()
#     this_im_permutations_ct <- list()
#     for (i in seq(1, n_perm)) {
#         message(i)
#         random_isoforms <- sapply(pull(im, `Ensembl Transcript ID`), function(x) {
#             sample(isoform_neighbors[[x]], 1)
#         })
#         for (asd_list in names(asd_lists)) {
#             this_im_permutations_asd[[asd_list]] <- c(this_im_permutations_asd[[asd_list]], length(intersect(pull(filter(im, `Ensembl Transcript ID` %in% random_isoforms), `Gene Symbol`), asd_lists[[asd_list]])))
#         }
#         for (cell_type in sort(unique(cell_types$cluster))) {
#             cell_type_list <- cell_types %>% 
#                 filter(cluster == cell_type) %>% pull(gene)
#             this_im_permutations_ct[[cell_type]] <- c(this_im_permutations_ct[[cell_type]], length(intersect(pull(filter(im, `Ensembl Transcript ID` %in% random_isoforms), `Gene Symbol`), cell_type_list)))
#         }
#     }
#     im_permutation_tests[[im_sheet]] <- list(
#         asd_empirical = this_im_empirical_asd,
#         asd_permutation = this_im_permutations_asd,
#         ct_empirical = this_im_empirical_ct,
#         ct_permutation = this_im_permutations_ct
#     )
# }
# saveRDS(im_permutation_tests,
#         paste0("data/list_permutation_results/im_permutations", slurm_task, ".rds"))

im_initial <- readRDS("data/list_permutation_results/im_permutations.rds")
im_initial_tbl <- lapply(names(im_initial), function(p) {
    asd_empirical <- im_initial[[p]][["asd_empirical"]]
    ct_empirical <- im_initial[[p]][["ct_empirical"]]
    bind_rows(
        tibble(
            contrast = p,
            list = names(asd_empirical),
            empirical = asd_empirical
        ),
        tibble(
            contrast = p,
            list = names(ct_empirical),
            empirical = ct_empirical
        )
    )
}) %>%
    bind_rows()
im_permutations <- list()
for (dat in list.files("data/list_permutation_results/", pattern="im_permutations\\d+.rds", full.names=TRUE)) {
    l <- readRDS(dat)
    for (my_contrast in names(l)) {
        asd_permutations <- l[[my_contrast]][["asd_permutation"]]
        ct_permutations <- l[[my_contrast]][["ct_permutation"]]
        for (asd_list in names(asd_permutations)) {
            im_permutations[[my_contrast]][[asd_list]] <- c(
                unlist(im_permutations[[my_contrast]][[asd_list]]),
                unlist(asd_permutations[[asd_list]], use.names=FALSE)
            )
        }
        for (ct_list in names(ct_permutations)) {
            im_permutations[[my_contrast]][[ct_list]] <- c(
                unlist(im_permutations[[my_contrast]][[ct_list]]),
                unlist(ct_permutations[[ct_list]], use.names=FALSE)
            )
        }
    }
}
im_pval <- im_initial_tbl %>%
    mutate(
        pval = apply(., 1, function(x) {
            my_contrast <- x[["contrast"]]
            my_list <- x[["list"]]
            my_empirical <- as.numeric(x[["empirical"]])
            1 - (sum(my_empirical > im_permutations[[my_contrast]][[my_list]]) / (length(im_permutations[[my_contrast]][[my_list]]) + 1))
        })
    ) %>%
    mutate(
        padj = p.adjust(pval, method = "fdr")
    )


pval_data <- bind_rows(
    mutate(deg_pval, data = 'DEG'),
    mutate(dei_pval, data = 'DEI'),
    mutate(gm_pval, data = 'GM'),
    mutate(im_pval, data = 'IM'),
)

write_csv(pval_data, "data/permutation_list_test_de_coex.csv")

all_deg <- lapply(
    setNames(nm=excel_sheets("data/SupplementaryTables/Supplementary Table 2.xlsx")[2:13]),
    function(s) {
        read_excel("data/SupplementaryTables/Supplementary Table 2.xlsx", sheet=s) %>%
            filter(FDR <= 0.05) %>%
            filter(abs(logFC) >= log2(1.5)) %>%
            mutate(contrast = s) %>%
            mutate(data = "DEG")
    }
) %>%
    bind_rows()
all_dei <- lapply(
    setNames(nm=excel_sheets("data/SupplementaryTables/Supplementary Table 3.xlsx")[2:13]),
    function(s) {
        read_excel("data/SupplementaryTables/Supplementary Table 3.xlsx", sheet=s) %>%
            filter(FDR <= 0.05) %>%
            filter(abs(logFC) >= log2(1.5)) %>%
            mutate(contrast = s) %>%
            mutate(data = "DEI")
    }
) %>%
    bind_rows()
all_gm <- lapply(
    setNames(nm=excel_sheets("data/SupplementaryTables/Supplementary Table 8.xlsx")[2:10]),
    function(s) {
        read_excel("data/SupplementaryTables/Supplementary Table 8.xlsx", sheet=s) %>%
            mutate(contrast = s) %>%
            mutate(data = "GM")
    }
) %>%
    bind_rows()
all_im <- lapply(
    setNames(nm=excel_sheets("data/SupplementaryTables/Supplementary Table 9.xlsx")[2:57]),
    function(s) {
        read_excel("data/SupplementaryTables/Supplementary Table 9.xlsx", sheet=s) %>%
            mutate(contrast = s) %>%
            mutate(data = "IM")
    }
) %>%
    bind_rows()
all_data <- bind_rows(all_deg, all_dei, all_gm, all_im)
pval_data <- read_csv("data/permutation_list_test_de_coex.csv") %>%
    mutate(
        jaccard = mapply(
            function(i, j) {
                contrast_genes <- all_data %>%
                    filter(contrast == i) %>%
                    pull(`Gene Symbol`)
                if (j %in% names(asd_lists)) {
                    list_genes <- asd_lists[[j]]
                } else {
                    list_genes <- cell_types %>% filter(cluster == j) %>% pull(gene)
                }
                return(length(intersect(contrast_genes, list_genes)) / length(union(contrast_genes, list_genes)))
            },
            contrast, list
        )
    ) %>%
    rename(
        gene_list = list
    )

plt <- ggplot(
    data = pval_data %>%
        mutate(padj_all = p.adjust(pval, method = "fdr")) %>%
        mutate(
            padj_star = ifelse(
                padj_all <= 0.05, 
                ifelse(
                    padj_all <= 0.01,
                    ifelse(
                        padj_all <= 0.005, 
                        '***', '**'
                    ), '*'
                ), ''
            )
        ),
    mapping = aes(x = gene_list, y = contrast, fill = jaccard, label = padj_star)
) +
    facet_grid(data ~ ., scales = "free_y") +
    geom_tile() +
    geom_text() +
    scale_fill_continuous(low = "white", high = "red") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

ggsave("data/figures/list_permutation_sig.pdf", plt, width = 8, height = 24)
