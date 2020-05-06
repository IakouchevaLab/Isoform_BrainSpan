# Functional enrichment and curated list enrichments

library(tidyverse)
library(biomaRt)
library(doParallel)

print(getwd())

dir.create("data/genes/Enrichments")
dir.create("data/isoforms/Enrichments")

registerDoParallel(cores = detectCores() - 1)

metadata <- read_tsv("data/metadata.tsv")
anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]

# Add extra information to annotation table
ensembl <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)
anno_expand <- anno %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "start_position", "end_position"),
            filters = "ensembl_gene_id",
            values = unique(.$ensembl_gene_id),
            mart = ensembl
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(gene_length = abs(start_position - end_position))

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

sv_gn <- 13
sv_tx <- 16

################################################################################
# Gene ASD List Enrichment                                                     #
################################################################################
print("Gene permutation analysis")
tt_genes <- readRDS(
    paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_gene_id")
# gene_neighbors <- mclapply(
#     setNames(nm = sort(unique(tt_genes$ensembl_gene_id))),
#     function(g) {
#         if (! g %in% anno_expand$ensembl_gene_id) {
#             return(g)
#         }
#         len <- anno_expand %>%
#             filter(ensembl_gene_id == g) %>%
#             distinct(gene_length) %>%
#             pull(gene_length) %>%
#             max()
#         gc <- anno_expand %>%
#             filter(ensembl_gene_id == g) %>%
#             distinct(percentage_gc_content) %>%
#             pull(percentage_gc_content) %>%
#             max()
#         anno_expand %>%
#             filter(
#                 gene_length >= 0.9 * len & gene_length >= 1.1 * len
#             ) %>%
#             filter(
#                 percentage_gc_content >= 0.9 * gc &
#                     percentage_gc_content >= 1.1 * gc
#             ) %>%
#             distinct(ensembl_gene_id) %>%
#             pull(ensembl_gene_id)
#     }
# )
# saveRDS(
#     gene_neighbors,
#     "data/genes/FeatureNeighbors_10pLengthGC.rds"
# )
# gene_neighbors <- readRDS("data/genes/FeatureNeighbors_10pLengthGC.rds")
# empirical_overlaps <- expand.grid(
#     List = names(asd_lists),
#     Contrast = sort(unique(tt_genes$Contrast)),
#     Direction = c(1, -1, 0)
# ) %>%
#     mutate(
#         overlap = apply(
#             ., 1, function(eo) {
#                 list_genes <- anno_expand$ensembl_gene_id[
#                     match(
#                         asd_lists[[as.character(eo[["List"]])]],
#                         anno_expand$external_gene_id
#                     )
#                 ]
#                 list_genes <- list_genes[!is.na(list_genes)]
#                 de_genes <- tt_genes %>%
#                     filter(Contrast == as.character(eo[["Contrast"]])) %>%
#                     filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
#                     filter(
#                         ifelse(
#                             as.numeric(eo[["Direction"]]) == 0,
#                             TRUE,
#                             sign(logFC) == as.numeric(eo[["Direction"]])
#                         )
#                     ) %>%
#                     pull(ensembl_gene_id)
#                 length(intersect(list_genes, de_genes))
#             }
#         )
#     )
# null_overlaps <- expand.grid(
#     Iteration = seq(1, 1000),
#     List = names(asd_lists),
#     Contrast = sort(unique(tt_genes$Contrast)),
#     Direction = c(1, -1, 0)
# )
# null_overlaps$overlap <- foreach(i = 1:nrow(null_overlaps)) %dopar% {
#     list_genes <- anno_expand$ensembl_gene_id[
#         match(
#             asd_lists[[as.character(null_overlaps[i, "List"])]],
#             anno_expand$external_gene_id
#         )
#         ]
#     list_genes <- list_genes[!is.na(list_genes)]
#     de_genes <- tt_genes %>%
#         filter(Contrast == as.character(null_overlaps[i, "Contrast"])) %>%
#         filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5))
#     if (as.numeric(null_overlaps[i, "Direction"]) != 0) {
#         de_genes <- de_genes %>%
#             filter(sign(logFC) == as.numeric(null_overlaps[i, "Direction"]))
#     }
#     de_genes <- pull(de_genes, ensembl_gene_id)
#     random_genes <- sapply(de_genes, function(gn) {
#         tryCatch({
#             sample(gene_neighbors[[gn]], 1)
#         }, error = function(e) {
#             return(gn)
#         })
#     })
#     random_genes_symbols <- anno_expand$external_gene_id[
#         match(random_genes, anno_expand$ensembl_gene_id)
#     ]
#     length(
#         intersect(
#             list_genes,
#             random_genes_symbols[!is.na(random_genes_symbols)]
#         )
#     )
# }
# saveRDS(
#     list(empirical = empirical_overlaps, null = null_overlaps),
#     "data/genes/Enrichments/ASDListPermutationTests.rds"
# )
# 
# permutation_tests_genes_results <- readRDS(
#     "data/genes/Enrichments/ASDListPermutationTests.rds"
# )

################################################################################
# Isoform ASD List Enrichment                                                  #
################################################################################
print("Isoform permutation analysis")
tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")
# iso_neighbors <- mclapply(
#     setNames(nm = sort(unique(tt_iso$ensembl_transcript_id))),
#     function(tx) {
#         if (! tx %in% anno_expand$ensembl_transcript_id) {
#             return(tx)
#         }
#         len <- anno_expand %>%
#             filter(ensembl_transcript_id == tx) %>%
#             distinct(transcript_length) %>%
#             pull(transcript_length) %>%
#             max()
#         gc <- anno_expand %>%
#             filter(ensembl_transcript_id == tx) %>%
#             distinct(percentage_gc_content) %>%
#             pull(percentage_gc_content) %>%
#             max()
#         anno_expand %>%
#             filter(
#                 transcript_length >= 0.9 * len & transcript_length >= 1.1 * len
#             ) %>%
#             filter(
#                 percentage_gc_content >= 0.9 * gc &
#                     percentage_gc_content >= 1.1 * gc
#             ) %>%
#             distinct(ensembl_transcript_id) %>%
#             pull(ensembl_transcript_id)
#     }
# )
# saveRDS(
#     iso_neighbors,
#     "data/isoforms/FeatureNeighbors_10pLengthGC.rds"
# )
iso_neighbors <- readRDS("data/isoforms/FeatureNeighbors_10pLengthGC.rds")
empirical_overlaps <- expand.grid(
    List = names(asd_lists),
    Contrast = sort(unique(tt_iso$Contrast)),
    Direction = c(1, -1, 0)
) %>%
    mutate(
        overlap = apply(
            ., 1, function(eo) {
                list_genes <- anno_expand$ensembl_gene_id[
                    match(
                        asd_lists[[as.character(eo[["List"]])]],
                        anno_expand$external_gene_id
                    )
                ]
                list_genes <- list_genes[!is.na(list_genes)]
                de_iso <- tt_iso %>%
                    filter(Contrast == as.character(eo[["Contrast"]])) %>%
                    filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
                    filter(
                        ifelse(
                            as.numeric(eo[["Direction"]]) == 0,
                            TRUE,
                            sign(logFC) == as.numeric(eo[["Direction"]])
                        )
                    ) %>%
                    pull(ensembl_gene_id)
                length(intersect(list_genes, de_iso))
            }
        )
    )
null_overlaps <- expand.grid(
    Iteration = seq(1, 1000),
    List = names(asd_lists),
    Contrast = sort(unique(tt_genes$Contrast)),
    Direction = c(1, -1, 0)
)
null_overlaps$overlap <- foreach(i = 1:nrow(null_overlaps)) %dopar% {
    list_genes <- anno_expand$ensembl_gene_id[
        match(
            asd_lists[[as.character(null_overlaps[i, "List"])]],
            anno_expand$external_gene_id
        )
        ]
    list_genes <- list_genes[!is.na(list_genes)]
    de_iso <- tt_iso %>%
        filter(Contrast == as.character(null_overlaps[i, "Contrast"])) %>%
        filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5))
    if(as.numeric(null_overlaps[i, "Direction"]) != 0) {
        de_iso <- de_iso %>%
            filter(
                sign(logFC) == as.numeric(null_overlaps[i, "Direction"])
            )
    }
    de_iso <- de_iso %>%
        pull(ensembl_transcript_id)
    random_iso <- sapply(de_iso, function(tx) {
        tryCatch({
            sample(iso_neighbors[[tx]], 1)
        }, error = function(e) {
            return(tx)
        })
    })
    random_genes_symbols <- anno_expand$external_gene_id[
        match(random_iso, anno_expand$ensembl_transcript_id)
    ]
    length(intersect(
        list_genes, random_genes_symbols[!is.na(random_genes_symbols)]
    ))
}
saveRDS(
    list(empirical = empirical_overlaps, null = null_overlaps),
    "data/isoforms/Enrichments/ASDListPermutationTests.rds"
)

permutation_tests_isoforms_results <- readRDS(
    "data/isoforms/Enrichments/ASDListPermutationTests.rds"
)