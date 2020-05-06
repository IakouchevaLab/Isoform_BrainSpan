library(tidyverse)
library(biomaRt)
source("scripts/utility/plotting.R")

dir.create("data/genes/Enrichments")
dir.create("data/isoforms/Enrichments")

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
tt_genes <- readRDS(
    paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_gene_id")
tt_iso <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")

################################################################################
# Functional Enrichment                                                        #
################################################################################

# go_enrichment <- setNames(nm = c("Gene", "Isoform")) %>%
#     lapply(
#         .,
#         function(data_type) {
#             if (data_type == "Gene") {
#                 tt <- tt_genes
#                 feature <- "ensembl_gene_id"
#             } else {
#                 tt <- tt_iso
#                 feature <- "ensembl_transcript_id"
#             }
#             setNames(nm = sort(unique(tt$Contrast))) %>%
#                 lapply(
#                     .,
#                     function(ctr) {
#                         gprofiler2::gost(
#                             query = tt %>%
#                                 filter(Contrast == ctr) %>%
#                                 filter(
#                                     adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)
#                                 ) %>%
#                                 arrange(adj.P.Val) %>%
#                                 pull(!!sym(feature)),
#                             organism = "hsapiens",
#                             ordered_query = TRUE,
#                             exclude_iea = TRUE,
#                             correction_method = "fdr",
#                             # hier_filtering = "none",
#                             custom_bg = tt[[feature]],
#                             # src_filter = c("GO:BP", "GO:MF")
#                             sources = c("GO:BP", "GO:MF")
#                         )$result
#                     }
#                 )
#         }
#     )
# 
# saveRDS(
#     go_enrichment$Gene,
#     "data/genes/Enrichments/DE_GOST.rds"
# )
# saveRDS(
#     go_enrichment$Isoform,
#     "data/isoforms/Enrichments/DE_GOST.rds"
# )
# 
# go_enrichment <- list(
#     Gene = readRDS("data/genes/Enrichments/DE_GOST.rds"),
#     Isoform = readRDS("data/isoforms/Enrichments/DE_GOST.rds")
# )
# 
# writexl::write_xlsx(
#     go_enrichment$Gene[!sapply(go_enrichment$Gene, is.null)],
#     "data/genes/Enrichments/DE_GOST.xlsx"
# )
# writexl::write_xlsx(
#     go_enrichment$Isoform[!sapply(go_enrichment$Isoform, is.null)],
#     "data/isoforms/Enrichments/DE_GOST.xlsx"
# )
# 
de_genes <- tt_genes %>%
    filter(
        adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)
    )
de_iso <- tt_iso %>%
    filter(
        adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)
    )
# 
# go_enrichment_spec <- setNames(nm = c("Gene", "Isoform")) %>%
#     lapply(
#         .,
#         function(data_type) {
#             if (data_type == "Gene") {
#                 tt <- tt_genes
#                 de <- de_genes
#                 other_de <- de_iso
#                 feature <- "ensembl_gene_id"
#             } else {
#                 tt <- tt_iso
#                 de <- de_iso
#                 other_de <- de_genes
#                 feature <- "ensembl_transcript_id"
#             }
#             setNames(nm = sort(unique(tt$Contrast))) %>%
#                 lapply(
#                     .,
#                     function(ctr) {
#                         deselector <- other_de %>%
#                             filter(Contrast == ctr) %>%
#                             pull(ensembl_gene_id)
#                         gprof_input <- de %>%
#                             filter(Contrast == ctr) %>%
#                             filter(! ensembl_gene_id %in% deselector) %>%
#                             arrange(adj.P.Val) %>%
#                             pull(!!sym(feature))
#                         gprofiler2::gost(
#                             query = gprof_input,
#                             organism = "hsapiens", 
#                             ordered_query = TRUE, 
#                             exclude_iea = TRUE, 
#                             correction_method = "fdr", 
#                             # hier_filtering = "none",
#                             custom_bg = tt[[feature]], 
#                             # src_filter = c("GO:BP", "GO:MF")
#                             sources = c("GO:BP", "GO:MF")
#                         )$result
#                     }
#                 )
#         }
#     )
# 
# saveRDS(
#     go_enrichment_spec$Gene,
#     "data/genes/Enrichments/DE_GOST_SPEC.rds"
# )
# saveRDS(
#     go_enrichment_spec$Isoform,
#     "data/isoforms/Enrichments/DE_GOST_SPEC.rds"
# )
# 
# go_enrichment_spec <- list(
#     Gene = readRDS("data/genes/Enrichments/DE_GOST_SPEC.rds"),
#     Isoform = readRDS("data/isoforms/Enrichments/DE_GOST_SPEC.rds")
# )
# 
# writexl::write_xlsx(
#     go_enrichment_spec$Gene[!sapply(go_enrichment_spec$Gene, is.null)],
#     "data/genes/Enrichments/DE_GOST_SPEC.xlsx"
# )
# writexl::write_xlsx(
#     go_enrichment_spec$Isoform[!sapply(go_enrichment_spec$Isoform, is.null)],
#     "data/isoforms/Enrichments/DE_GOST_SPEC.xlsx"
# )

