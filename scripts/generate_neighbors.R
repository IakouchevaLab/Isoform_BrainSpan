library(tidyverse)
library(biomaRt)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

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

sv_gn <- 13
sv_tx <- 16

tt_genes <- readRDS(
    paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_gene_id")

gene_neighbors <- mclapply(
    setNames(nm = sort(unique(tt_genes$ensembl_gene_id))),
    function(g) {
        if (! g %in% anno_expand$ensembl_gene_id) {
            # return(g)
            return(
                data.frame(
                    source = g,
                    neighbors = g
                )
            )
        }
        len <- anno_expand %>%
            filter(ensembl_gene_id == g) %>%
            distinct(gene_length) %>%
            pull(gene_length) %>%
            max()
        gc <- anno_expand %>%
            filter(ensembl_gene_id == g) %>%
            distinct(percentage_gc_content) %>%
            pull(percentage_gc_content) %>%
            max()
        g_neighbors <- anno_expand %>%
            filter(
                gene_length >= 0.9 * len & gene_length >= 1.1 * len
            ) %>%
            filter(
                percentage_gc_content >= 0.9 * gc &
                    percentage_gc_content >= 1.1 * gc
            ) %>%
            distinct(ensembl_gene_id) %>%
            pull(ensembl_gene_id)
        if (length(g_neighbors) == 0) {
            return(
                data.frame(
                    source = g,
                    neighbors = g
                )
            )
        } else {
            return(
                data.frame(
                    source = g,
                    neighbors = g_neighbors
                )
            )
        }
        
    }
) %>% 
    bind_rows()

saveRDS(
    gene_neighbors,
    "data/genes/FeatureNeighbors_10pLengthGC_asDataFrame.rds"
)

tt_isoforms <- readRDS(
    paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
) %>%
    left_join(anno_expand, by = "ensembl_transcript_id")

isoform_neighbors <- mclapply(
    setNames(nm = sort(unique(tt_isoforms$ensembl_transcript_id))),
    function(g) {
        if (! g %in% anno_expand$ensembl_transcript_id) {
            # return(g)
            return(
                data.frame(
                    source = g,
                    neighbors = g
                )
            )
        }
        len <- anno_expand %>%
            filter(ensembl_transcript_id == g) %>%
            distinct(transcript_length) %>%
            pull(transcript_length) %>%
            max()
        gc <- anno_expand %>%
            filter(ensembl_transcript_id == g) %>%
            distinct(percentage_gc_content) %>%
            pull(percentage_gc_content) %>%
            max()
        g_neighbors <- anno_expand %>%
            filter(
                transcript_length >= 0.9 * len & transcript_length >= 1.1 * len
            ) %>%
            filter(
                percentage_gc_content >= 0.9 * gc &
                    percentage_gc_content >= 1.1 * gc
            ) %>%
            distinct(ensembl_transcript_id) %>%
            pull(ensembl_transcript_id)
        if (length(g_neighbors) == 0) {
            return(
                data.frame(
                    source = g,
                    neighbors = g
                )
            )
        } else {
            return(
                data.frame(
                    source = g,
                    neighbors = g_neighbors
                )
            )
        }
        
    }
) %>%
    bind_rows()

saveRDS(
    isoform_neighbors,
    "data/isoforms/FeatureNeighbors_10pLengthGC_asDataFrame.rds"
)