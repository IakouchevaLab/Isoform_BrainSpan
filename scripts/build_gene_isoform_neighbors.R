library(tidyverse)
library(biomaRt)

mart <- useMart(biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl",
                host = "GRCh37.ensembl.org")

# Build neighbor set for filtered isoforms and genes

filtered_isoforms <- rownames(readRDS("data/iso_tpm_filter.rds"))
filtered_genes <- rownames(readRDS("data/gene_tpm_filter.rds"))

filtered_isoforms_annotations <- getBM(attributes = c("ensembl_gene_id",
                                                      "ensembl_transcript_id",
                                                      "percentage_gene_gc_content",
                                                      "transcript_start",
                                                      "transcript_end"),
                                       filters = c("ensembl_transcript_id"),
                                       mart = mart,
                                       values = filtered_isoforms) %>%
    mutate(transcript_length = abs(transcript_start - transcript_end))
filtered_genes_annotations <- getBM(attributes = c("ensembl_gene_id",
                                                   "percentage_gene_gc_content",
                                                   "start_position",
                                                   "end_position"),
                                       filters = c("ensembl_gene_id"),
                                       mart = mart,
                                       values = filtered_genes) %>%
    mutate(gene_length = abs(start_position - end_position))

counter <- 0
isoform_neighbors <- lapply(
    setNames(nm=pull(filtered_isoforms_annotations, ensembl_transcript_id)),
    function(x) {
        counter <<- counter + 1
        if (counter %% 10000 == 0) { 
            message(Sys.time())
            message(counter, ' ', x) 
        }
        x_gc <- filtered_isoforms_annotations %>%
            filter(ensembl_transcript_id == x) %>%
            pull(percentage_gene_gc_content) %>%
            as.numeric() %>% unique()
        x_ln <- filtered_isoforms_annotations %>%
            filter(ensembl_transcript_id == x) %>%
            pull(transcript_length) %>%
            as.numeric() %>% unique()
        filtered_isoforms_annotations  %>%
            filter(between(percentage_gene_gc_content, x_gc * 0.9, x_gc * 1.1)) %>%
            filter(between(transcript_length, x_ln * 0.9, x_ln * 1.1)) %>%
            pull(ensembl_transcript_id)
    }
)
saveRDS(isoform_neighbors, "data/filtered_isoform_neighbors")
counter <- 0
gene_neighbors <- lapply(
    setNames(nm=pull(filtered_genes_annotations, ensembl_gene_id)),
    function(x) {
        counter <<- counter + 1
        if (counter %% 10000 == 0) { 
            message(Sys.time())
            message(counter, ' ', x) 
        }
        x_gc <- filtered_genes_annotations %>%
            filter(ensembl_gene_id == x) %>%
            pull(percentage_gene_gc_content) %>%
            as.numeric() %>% unique()
        x_ln <- filtered_genes_annotations %>%
            filter(ensembl_gene_id == x) %>%
            pull(gene_length) %>%
            as.numeric() %>% unique()
        filtered_genes_annotations  %>%
            filter(between(percentage_gene_gc_content, x_gc * 0.9, x_gc * 1.1)) %>%
            filter(between(gene_length, x_ln * 0.9, x_ln * 1.1)) %>%
            pull(ensembl_gene_id)
    }
)
saveRDS(gene_neighbors, "data/filtered_gene_neighbors")