library(tidyverse)
library(biomaRt)

mart <- useMart(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt",
    header = TRUE, sep = ",", row.names = 1
)

genes <- rownames(readRDS("data/RegressGeneCounts.rds"))
isoforms <- rownames(readRDS("data/RegressIsoformCounts.rds"))

biotypes <- annotations %>%
    filter(ensembl_transcript_id %in% isoforms) %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "gene_biotype"),
            filters = c("ensembl_gene_id"),
            values = unique(.$ensembl_gene_id),
            mart = mart
        ),
        by = "ensembl_gene_id"
    ) %>%
    dplyr::select(
        ensembl_gene_id, ensembl_transcript_id, 
        external_gene_id, external_transcript_id, 
        gene_biotype, transcript_biotype
    )

write_csv(biotypes, "data/biotypes.csv")
