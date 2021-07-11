library(tidyverse)
library(ggsignif)

metadata <- read_tsv("data/metadata.tsv")
iso_tpm <- readRDS("data/iso_tpm_filter.rds")
gn_tpm <- readRDS("data/gene_tpm_filter.rds")
splice_variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
) %>%
    filter(ifelse(is.na(str_match(Consequence, "splice_region")), TRUE, FALSE)) %>%
    filter(ifelse(is.na(str_match(Consequence, "splice")), FALSE, TRUE))

splice_transcripts <- splice_variants %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Transcript ID`) %>%
    unique()
splice_genes <- splice_variants %>%
    filter(`Affected status` == 1) %>%
    pull(`Ensembl Gene ID`) %>%
    unique()

splice_variants %>%
    distinct(Chromosome, `Variant start`, `Variant end`, Reference, Alternate, `Affected status`) %>%
    count(`Affected status`)
