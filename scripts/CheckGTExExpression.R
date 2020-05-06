library(tidyverse)


gtex_cor_perm <- readRDS("data/gtex_permutations.rds")

gtex_expr <- read_tsv("data/source/GTEx_V8/gtex_brainspan_20-40.txt", skip = 1) %>%
    mutate(ensembl_transcript_id = gsub("\\.[0-9]+", "", transcript_id))
itpm <- readRDS("data/iso_tpm_filter.rds")
metadata <- read_tsv("data/metadata.tsv")
p13samples <- metadata$Sample[metadata$Period == 13]
intersect_isoforms <- intersect(rownames(itpm), gtex_expr$ensembl_transcript_id)
gtex_brain_samples <- c(
    "AMY" = "Brain - Amygdala",
    "CBC" = "Brain - Cerebellar Hemisphere",
    "HIP" = "Brain - Hippocampus"
)
amy_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "AMY") %>%
    pull(Sample)
cbc_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "CBC") %>%
    pull(Sample)
hip_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "HIP") %>%
    pull(Sample)
amy_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Amygdala") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))
cbc_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Cerebellar Hemisphere") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))
hip_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Hippocampus") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))

random_check_isoforms <- sample(intersect_isoforms, size = 10, replace = FALSE)

gtex_expr[]