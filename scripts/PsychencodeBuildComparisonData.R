# Subset and consolidate overlapping samples between our BrainSpan data and the Psychencode cmc data
# Normal controls

library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

# Load sample metadata and map to each other
cmc_metadata <- read.csv('data/psychencode/psychencode_isoform_clean_input/Capstone_datMeta.csv', row.names=1) %>%
    mutate(Sample = rownames(.)) %>%
    filter(diagnosis == 'Control') %>%
    filter(study == 'CMC_HBCC') %>% 
    mutate(RegionMatch = as.character(tissue)) %>%
    mutate(AgeMatch = as.character(cut(ageDeath, seq(-5, 100, by = 5)))) %>%
    mutate(Group = paste0(RegionMatch, AgeMatch))

brainspan_metadata <- read_tsv('data/metadata.tsv') %>%
    mutate(RegionMatch = sapply(Regioncode, function(r) {
        if (r %in% c('OFC', 'DFC', 'VFC', 'MFC', 'M1C')) {
            return('frontal cortex')
        } else if (r %in% c('A1C', 'STC', 'ITC')) {
            return('temporal cortex')
        } else {
            return(NA)
        }
    })) %>%
    mutate(AgeMatch = as.character(cut(Days / 365, seq(-5, 100, by = 5)))) %>%
    mutate(Group = paste0(RegionMatch, AgeMatch))

load('data/psychencode/psychencode_isoform_clean_input/Isoform.TPM.noOutliersamples.AllIsoforms.RData')
colnames(tpm) = str_replace(colnames(tpm), pattern = "\\.", replacement = "-")
cmc_tpm <- tpm[, str_detect(colnames(tpm), "CMC_HBCC.+")]

load('data/psychencode/psychencode_isoform_clean_input/Isoform.Counts.noOutliersamples.AllIsoforms.RData')
colnames(counts) = str_replace(colnames(counts), pattern = "\\.", replacement = "-")
cmc_cts <- counts[, str_detect(colnames(counts), "CMC_HBCC.+")]

brainspan_tpm <- readRDS('data/iso_tpm_filter.rds')
brainspan_cts <- readRDS('data/iso_cts_filter.rds')

# Subsetting for matched data
# TPM
groups <- intersect(cmc_metadata$Group, brainspan_metadata$Group)
overlap_transcripts <- intersect(rownames(cmc_tpm), rownames(brainspan_tpm))
cmc_matched <- lapply(groups, function(g) {
    rowMeans(cmc_tpm[overlap_transcripts, filter(cmc_metadata, Group == g) %>% pull(Sample), drop = FALSE])
}) %>% 
    bind_cols() %>%
    setNames(groups)
rownames(cmc_matched) <- overlap_transcripts
bspan_matched <- lapply(groups, function(g) {
    rowMeans(brainspan_tpm[overlap_transcripts, filter(brainspan_metadata, Group == g) %>% pull(Sample), drop = FALSE])
}) %>% 
    bind_cols() %>%
    setNames(groups)
rownames(bspan_matched) <- overlap_transcripts

# # COUNTS
# cmc_matched_cts <- lapply(groups, function(g) {
#     rowMeans(cmc_cts[overlap_transcripts, filter(cmc_metadata, Group == g) %>% pull(Sample), drop = FALSE])
# }) %>% 
#     bind_cols() %>%
#     setNames(groups)
# rownames(cmc_matched_cts) <- overlap_transcripts
# bspan_matched_cts <- lapply(groups, function(g) {
#     rowMeans(brainspan_cts[overlap_transcripts, filter(brainspan_metadata, Group == g) %>% pull(Sample), drop = FALSE])
# }) %>% 
#     bind_cols() %>%
#     setNames(groups)
# rownames(bspan_matched_cts) <- overlap_transcripts

empirical_correlations <- sapply(setNames(nm = overlap_transcripts), function(x) {
    cor(as.numeric(cmc_matched[x, groups]), as.numeric(bspan_matched[x, groups]))
})
empirical_correlations <- empirical_correlations[!is.na(empirical_correlations)]
# 
# saveRDS(empirical_correlations, file = "data/psychencode/BrainSpan-CMC-HBCC_empirical_correlations.rds")
# saveRDS(bspan_matched, file = "data/psychencode/CMC-HBCC_matched_transcript_tpm.rds")
# saveRDS(cmc_matched, file = "data/psychencode/BrainSpan-CMC-HBCC_matched_transcript_tpm.rds")

# null_permutations <- mclapply(1:10, function(i) {
#     tibble(
#         IterGroup = Sys.getenv("SLURM_ARRAY_TASK_ID"),
#         Iter = i,
#         corr = sapply(setNames(nm = overlap_transcripts), function(x) {
#             cor(as.numeric(cmc_matched[sample(overlap_transcripts, size = 1), groups]), as.numeric(bspan_matched[x, groups]))
#         })
#     )
# }) %>%
#     bind_rows() %>%
#     mutate(
#         x = cut(corr, seq(-1, 1, by = 0.1))
#     ) %>%
#     group_by(IterGroup, Iter) %>%
#     count(x)
# saveRDS(null_permutations, file = paste0("data/psychencode/BrainSpan-CMC-HBCC_null_correlations_", Sys.getenv("SLURM_ARRAY_TASK_ID"), ".rds"))