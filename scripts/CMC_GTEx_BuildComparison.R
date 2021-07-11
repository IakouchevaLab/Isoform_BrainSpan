# Subset and consolidate overlapping samples between our BrainSpan data and the Psychencode cmc data
# Normal controls

library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

bspan_tx <- rownames(readRDS("data/iso_tpm_filter.rds"))

# Set up GTEx
gtex_con <- file('data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct', 'r')
gtex_header <- str_split(readLines(gtex_con, n = 3)[3], '\t')[[1]]
close(gtex_con)
gtex_sample_attributes <- read_tsv(
    "data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    col_types = c('.default' = 'c')
) %>%
    dplyr::select(SAMPID, SMTS, SMTSD) %>%
    filter(SMTS == 'Brain') %>%
    mutate(SUBJID = sapply(str_split(SAMPID, pattern = "-"), function(x) paste(x[1:2], collapse = '-'))) %>%
    left_join(
        read_tsv("data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
    ) %>%
    mutate(SEX = ifelse(SEX == 1, 'male', 'female')) %>%
    mutate(Group = mapply(function(r, a, s) { if(str_detect(r, "Cortex")) { return(paste0("cortex", a, s)) } else { return(NA) } }, SMTSD, AGE, SEX)) %>%
    filter(toupper(SAMPID) %in% toupper(gtex_header))
gtex_brain_samples_indexes <- paste(which(gtex_header %in% gtex_sample_attributes[["SAMPID"]]), collapse = ',')
#system(command = paste0("tail -n +3 data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | cut -f1,2,", gtex_brain_samples_indexes, " > data/gtex_analysis/brain_samples.tsv"))
gtex_tpm_prefilt <- read_tsv("data/gtex_analysis/brain_samples.tsv") %>%
    mutate(ensembl_transcript_id = str_split(transcript_id, "\\.", simplify = TRUE)[, 1]) %>%
    mutate(ensembl_transcript_id_version = as.numeric(str_split(transcript_id, "\\.", simplify = TRUE)[, 2])) %>%
    mutate(ensembl_transcript_id_version = ifelse(is.na(ensembl_transcript_id_version), 1, ensembl_transcript_id_version)) %>%
    group_by(ensembl_transcript_id)
gtex_tpm <- gtex_tpm_prefilt %>%
    slice(which.max(ensembl_transcript_id_version))
gtex_txid <- gtex_tpm[["ensembl_transcript_id"]]
gtex_tpm <- gtex_tpm[, which(! colnames(gtex_tpm) %in% c("transcript_id", "gene_id", "ensembl_transcript_id", "ensembl_transcript_id_version"))]
colnames(gtex_tpm) <- toupper(colnames(gtex_tpm))
rownames(gtex_tpm) <- gtex_txid

# Setup CMC
cmc_metadata <- read.csv('data/psychencode/psychencode_isoform_clean_input/Capstone_datMeta.csv', row.names=1) %>%
    mutate(Sample = rownames(.)) %>%
    filter(diagnosis == 'Control') %>%
    filter(study == 'CMC_HBCC') %>% 
    mutate(RegionMatch = as.character(tissue)) %>%
    mutate(REGION = ifelse(str_detect(RegionMatch, 'frontal cortex'), 'frontal cortex', 'NA')) %>%
    mutate(AgeMatch = as.character(cut(ageDeath, seq(-11, 100, by = 10)))) %>%
    mutate(AGE_MIN = as.numeric(str_match(string = AgeMatch, pattern = "^\\((-*\\d+)")[, 2]) + 1) %>%
    mutate(AGE_MAX = as.numeric(str_match(string = AgeMatch, pattern = "(-*\\d+)\\]$")[, 2])) %>%
    mutate(AGE = paste(AGE_MIN, AGE_MAX, sep = '-')) %>%
    mutate(SEX = ifelse(sex == 'M', 'male', 'female')) %>%
    mutate(Group = paste0(REGION, AGE, SEX))
load('data/psychencode/psychencode_isoform_clean_input/Isoform.TPM.noOutliersamples.AllIsoforms.RData')
colnames(tpm) = str_replace(colnames(tpm), pattern = "\\.", replacement = "-")
cmc_tpm <- tpm[, str_detect(colnames(tpm), "CMC_HBCC.+")]

# Subsetting for matched data
# TPM
groups <- sort(intersect(cmc_metadata$Group, gtex_sample_attributes$Group))
groups <- groups[1:6]
overlap_transcripts <- intersect(intersect(rownames(cmc_tpm), rownames(gtex_tpm)), bspan_tx)
cmc_matched <- lapply(groups, function(g) {
    rowMeans(cmc_tpm[overlap_transcripts, filter(cmc_metadata, Group == g) %>% pull(Sample), drop = FALSE])
}) %>% 
    bind_cols() %>%
    setNames(groups)
rownames(cmc_matched) <- overlap_transcripts
gtex_matched <- lapply(groups, function(g) {
    rowMeans(gtex_tpm[overlap_transcripts, filter(gtex_sample_attributes, Group == g) %>% filter(toupper(SAMPID) %in% toupper(colnames(gtex_tpm))) %>% pull(SAMPID) %>% toupper(), drop = FALSE])
}) %>%
    bind_cols() %>%
    setNames(groups)
rownames(gtex_matched) <- overlap_transcripts

empirical_correlations <- sapply(setNames(nm = overlap_transcripts), function(x) {
    cor(as.numeric(cmc_matched[x, groups]), as.numeric(gtex_matched[x, groups]))
})
empirical_correlations <- empirical_correlations[!is.na(empirical_correlations)]

ggplot(
    data = tibble(
        interval = cut(empirical_correlations, seq(-1, 1, by = 0.1))
    ) %>%
        count(interval),
    aes(x = interval, y = n)
) +
    geom_point() +
    theme(
        text = element_text(size = 24)
    )

# saveRDS(empirical_correlations, file = "data/gtex_analysis/GTEx-CMC-HBCC_empirical_correlations.rds")
# saveRDS(gtex_matched, file = "data/gtex_analysis/GTEx-CMC-HBCC_matched_transcript_tpm.rds")
# saveRDS(cmc_matched, file = "data/gtex_analysis/CMC-HBCC-GTEx_matched_transcript_tpm.rds")

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