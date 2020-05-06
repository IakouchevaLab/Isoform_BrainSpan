library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores())

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
    filter(SAMPID %in% gtex_header)
gtex_sample_phenotypes <- read_tsv(
    "data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
    col_types = c('.default' = 'c')
) %>%
    separate(AGE, into = c("AGE_MIN", "AGE_MAX"), sep = '-') %>%
    filter(AGE_MIN <= 40)
gtex_sample_metadata <- inner_join(gtex_sample_attributes, gtex_sample_phenotypes, by = "SUBJID")
brainspan_sample_metadata <- read_tsv("data/metadata.tsv")
adult_samples <- brainspan_sample_metadata %>%
    filter(Days >= 7931) %>%
    mutate(AgeYears = as.numeric(str_replace(Age, 'Y', '')))

SAMPLE_MAP <- tibble(
    BRAINSPAN = c("A1C",
                  "AMY",
                  "CBC",
                  "DFC",
                  "HIP",
                  "IPC",
                  "ITC",
                  "M1C",
                  "MD" ,
                  "MFC",
                  "OFC",
                  "S1C",
                  "STC",
                  "STR",
                  "STR",
                  "V1C",
                  "VFC"),
    GTEX = c("Brain - Cortex",
             "Brain - Amygdala",
             "Brain - Cerebellar Hemisphere",
             "Brain - Frontal Cortex (BA9)",
             "Brain - Hippocampus",
             "Brain - Cortex",
             "Brain - Cortex",
             "Brain - Cortex",
             "Brain - Hypothalamus",
             "Brain - Cortex",
             "Brain - Frontal Cortex (BA9)",
             "Brain - Cortex",
             "Brain - Cortex",
             "Brain - Caudate (basal ganglia)",
             "Brain - Putamen (basal ganglia)",
             "Brain - Cortex",
             "Brain - Frontal Cortex (BA9)")
)

# Map samples to groups based on region, age, and sex (Using broadest classifiers)

sample_groups <- expand.grid(
    AgeMin = c(20, 30, 40),
    AgeMax = c(29, 39, 49),
    Region = unique(SAMPLE_MAP$BRAINSPAN),
    Sex = c("M", "F")
) %>%
    as_tibble() %>%
    mutate(SAMPLENAME = paste(Region, AgeMin, AgeMax, Sex, sep = "_")) %>%
    mutate(
        BrainSpanSamples = apply(
            .,
            1,
            function(x) {
                adult_samples %>%
                    filter(AgeYears >= x[["AgeMin"]]) %>%
                    filter(AgeYears <= x[["AgeMax"]]) %>%
                    filter(Regioncode == x[["Region"]]) %>%
                    filter(Sex == x[["Sex"]]) %>%
                    pull(Sample)
            }
        )
    ) %>%
    mutate(GTExSamples = apply(
        ., 1, function(x) {
            gtex_sample_metadata %>%
                filter(AGE_MIN == x[['AgeMin']]) %>%
                filter(AGE_MAX == x[['AgeMax']]) %>%
                filter(SMTSD %in% SAMPLE_MAP[SAMPLE_MAP$BRAINSPAN == x[['Region']], 'GTEX']) %>%
                filter(as.numeric(SEX) == as.numeric((x[['Sex']] == 'M'))) %>%
                pull(SAMPID)
        }
    )) %>%
    filter(
        mapply(
            function(b, g) {
                length(b) > 0 && length(g) > 0
            },
            BrainSpanSamples, GTExSamples
        )
    )

# Subset data and convert to grouped identifiers
all_brainspan_samples <- unlist(sample_groups$BrainSpanSamples)
all_gtex_samples <- unlist(sample_groups$GTExSamples)

subset_brainspan <- readRDS("data/iso_tpm_filter.rds")[, all_brainspan_samples]
header_idx <- paste(which(gtex_header %in% unlist(sample_groups$GTExSamples)), collapse = ',')
#system(command = paste0("tail -n +3 data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | cut -f1,2,", header_idx, " > data/gtex_analysis/adult_brain_samples.tsv"))
subset_gtex <- read_tsv("data/gtex_analysis/adult_brain_samples.tsv")

transcripts <- intersect(
    rownames(subset_brainspan),
    str_replace(subset_gtex$transcript_id, "\\.\\d+", "")
)

brainspan_group <- lapply(sample_groups$SAMPLENAME, function(x) {
    samples <- sample_groups %>%
        filter(SAMPLENAME == x) %>%
        pull(BrainSpanSamples) %>%
        unlist()
    row_means <- rowMeans(subset_brainspan[, samples, drop = FALSE])
}) %>%
    bind_cols() %>%
    setNames(sample_groups$SAMPLENAME) %>%
    mutate(transcript_id = rownames(subset_brainspan)) %>%
    filter(transcript_id %in% transcripts)

gtex_group <- lapply(sample_groups$SAMPLENAME, function(x) {
    samples <- sample_groups %>%
        filter(SAMPLENAME == x) %>%
        pull(GTExSamples) %>%
        unlist()
    row_means <- rowMeans(subset_gtex[, samples, drop = FALSE])
}) %>%
    bind_cols() %>%
    setNames(sample_groups$SAMPLENAME) %>%
    mutate(transcript_id = str_replace(subset_gtex$transcript_id, "\\.\\d+", "")) %>%
    filter(transcript_id %in% transcripts)

saveRDS(brainspan_group, "data/gtex_analysis/BrainSpanTPM_matchGTEx.rds")
saveRDS(gtex_group, "data/gtex_analysis/GTExTPM_matchBrainSpan.rds")

# Create empirical correlation data
correlation_data <- mclapply(setNames(nm = transcripts), function(tx) {
    samples <- sample_groups$SAMPLENAME
    brainspan_data <- as.numeric(brainspan_group[brainspan_group$transcript_id == tx, samples])
    gtex_data <- as.numeric(gtex_group[gtex_group$transcript_id == tx, samples])
    correlation <- cor.test(brainspan_data, gtex_data, method = "pearson") %>%
        broom::tidy() %>%
        mutate(transcript_id = tx)
    return(correlation)
}) %>%
    bind_rows()

write_tsv(correlation_data, "data/gtex_analysis/empirical_correlations.tsv")
