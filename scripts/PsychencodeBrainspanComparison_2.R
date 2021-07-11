library(tidyverse)
brainspan_cortex_regions <- c(
    "OFC", "DFC", "VCF", "MFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "OCX"
)
brainspan_tpm <- readRDS("data/iso_tpm_filter.rds")
brainspan_metadata <- read_tsv("data/metadata.tsv") %>%
    mutate(standardized_region = ifelse(Regioncode %in% brainspan_cortex_regions, "cortex", NA)) %>%
    mutate(year = floor(Days / 365)) %>%
    mutate(age_interval = as.character(cut(year, seq(-11, 100, by = 10)))) %>%
    mutate(standardized_age = sapply(age_interval, function(i) {
        paste0(
            as.numeric(gsub("^\\(([-0-9]+),.+", "\\1", i)) + 1,
            "-",
            as.numeric(gsub(".+,([0-9]+)\\]$", "\\1", i))
        )
    })) %>%
    mutate(standardized_sex = ifelse(Sex == "M", "male", "female")) %>%
    mutate(standardized_groups = paste0(standardized_region, standardized_age, standardized_sex)) %>%
    filter(!is.na(standardized_region)) %>%
    dplyr::select(Sample, starts_with("standardized"))

load("data/psychencode/psychencode_isoform_clean_input/Isoform.TPM.noOutliersamples.AllIsoforms.RData")
psychencode_metadata <- read.csv("data/psychencode/psychencode_isoform_clean_input/Capstone_datMeta.csv", row.names = 1) %>%
    mutate(SampleName = rownames(.)) %>%
    filter(diagnosis != "Control") %>%
    mutate(standardized_region = "cortex") %>%
    mutate(age_interval = as.character(cut(ageDeath, seq(-11, 100, by = 10)))) %>%
    mutate(standardized_age = sapply(age_interval, function(i) {
        paste0(
            as.numeric(gsub("^\\(([-0-9]+),.+", "\\1", i)) + 1,
            "-",
            as.numeric(gsub(".+,([0-9]+)\\]$", "\\1", i))
        )
    })) %>%
    mutate(standardized_sex = ifelse(sex == "M", "male", "female")) %>%
    mutate(standardized_groups = paste0(standardized_region, standardized_age, standardized_sex)) %>%
    filter(!is.na(standardized_region)) %>%
    dplyr::select(SampleName, starts_with("standardized"), study)
psychencode_tpm <- tpm[, which(colnames(tpm) %in% pull(psychencode_metadata, SampleName))]
rm(tpm)
psychencode_metadata <- psychencode_metadata %>%
    filter(SampleName %in% colnames(psychencode_tpm))

gtex_con <- file('data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct', 'r')
gtex_header <- toupper(str_split(readLines(gtex_con, n = 3)[3], '\t')[[1]])
close(gtex_con)
gtex_metadata <- read_tsv(
    "data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    col_types = c('.default' = 'c')
) %>%
    mutate(SAMPID = toupper(SAMPID)) %>%
    dplyr::select(SAMPID, SMTS, SMTSD) %>%
    filter(SMTS == 'Brain') %>%
    mutate(SUBJID = sapply(str_split(SAMPID, pattern = "-"), function(x) paste(x[1:2], collapse = '-'))) %>%
    left_join(
        read_tsv("data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
    ) %>%
    mutate(standardized_sex = ifelse(SEX == 1, 'male', 'female')) %>%
    mutate(standardized_region = ifelse(str_detect(SMTSD, "Cortex"), "cortex", NA)) %>%
    mutate(standardized_age = AGE) %>%
    mutate(standardized_groups = paste0(standardized_region, standardized_age, standardized_sex)) %>%
    filter(toupper(SAMPID) %in% toupper(gtex_header)) %>%
    filter(!is.na(standardized_region)) %>%
    dplyr::select(SAMPID, starts_with("standardized"))
gtex_brain_samples_indexes <- paste(which(toupper(gtex_header) %in% toupper(gtex_metadata[["SAMPID"]])), collapse = ',')
#system(command = paste0("tail -n +3 data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | cut -f1,2,", gtex_brain_samples_indexes, " > data/gtex_analysis/brain_samples.tsv"))
gtex_tpm <- read_tsv("data/gtex_analysis/brain_samples.tsv") %>%
    mutate(ensembl_transcript_id = str_split(transcript_id, "\\.", simplify = TRUE)[, 1]) %>%
    mutate(ensembl_transcript_id_version = as.numeric(str_split(transcript_id, "\\.", simplify = TRUE)[, 2])) %>%
    mutate(ensembl_transcript_id_version = ifelse(is.na(ensembl_transcript_id_version), 1, ensembl_transcript_id_version)) %>%
    group_by(ensembl_transcript_id) %>%
    slice(which.max(ensembl_transcript_id_version))
gtex_txid <- gtex_tpm[["ensembl_transcript_id"]]
gtex_tpm <- gtex_tpm[, which(! colnames(gtex_tpm) %in% c("transcript_id", "gene_id", "ensembl_transcript_id", "ensembl_transcript_id_version"))]
colnames(gtex_tpm) <- toupper(colnames(gtex_tpm))
rownames(gtex_tpm) <- gtex_txid

brainspan_gtex_groups <- intersect(brainspan_metadata$standardized_groups, gtex_metadata$standardized_groups)
psychencode_studies <- as.character(sort(unique(psychencode_metadata$study)))

intersect_with_studies <- lapply(setNames(nm = psychencode_studies), function(ps) {
    tibble(
        PsychencodeStudy = ps,
        GroupsCommonBrainspanGTEx = sort(intersect(brainspan_gtex_groups, pull(filter(psychencode_metadata, study == ps), standardized_groups)))
    )
}) %>%
    bind_rows()

writexl::write_xlsx(
    x = list(
        "BrainSpanStandardizedGroups" = brainspan_metadata, 
        "GTExStandardizedGroups" = gtex_metadata, 
        "PsychencodeStandardizedGroups" = psychencode_metadata,
        "IntersectThreeSets" = intersect_with_studies
    ),
    path = "data/TranscriptCorrelationStandardizedSamples.xlsx"
)

overlap_transcripts <- intersect(intersect(rownames(brainspan_tpm), rownames(gtex_tpm)), rownames(psychencode_tpm))
# overlap_transcripts <- sample(
#     intersect(intersect(rownames(brainspan_tpm), rownames(gtex_tpm)), rownames(psychencode_tpm)),
#     size = 1000
# )

standardize_tpm <- function(tpm, sample_maps, sample_names) {
    df <- lapply(sample_names, function(s) {
        rowMeans(tpm[, pull(filter(sample_maps, standardized_group == s), sample), drop = FALSE])
    }) %>%
        bind_cols() %>%
        setNames(sample_names)
    rownames(df) <- rownames(tpm)
    return(df)
}

brainspan_tpm_standardized <- standardize_tpm(
    brainspan_tpm,
    tibble(sample = brainspan_metadata$Sample, standardized_group = brainspan_metadata$standardized_groups),
    sort(unique(brainspan_metadata$standardized_groups))
)
gtex_tpm_standardized <- standardize_tpm(
    gtex_tpm,
    tibble(sample = gtex_metadata$SAMPID, standardized_group = gtex_metadata$standardized_groups),
    sort(unique(gtex_metadata$standardized_groups))
)

psychencode_brainspan_correlation_data <- list()
psychencode_gtex_correlation_data <- list()
brainspan_gtex_correlation_data <- list()
for (ps in psychencode_studies) {
    sample_maps <- psychencode_metadata %>%
        filter(study == ps) %>%
        dplyr::select(SampleName, standardized_groups) %>%
        rename(sample = SampleName, standardized_group = standardized_groups)
    psychencode_tpm_standardized <- standardize_tpm(psychencode_tpm, sample_maps, sort(unique(sample_maps$standardized_group)))
    group_intersect <- intersect(
        intersect(
            colnames(brainspan_tpm_standardized), colnames(gtex_tpm_standardized)
        ),
        colnames(psychencode_tpm_standardized)
    )
    psychencode_brainspan_correlation_data[[ps]] <- sapply(setNames(nm = overlap_transcripts), function(x) {
        tryCatch({
            cor(
                as.numeric(brainspan_tpm_standardized[x, group_intersect]), 
                as.numeric(psychencode_tpm_standardized[x, group_intersect])
            )
        },
        error = function(e) { e },
        warning = function(w) { 
            0
        })
    })
    psychencode_gtex_correlation_data[[ps]] <- sapply(setNames(nm = overlap_transcripts), function(x) {
        tryCatch({
            cor(
                as.numeric(gtex_tpm_standardized[x, group_intersect]),
                as.numeric(psychencode_tpm_standardized[x, group_intersect])
            )
        },
        warning = function(w) {
            0
        })
        
    })
    brainspan_gtex_correlation_data[[ps]] <- sapply(setNames(nm = overlap_transcripts), function(x) {
        tryCatch({
            cor(
                as.numeric(gtex_tpm_standardized[x, group_intersect]),
                as.numeric(brainspan_tpm_standardized[x, group_intersect])
            )
        },
        warning = function(w) {
            0
        })
    })
}

plot_data <- bind_rows(
    data.frame(table(cut(psychencode_brainspan_correlation_data$LIBD_szControl, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsBrainSpan") %>%
        mutate(Group = "LIBD_szControl") %>%
        mutate(Column = "PsychencodeVs"),
    data.frame(table(cut(psychencode_gtex_correlation_data$LIBD_szControl, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsGTEx") %>%
        mutate(Group = "LIBD_szControl") %>%
        mutate(Column = "PsychencodeVs"),
    data.frame(table(cut(psychencode_brainspan_correlation_data$CMC, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsBrainSpan") %>%
        mutate(Group = "CMC") %>%
        mutate(Column = "PsychencodeVs"),
    data.frame(table(cut(psychencode_gtex_correlation_data$CMC, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsGTEx") %>%
        mutate(Group = "CMC") %>%
        mutate(Column = "PsychencodeVs"),
    data.frame(table(cut(psychencode_brainspan_correlation_data$CMC_HBCC, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsBrainSpan") %>%
        mutate(Group = "CMC_HBCC") %>%
        mutate(Column = "PsychencodeVs"),
    data.frame(table(cut(psychencode_gtex_correlation_data$CMC_HBCC, seq(-1, 1, by = 0.1)))) %>%
        mutate(Comparison = "vsGTEx") %>%
        mutate(Group = "CMC_HBCC") %>%
        mutate(Column = "PsychencodeVs")
)
ggplot(
    data = plot_data,
    aes(
        x = Var1, y = Freq, colour = Comparison
    )
) +
    facet_grid(Group ~ .) +
    geom_point(size = 3) +
    labs(
        x = "Correlation Bin", y = "Frequency"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

brainspan_gtex <- sapply(setNames(nm = overlap_transcripts), function(x) {
    tryCatch({
        cor(
            as.numeric(gtex_tpm_standardized[x, c("cortex20-29female", "cortex20-29male", "cortex30-39female", "cortex30-39male", "cortex40-49female")]), 
            as.numeric(brainspan_tpm_standardized[x, c("cortex20-29female", "cortex20-29male", "cortex30-39female", "cortex30-39male", "cortex40-49female")])
        )
    },
    warning = function(w) { 0 })
})

ggplot() +
    geom_point(
        data = data.frame(table(cut(brainspan_gtex, seq(-1, 1, by = 0.1)))),
        aes(x = Var1, y = Freq)
    )

mean_brainspan_tpm <- rowMeans(brainspan_tpm_standardized)
mean_gtex_tpm <- rowMeans(gtex_tpm_standardized)
correlation_plts <- list()
ggplot(
    data = tibble(
        ensembl_transcript_id = names(brainspan_gtex_correlation_data[["CMC"]]),
        correlation = brainspan_gtex_correlation_data[["CMC"]]
    ) %>%
        filter(!is.na(correlation)) %>%
        mutate(interval = cut(correlation, breaks = seq(-1, 1, by = 0.1))) %>%
        mutate(mean_brainspan_tpm = mean_brainspan_tpm[ensembl_transcript_id]) %>%
        mutate(mean_gtex_tpm = mean_gtex_tpm[ensembl_transcript_id]) %>%
        pivot_longer(cols = c(mean_brainspan_tpm, mean_gtex_tpm), names_to = "dataset", values_to = "mean_tpm"),
    aes(
        x = interval, y = log2(mean_tpm), fill = dataset
    )
) +
    # scale_y_continuous(limits = c(0, 10)) +
    geom_boxplot(position = "dodge") +
    labs(x = "correlation_bin") +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 1)
    )




# sample_maps_CMC <- psychencode_metadata %>%
#     filter(study == "CMC") %>%
#     dplyr::select(SampleName, standardized_groups) %>%
#     rename(sample = SampleName, standardized_group = standardized_groups)
# CMC_tpm_standardized <- standardize_tpm(psychencode_tpm, sample_maps_CMC, sort(unique(sample_maps_CMC$standardized_group)))
# sample_maps_LIBD <- psychencode_metadata %>%
#     filter(study == "LIBD_szControl") %>%
#     dplyr::select(SampleName, standardized_groups) %>%
#     rename(sample = SampleName, standardized_group = standardized_groups)
# LIBD_tpm_standardized <- standardize_tpm(psychencode_tpm, sample_maps_LIBD, sort(unique(sample_maps_LIBD$standardized_group)))
# CMC_LIBD_correlations <- sapply(setNames(nm = overlap_transcripts), function(x) {
#     cor(
#         as.numeric(CMC_tpm_standardized[x, c("cortex20-29female", "cortex20-29male", "cortex30-39female", "cortex30-39male", "cortex40-49female")]), 
#         as.numeric(LIBD_tpm_standardized[x, c("cortex20-29female", "cortex20-29male", "cortex30-39female", "cortex30-39male", "cortex40-49female")])
#     )
# })
# ggplot() +
#     geom_point(
#         data = data.frame(table(cut(CMC_LIBD_correlations, seq(-1, 1, by = 0.1)))),
#         aes(x = Var1, y = Freq)
#     )
# median(CMC_LIBD_correlations[!is.na(CMC_LIBD_correlations)])

