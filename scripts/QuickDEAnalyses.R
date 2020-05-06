library(tidyverse)

tt_combined <- read_tsv(
    "data/limmaResultsAll_LoFTargets.tsv",
    col_types = cols(.default = "c")
)

metadata <- read_tsv("data/metadata.tsv")

nested_significant_hits <- tt_combined %>%
    filter(as.logical(Significant)) %>%
    mutate(
        Feature = ifelse(
            Data == "Gene", ensembl_gene_id, ensembl_transcript_id
        )
    ) %>%
    dplyr::select(Data, Contrast, Feature) %>%
    group_by(Data, Contrast) %>%
    nest(Feature, .key = "Feature") %>%
    ungroup() %>%
    mutate(
        TotalFeatures = ifelse(Data == "Gene", 26307, 100754)
    ) %>%
    mutate(
        ProportionFeatures = mapply(
            function(feat, tot) nrow(feat) / tot,
            Feature, TotalFeatures
        )
    )
distinct_gene_prenatal_de <- nested_significant_hits %>%
    filter(Data == "Gene") %>%
    filter(Contrast %in% c(
        "P02P03", "P03P04", "P04P05", "P05P06", "P06P07"
    )) %>%
    pull(Feature) %>%
    lapply(function(f) pull(f, Feature)) %>%
    unlist() %>%
    unique() %>%
    length()
distinct_gene_prenatal_de / 26307
distinct_gene_postnatal_de <- nested_significant_hits %>%
    filter(Data == "Gene") %>%
    filter(Contrast %in% c(
        "P08P09", "P09P10", "P10P11", "P11P12", "P12P13"
    )) %>%
    pull(Feature) %>%
    lapply(function(f) pull(f, Feature)) %>%
    unlist() %>%
    unique() %>%
    length()
distinct_gene_postnatal_de / 26307

distinct_isoform_prenatal_de <- nested_significant_hits %>%
    filter(Data == "Isoform") %>%
    filter(Contrast %in% c(
        "P02P03", "P03P04", "P04P05", "P05P06", "P06P07"
    )) %>%
    pull(Feature) %>%
    lapply(function(f) pull(f, Feature)) %>%
    unlist() %>%
    unique() %>%
    length()
distinct_isoform_prenatal_de / 26307
distinct_isoform_postnatal_de <- nested_significant_hits %>%
    filter(Data == "Isoform") %>%
    filter(Contrast %in% c(
        "P08P09", "P09P10", "P10P11", "P11P12", "P12P13"
    )) %>%
    pull(Feature) %>%
    lapply(function(f) pull(f, Feature)) %>%
    unlist() %>%
    unique() %>%
    length()
distinct_isoform_postnatal_de / 26307



overlaps <- t(combn(
    sort(unique(nested_significant_hits$Contrast)), 2
)) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rename(CTR1 = V1, CTR2 = V2) %>%
    filter(CTR1 != CTR2) %>%
    mutate(Data = "Gene") %>%
    bind_rows(
        mutate(., Data = "Isoform")
    )
