library(tidyverse)
library(biomaRt)

iso_tpm <- readRDS("data/iso_tpm_filter.rds")
metadata <- read_tsv("data/metadata.tsv")
ensembl_transcript_id <- rownames(iso_tpm)
mart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)
transcript_coordinates <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_exon_id",
                   "chromosome_name", "exon_chrom_start", "exon_chrom_end"),
    filters = c("ensembl_transcript_id"),
    values = ensembl_transcript_id,
    mart = mart
) %>%
    as_tibble()

missense_variants <- readxl::read_xlsx(
    "data/Satterstrom_DNMs_filtered.xlsx", sheet = 1
) %>%
    dplyr::select(
        Chromose_number, Position, 
        Reference_allele, Alternate_allele, 
        Child_Sex, Affected_Status, VEP_functional_class_canonical_simplified
    ) %>%
    rename(
        chromosome_name = Chromose_number,
        position = Position,
        reference = Reference_allele,
        alternate = Alternate_allele,
        sex = Child_Sex,
        affected_status = Affected_Status,
        consequence = VEP_functional_class_canonical_simplified
    ) %>%
    mutate(
        affected_transcripts = apply(., 1, function(v) {
            transcript_coordinates %>%
                filter(chromosome_name == v[["chromosome_name"]]) %>%
                filter(exon_chrom_start <= v[["position"]]) %>%
                filter(exon_chrom_end >= v[["position"]]) %>%
                pull(ensembl_transcript_id) %>%
                unique()
        })
    )

missense_asd_transcripts <- unique(unlist(pull(filter(missense_variants, affected_status == 2), affected_transcripts)))
missense_ctr_transcripts <- unique(unlist(pull(filter(missense_variants, affected_status == 1), affected_transcripts)))
nonimpact_transcripts <- rownames(iso_tpm)[which(!rownames(iso_tpm) %in% c(missense_asd_transcripts, missense_ctr_transcripts))]

source("scripts/utility/plot_expression.R")

missense_asd_tpm_z <- expr_to_z(iso_tpm[missense_asd_transcripts, ]) %>%
    as.data.frame() %>%
    mutate(ensembl_transcript_id = rownames(.)) %>%
    gather(
        key = "Sample", value = "TPM Z-Score", 
        starts_with("BrainSpan")
    ) %>%
    left_join(metadata, by = "Sample")
missense_asd_plt <- ggplot(
    data = missense_asd_tpm_z,
    mapping = aes(
        x = Days, y = `TPM Z-Score`
    )
) +
    geom_line(
        stat = "smooth", method = "loess", size = 2
    )

iso_tpm_byPeriod <- lapply(
    sort(unique(metadata$Period)),
    function(p) {
        samples <- metadata %>%
            filter(Period == p) %>%
            pull(Sample)
        rowMeans(iso_tpm[, samples, drop = FALSE])
    }
) %>%
    bind_cols() %>%
    setNames(nm = sprintf("P%02d", sort(unique(metadata$Period)))) %>%
    mutate(ensembl_transcript_id = ensembl_transcript_id)


