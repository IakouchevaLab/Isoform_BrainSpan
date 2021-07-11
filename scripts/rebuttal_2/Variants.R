library(tidyverse)

variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
)

splice_variants <- variants %>%
    filter(ifelse(is.na(str_match(Consequence, "splice_region")), TRUE, FALSE)) %>%
    filter(ifelse(is.na(str_match(Consequence, "splice")), FALSE, TRUE))

splice_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`) %>%
    count(`Affected status` == 1)
splice_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, `Ensembl Transcript ID`) %>%
    filter(`Affected status` == 1) %>%
    distinct(`Ensembl Transcript ID`)
splice_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, `Ensembl Transcript ID`) %>%
    filter(`Affected status` == 2) %>%
    distinct(`Ensembl Transcript ID`)


spl_variants <- variants %>%
    separate_rows(Consequence, sep = ',') %>%
    filter(
        Consequence %in% c(
            "splice_donor_variant", "splice_acceptor_variant"
        )
    ) %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, Consequence, `Ensembl Transcript ID`)
spl_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`) %>%
    count(`Affected status`)


lof_variants <- variants %>%
    separate_rows(Consequence, sep = ',') %>%
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained",
            "splice_donor_variant", "splice_acceptor_variant"
        )
    ) %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, Consequence, `Ensembl Transcript ID`)
lof_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`) %>%
    count(`Affected status`)
lof_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, `Ensembl Transcript ID`) %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Transcript ID`) %>%
    unique() %>%
    length()
lof_variants %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, `Ensembl Transcript ID`) %>%
    filter(`Affected status` == 1) %>%
    pull(`Ensembl Transcript ID`) %>%
    unique() %>%
    length()

variants %>%
    filter(`Variant start` == 51585462) %>%
    View()

lof_spl_variants <- variants %>%
    separate_rows(Consequence, sep = ',') %>%
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained",
            "splice_donor_variant", "splice_acceptor_variant"
        )
    ) %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, Consequence) %>%
    mutate(
        spl = ifelse(is.na(str_match(Consequence, "splice")), FALSE, TRUE)
    ) %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, spl) %>%
    count(Chromosome, `Variant start`, Reference, Alternate, `Affected status`)
    count(Chromosome, `Variant start`, Reference, Alternate, Consequence) %>%
    filter(n > 1) %>%
    left_join(variants, by = c("Chromosome", "Variant start", "Reference", "Alternate")) %>%
    distinct(Chromosome, `Variant start`, Reference, Alternate, `Affected status`, Consequence)


library(bioMart)
ensembl_transcript_id <- rownames(readRDS("data/iso_tpm_filter.rds"))
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
missense_variants %>%
    distinct(chromosome_name, position, reference, alternate, affected_status) %>%
    count(affected_status)
missense_variants %>%
    filter(affected_status == 2) %>%
    pull(affected_transcripts) %>%
    unlist() %>%
    unique() %>%
    length()
missense_variants %>%
    filter(affected_status == 1) %>%
    pull(affected_transcripts) %>%
    unlist() %>%
    unique() %>%
    length()
