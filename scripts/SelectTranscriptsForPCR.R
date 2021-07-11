library(tidyverse)

iso_tpm <- readRDS("data/iso_tpm_filter.rds")
metadata <- read_tsv("data/metadata.tsv")
annotations <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]

female_22pcw_frontal_lobe <- metadata %>%
    filter(Regioncode %in% c("OFC", "DFC", "VFC", "MFC", "M1C")) %>%
    filter(Period == 6) %>%
    filter(Sex == "F") %>%
    pull(Sample)
female_27y_frontal_cortex <- metadata %>%
    filter(Regioncode %in% c("OFC", "DFC", "VFC", "MFC", "M1C")) %>%
    filter(Period == 13) %>%
    filter(Sex == "F") %>%
    pull(Sample)
female_22pcw_frontal_lobe_tpm <- rowMeans(iso_tpm[, female_22pcw_frontal_lobe])
female_27y_frontal_cortex_tpm <- rowMeans(iso_tpm[, female_27y_frontal_cortex])

isoforms <- rownames(iso_tpm)
delta_iso <- log2(female_22pcw_frontal_lobe_tpm[isoforms]/(female_27y_frontal_cortex_tpm+0.0001))

isoforms_for_qpcr <- tibble(
    ensembl_transcript_id = isoforms
) %>%
    left_join(
        dplyr::select(annotations, ensembl_transcript_id, ensembl_gene_id), 
        by = "ensembl_transcript_id"
    ) %>%
    mutate(
        female_22pcw_frontal_lobe_tpm = female_22pcw_frontal_lobe_tpm,
        female_27y_frontal_cortex_tpm = female_27y_frontal_cortex_tpm,
        log2fc = delta_iso
    ) %>%
    filter(female_22pcw_frontal_lobe_tpm != 0) %>%
    filter(female_27y_frontal_cortex_tpm != 0) %>%
    mutate(abs_log2fc = abs(log2fc)) %>%
    arrange(desc(abs_log2fc))
write_csv(isoforms_for_qpcr, "data/isoforms_for_qpcr.csv")

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "GRCh37.ensembl.org")
exons <- getBM(
    attributes = c(
        "ensembl_gene_id", "ensembl_transcript_id", "transcript_biotype",
        "start_position", "end_position",
        "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "strand"
    ),
    filters = c("ensembl_transcript_id"),
    values = read.csv("data/source/annotation.transcript.ensg75.txt")$ensembl_transcript_id,
    mart = mart
)
 
calculate_overlap <- function(start1, end1, start2, end2) {
    return(max(0, min(c(end1, end2)) - max(c(start1, start2))))
}
isoforms_for_qpcr <- isoforms_for_qpcr %>%
    left_join(
        dplyr::select(exons, ensembl_transcript_id, transcript_biotype, ensembl_gene_id, start_position, end_position)
    ) %>%
    distinct() %>%
    filter(transcript_biotype == "protein_coding")
tx_unique_region <- c()
ex_unique_region <- c()
ln_unique_region <- c()
for (row in 1:nrow(isoforms_for_qpcr)) {
    tx <- pull(isoforms_for_qpcr[row, ], ensembl_transcript_id)
    gn <- unique(pull(filter(exons, ensembl_transcript_id == tx), ensembl_gene_id))
    
    print(tx)
    
    my_exons <- exons %>%
        filter(ensembl_transcript_id == tx) %>%
        pull(ensembl_exon_id)
    other_tx <- annotations %>%
        filter(ensembl_gene_id == gn) %>%
        pull(ensembl_transcript_id) %>%
        unique()
    if (length(other_tx) < 1) { next }
    for (my_exon in my_exons) {
        my_exon_start <- exons %>%
            filter(ensembl_exon_id == my_exon) %>%
            pull(exon_chrom_start) %>%
            unique()
        my_exon_end <- exons %>%
            filter(ensembl_exon_id == my_exon) %>%
            pull(exon_chrom_end) %>%
            unique()
        message(my_exon_start, ":", my_exon_end)
        found_exons <- exons %>%
            filter(ensembl_exon_id != my_exon) %>%
            filter(
                (exon_chrom_start <= my_exon_start & my_exon_start <= exon_chrom_end) |
                    (exon_chrom_start <= my_exon_end & my_exon_end <= exon_chrom_end)
            ) %>%
            nrow()
        if (found_exons == 0) {
            message("Found transcript: ", tx)
            tx_unique_region <- c(tx_unique_region, tx)
            ex_unique_region <- c(ex_unique_region, my_exon)
            ln_unique_region <- c(ln_unique_region, my_exon_end - my_exon_start)
            break
        }
    }
}
unique_regions <- tibble(
    ensembl_transcript_id = tx_unique_region,
    ensembl_exon_id = ex_unique_region,
    unique_region_length = ln_unique_region
)
exon_seq <- getBM(
    attributes = c("ensembl_exon_id", "gene_exon", "exon_chrom_start", "exon_chrom_end", "rank"),
    filters = c("ensembl_exon_id"),
    mart = mart,
    values = unique_regions$ensembl_exon_id
)
unique_regions_seq <- left_join(unique_regions, exon_seq, by = "ensembl_exon_id") %>%
    left_join(isoforms_for_qpcr, by = "ensembl_transcript_id") %>%
    filter(unique_region_length > 200)
write_csv(unique_regions_seq, "data/TranscriptsWithUniqueExons.csv")

unique_regions_seq <- read_csv("data/TranscriptsWithUniqueExons.csv") %>%
    mutate(female_22pcw_frontal_lobe_tpm = female_22pcw_frontal_lobe_tpm[ensembl_transcript_id])
