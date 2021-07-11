library(tidyverse)
library(biomaRt)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)
set.seed(1)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", 
    header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE
)

mart <- useMart(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    host = "GRCh37.ensembl.org"
)
all_exons <- getBM(
    attributes = c(
        "ensembl_gene_id", "ensembl_transcript_id", 
        "transcript_length", "percentage_gene_gc_content"
    ),
    filters = c("ensembl_transcript_id"),
    values = rownames(readRDS("data/RegressIsoformCounts.rds")),
    mart = mart
) %>%
    left_join(
        getBM(
            attributes = c(
                "ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"
            ),
            filters = c("ensembl_transcript_id"),
            values = rownames(readRDS("data/RegressIsoformCounts.rds")),
            mart = mart
        ),
        by = "ensembl_transcript_id"
    ) %>%
    mutate(exon_size = abs(exon_chrom_start - exon_chrom_end))
microexons <- all_exons %>%
    mutate(microexon = exon_size <= 27 & exon_size >= 3)
variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
) %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained", 
            "splice_donor_variant", "splice_acceptor_variant"
        )
    ) %>%
    mutate(
        microexon = mapply(
            function(s, e) {
                (microexons %>%
                    filter(microexon) %>%
                    filter(
                        (s <= exon_chrom_end & s >= exon_chrom_start) |
                            (e <= exon_chrom_end & e >= exon_chrom_start)
                    ) %>%
                    nrow()) > 0
            },
            `Variant start`, `Variant end`
        )
    )

fisher.test(
    x = matrix(
        c(
            nrow(variants %>% filter(`Affected status` == 2 & microexon)),
            nrow(variants %>% filter(`Affected status` == 1 & microexon)),
            nrow(variants %>% filter(`Affected status` == 2 & ! microexon)),
            nrow(variants %>% filter(`Affected status` == 1 & ! microexon))
        ),
        nrow = 2
    )
)
