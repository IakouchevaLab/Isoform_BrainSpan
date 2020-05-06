library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

PCC_THRESHOLD <- 0.9

iexpr <- readRDS("data/RegressIsoformCounts.rds")
metadata <- read_tsv("data/metadata.tsv")

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt",
    header = TRUE, sep = ",", row.names = 1,
    stringsAsFactors = FALSE
)
rbp_mele <- c(
    "HNRNPD", "SRSF4", "HNRNPA0", "KHDRBS1", "KHSRP", "SRSF7", "ZRANB2", 
    "SRSF11", "RBM25", "HNRNPH1", "SRSF1", "SFPQ", "HNRNPH3", "TRA2A", "SRSF6", 
    "RBFOX2", "SRSF2", "SF1", "RBM5", "SF3B1", "HNRNPA3", "TARDBP", "PCBP2", 
    "HNRNPM", "HNRNPA1", "HNRNPU", "SRSF3", "SRSF9", "HNRNPC", "ELAVL3", "QKI", 
    "TIA1", "SRRM1", "TIAL1", "CELF1", "PTBP1", "SYNCRIP", "MBNL1", "ELAVL1", 
    "PTBP2", "FMR1", "TRA2B", "RBMX", "DAZAP1", "HNRNPF", "CELF2", "RBFOX1", 
    "HNRNPK", "HNRNPL", "HNRPDL", "SRSF5", "FUS", "YBX1", "PCBP1", "HNRNPA2B1", 
    "SRSF10", "HNRPLL", "NOVA1", "KHDRBS3", "NOVA2", "ELAVL2", "ELAVL4", 
    "HNRNPH2", "KHDRBS2", "RBM4", "ESRP1", "ESRP2"
)

rbp_ensembl <- annotations %>%
    filter(external_gene_id %in% rbp_mele) %>%
    filter(ensembl_transcript_id %in% rownames(iexpr))

variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
) %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(`Affected status` == 2) %>% 
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained", 
            "splice_donor_variant", "splice_acceptor_variant"
        )
    )

all_observed_counts <- tibble()
for (period in c(sort(unique(metadata$Period)), "Prenatal", "Postnatal")) {
    samples <- metadata %>%
        filter(
            if (period == "Prenatal") {
                Prenatal
            } else if (period == "Postnatal") {
                !Prenatal
            } else {
                Period == period
            }
        ) %>%
        pull(Sample)
    # Filter for expression data pertaining to only these samples
    this_expr <- iexpr[, samples]
    # Select only transcripts that are expressed in this sample set
    expressed_transcripts <- rownames(this_expr)[apply(this_expr, 1, function(x) sum(x >= 0.1) >= 0.8 * length(samples))]
    
    impact_isoforms <- intersect(unique(pull(variants, `Ensembl Transcript ID`)), expressed_transcripts)
    
    observed <- expand.grid(
        period = period,
        rbp_ensembl_transcript_id = rbp_ensembl$ensembl_transcript_id,
        target_ensembl_transcript_id = impact_isoforms,
        stringsAsFactors = FALSE
    ) %>%
        filter(rbp_ensembl_transcript_id != target_ensembl_transcript_id) %>%
        distinct() %>%
        mutate(
            rbp_ensembl_gene_id = rbp_ensembl$ensembl_gene_id[match(
                rbp_ensembl_transcript_id, rbp_ensembl$ensembl_transcript_id
            )]
        ) %>%
        mutate(
            rbp_gene_symbol = rbp_ensembl$external_gene_id[match(
                rbp_ensembl_gene_id, rbp_ensembl$ensembl_gene_id
            )]
        ) %>%
        mutate(
            pcc = mcmapply(
                function(r, t) {
                    cor(
                        this_expr[r, samples], 
                        this_expr[t, samples], 
                        method = "pearson"
                    )
                }, 
                rbp_ensembl_transcript_id, target_ensembl_transcript_id
            )
        ) %>%
        group_by(period, rbp_ensembl_gene_id, target_ensembl_transcript_id) %>%
        summarize(pcc_max = max(pcc)) %>%
        filter(pcc_max >= PCC_THRESHOLD) %>%
        count(period, rbp_ensembl_gene_id)
    
    all_observed_counts <- bind_rows(all_observed_counts, observed)
}