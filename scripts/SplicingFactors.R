library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

iexpr <- readRDS("data/RegressIsoformCounts.rds")
gexpr <- readRDS("data/RegressGeneCounts.rds")
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
impact_isoforms <- intersect(unique(pull(variants, `Ensembl Transcript ID`)), rownames(iexpr))
impact_genes <- unique(
    annotations$ensembl_gene_id[match(impact_isoforms, annotations$ensembl_transcript_id)]
)
nonimpact_isoforms <- annotations %>%
    filter(! ensembl_transcript_id %in% impact_isoforms) %>%
    filter(ensembl_transcript_id %in% rownames(iexpr)) %>%
    mutate(
        ln_min = transcript_length * 0.9,
        ln_max = transcript_length * 1.1,
        gc_min = percentage_gc_content * 0.9,
        gc_max = percentage_gc_content * 1.1
    ) %>%
    distinct() %>% {
        mclapply(
            setNames(nm = impact_isoforms),
            function(tx) {
                this_ln <- annotations$transcript_length[annotations$ensembl_transcript_id == tx]
                this_gc <- annotations$percentage_gc_content[annotations$ensembl_transcript_id == tx]
                neighbors <- filter(., ln_min <= this_ln & this_ln <= ln_max) %>%
                    filter(gc_min <= this_gc & this_gc <= gc_max) %>%
                    pull(ensembl_transcript_id)
                if (length(neighbors) == 0) {
                    return(NA)
                } else {
                    return(neighbors)
                }
            }
        )
    }

PCC_THRESHOLD <- 0.90

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
    observed <- as_tibble(expand.grid(
        period = period,
        rbp_ensembl_transcript_id = rbp_ensembl$ensembl_transcript_id,
        target_ensembl_transcript_id = impact_isoforms,
        stringsAsFactors = FALSE
    )) %>%
        distinct() %>%
        mutate(
            rbp_ensembl_gene_id = rbp_ensembl$ensembl_gene_id[match(
                rbp_ensembl_transcript_id, rbp_ensembl$ensembl_transcript_id
            )]
        ) %>%
        mutate(
            pcc = mcmapply(
                function(r, t) {
                    cor(iexpr[r, samples], iexpr[t, samples], method = "pearson")
                }, rbp_ensembl_transcript_id, target_ensembl_transcript_id
            )
        )
    observed_counts <- observed %>%
        filter(rbp_ensembl_transcript_id != target_ensembl_transcript_id) %>%
        group_by(period, rbp_ensembl_gene_id, target_ensembl_transcript_id) %>%
        summarize(pcc_max = max(pcc)) %>%
        filter(pcc_max >= PCC_THRESHOLD) %>%
        count(period, rbp_ensembl_gene_id) %>%
        rename(observed_n = n)

    neighbors_permutations <- lapply(
        1:10,
        function(x) {
            sapply(nonimpact_isoforms, sample, size = 1, USE.NAMES = FALSE)
        }
    )

    expected_counts <- expand.grid(
        period = period,
        rbp_ensembl_transcript_id = rbp_ensembl$ensembl_transcript_id,
        target_ensembl_transcript_id = neighbors_permutations[[i]],
        permutation = i,
        stringsAsFactors = FALSE
    ) %>%
        as_tibble() %>%
        distinct() %>%
        mutate(
            rbp_ensembl_gene_id = rbp_ensembl$ensembl_gene_id[match(
                rbp_ensembl_transcript_id, rbp_ensembl$ensembl_transcript_id
            )]
        ) %>%
        mutate(
            pcc = mcmapply(
                function(r, t) {
                    cor(iexpr[r, samples], iexpr[t, samples], method = "pearson")
                }, rbp_ensembl_transcript_id, target_ensembl_transcript_id
            )
        ) %>%
        filter(rbp_ensembl_transcript_id != target_ensembl_transcript_id) %>%
        group_by(period, rbp_ensembl_gene_id, target_ensembl_transcript_id) %>%
        summarize(pcc_max = max(pcc)) %>%
        filter(pcc_max >= PCC_THRESHOLD) %>%
        count(period, rbp_ensembl_gene_id) %>%
        rename(expected_n = n) %>%
        inner_join(observed_counts, by = c("period", "rbp_ensembl_gene_id")) %>%
        mutate(Permutation = i)
    saveRDS(
        expected_counts,
        paste0(
            "data/rbp/rbp_expected_Period", period, "_",
            i, "_Threshold", PCC_THRESHOLD, ".rds"
        )
    )
}

# rbp_permutations <- list.files("data/rbp/", full.names = TRUE)
# expected_counts <- readRDS(rbp_permutations[[1]])


