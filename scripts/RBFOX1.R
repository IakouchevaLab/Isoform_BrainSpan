library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

iexpr <- readRDS("data/RegressIsoformCounts.rds")
metadata <- read_tsv("data/metadata.tsv")

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt",
    header = TRUE, sep = ",", row.names = 1,
    stringsAsFactors = FALSE
)
rbfox1_ensembl <- annotations %>%
    filter(external_gene_id == "RBFOX1")
# mart <- biomaRt::useMart(biomart = "ensembl",
#                          dataset = "hsapiens_gene_ensembl",
#                          host = "GRCh37.ensembl.org")
# Exon 7743317-7743369
# biomartCacheClear()
# rbfox1_ensembl <- biomaRt::getBM(
#     attributes = c("ensembl_gene_id", "ensembl_transcript_id",
#                    "external_gene_name", "ensembl_exon_id", "rank",
#                    "exon_chrom_start", "exon_chrom_end"),
#     filters = c("ensembl_transcript_id"),
#     values = annotations$ensembl_transcript_id[annotations$external_gene_id == "RBFOX1"],
#     mart = mart
# ) %>%
#     mutate(exon_ranges = paste(exon_chrom_start, exon_chrom_end, sep = '-')) %>%
#     group_by(ensembl_gene_id, ensembl_transcript_id, external_gene_name) %>%
#     mutate(exon_ranges = paste(exon_ranges, collapse = ' ')) %>%
#     distinct(ensembl_gene_id, ensembl_transcript_id, external_gene_name, exon_ranges) %>%
#     mutate(Localization = sapply(
#         exon_ranges, 
#         function(er) {
#             if(grepl("7743317-7743369", er)) { 
#                 return("cytoplasmic") 
#             } else { 
#                 return("nuclear") 
#             }
#         }
#     ))

cytoplasmic_all <- c("ENST00000355637", "ENST00000535565", "ENST00000547372", 
                     "ENST00000552089")
nuclear_all <- c("ENST00000311745", "ENST00000340209", "ENST00000422070",
                 "ENST00000436368", "ENST00000547338", "ENST00000547427",
                 "ENST00000547605", "ENST00000548749", "ENST00000550418", 
                 "ENST00000551752", "ENST00000553186", "ENST00000564850",
                 "ENST00000567470", "ENST00000569889", "ENST00000570188", 
                 "ENST00000570626")

cytoplasmic <- intersect(cytoplasmic_all, rownames(iexpr))
nuclear <- intersect(nuclear_all, rownames(iexpr))

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

periods <- c(sort(unique(metadata$Period)), "Prenatal", "Postnatal")

empirical_data <- foreach(i = seq(1, length(periods))) %dopar% {
    period <- periods[[i]]
    samples <- metadata %>%
        filter(
            if (period == "Prenatal") {
                Prenatal
            } else if (period == "Postnatal") {
                !Prenatal
            } else {
                Period == as.numeric(period)
            }
        ) %>%
        pull(Sample)
    # Filter for expression data pertaining to only these samples
    this_expr <- iexpr[, samples]
    # Select only transcripts that are expressed in this sample set
    expressed_transcripts <- rownames(this_expr)[apply(this_expr, 1, function(x) sum(x >= 0.1) >= 0.8 * length(samples))]
    impact_isoforms <- intersect(unique(pull(variants, `Ensembl Transcript ID`)), expressed_transcripts)
    this_dat <- bind_rows(
        expand.grid(
            rbfox1 = cytoplasmic,
            target = sapply(impact_isoforms, sample, size = 1),
            localization = "cytoplasmic"
        ),
        expand.grid(
            rbfox1 = nuclear,
            target = sapply(impact_isoforms, sample, size = 1),
            localization = "nuclear"
        )
    ) %>%
        mutate(cor = mapply(
            function(src, tar) {
                cor(
                    this_expr[src, samples], this_expr[tar, samples],
                    method = "pearson"
                )
            }, rbfox1, target
        )) %>%
        mutate(Period = period)
} %>%
    bind_rows()

grouped_empirical_data <- empirical_data %>%
    group_by(target, localization, Period) %>%
    summarise(mean_cor = mean(cor))

empirical_data_counts <- empirical_data %>%
    mutate(intervals = cut(cor, breaks = seq(-1, 1, by = 0.1))) %>%
    group_by(localization, Period) %>%
    count(localization, Period, intervals) %>% {
        ggplot(
            data = .,
            mapping = aes(
                x = intervals, y = n, colour = localization
            )
        ) +
            facet_wrap(~ Period) +
            geom_point(size = 5) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 30, hjust = 1)
            )
    }

nuclear_cor <- arrange(filter(grouped_empirical_data, localization == "nuclear"), Period, target)
cytoplasmic_cor <- arrange(filter(grouped_empirical_data, localization == "cytoplasmic"), Period, target)
empirical_cor_tests <- lapply(periods, function(p) {
    wilcox.test(
        nuclear_cor$mean_cor[nuclear_cor$Period == p],
        cytoplasmic_cor$mean_cor[cytoplasmic_cor$Period == p],
        alternative = "g"
    ) %>%
        broom::tidy() %>%
        mutate(Period = p)
}) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))
