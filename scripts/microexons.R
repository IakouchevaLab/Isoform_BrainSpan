library(tidyverse)
library(biomaRt)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)
set.seed(1)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", 
    header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE
)
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
impacted_isoforms <- unique(variants$`Ensembl Transcript ID`)
impacted_genes <- unique(variants$`Ensembl Gene ID`)

mart <- useMart(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    host = "GRCh37.ensembl.org"
)
all_isoforms <- rownames(readRDS("data/RegressIsoformCounts.rds"))
all_genes <- rownames(readRDS("data/RegressGeneCounts.rds"))
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
# microexons_pos <- microexons %>%
#     filter(microexon) %>%
#     arrange(ensembl_gene_id, ensembl_transcript_id) %>%
#     left_join(gene_iso_count, by = "ensembl_gene_id") %>%
#     count(
#         ensembl_gene_id, percentage_gene_gc_content,
#         exon_chrom_start, exon_chrom_end, microexon, number_isoforms,
#         name = "number_microexon"
#     ) %>%
#     mutate(altRegMicro = number_isoforms != number_microexon)
# microexons_distinct <- microexons %>% 
#     distinct(
#         ensembl_gene_id, ensembl_transcript_id, 
#         transcript_length, percentage_gene_gc_content, microexon
#     ) %>% {
#         m <- filter(., microexon)
#         n <- filter(
#             ., !microexon & ! ensembl_transcript_id %in% m$ensembl_transcript_id
#         )
#         bind_rows(m, n)
#     } %>%
#     distinct()

microexon_isoforms <- microexons %>%
    filter(microexon) %>%
    pull(ensembl_transcript_id)

asd <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

impact_mx <- tibble(
    ensembl_gene_id = all_genes
) %>%
    left_join(
        annotations[, c("ensembl_gene_id", "ensembl_transcript_id")],
        by = "ensembl_gene_id"
    ) %>%
    mutate(ASD = ensembl_transcript_id %in% asd) %>%
    filter(ensembl_transcript_id %in% all_isoforms) %>%
    mutate(Impact = ensembl_transcript_id %in% impacted_isoforms) %>%
    mutate(GeneImpact = ensembl_gene_id %in% impacted_genes) %>%
    mutate(Microexon = ensembl_transcript_id %in% microexon_isoforms) %>%
    mutate(
        Background = ensembl_transcript_id %in% (
            all_exons %>%
                filter(exon_size <= 185 * 1.1 & exon_size >= 185 * 0.9) %>%
                pull(ensembl_transcript_id)
        )
    ) %>%
    left_join(
        annotations[, c("ensembl_transcript_id", "transcript_length", "percentage_gc_content")],
        by = "ensembl_transcript_id"
    )

# Impact v nonimpact transcript_length
wilcox.test(
    mean(filter(filter(impact_mx, GeneImpact), Impact)$transcript_length), 
    mean(filter(filter(impact_mx, GeneImpact), !Impact)$transcript_length),
    alternative = "greater"
)

# Permutation test impact v nonimpact
observed_microexons <- impact_mx %>%
    filter(Impact) %>% {
        mutate(
            ., 
            transcript_length_min = .$transcript_length * 0.85,
            transcript_length_max = .$transcript_length * 1.15,
            gc_min = .$percentage_gc_content * 0.85,
            gc_max = .$percentage_gc_content * 1.15
        )
    }
neighbor_sets <- lapply(
    setNames(nm = observed_microexons$ensembl_transcript_id),
    function(tx) {
        observed_microexons %>%
            filter(ensembl_transcript_id == tx) %>% {
                tx_min <- pull(., transcript_length_min)
                tx_max <- pull(., transcript_length_max)
                gc_min <- pull(., gc_min)
                gc_max <- pull(., gc_max)
                sample_set <- impact_mx %>%
                    filter(! Impact) %>%
                    filter(transcript_length <= tx_max & transcript_length >= tx_min) %>%
                    filter(percentage_gc_content <= gc_max & percentage_gc_content >= gc_min) %>%
                    pull(ensembl_transcript_id)
                if (length(sample_set) == 0) { return(NA) }
                else { return(sample_set) }
            }
    }
)
expected_microexons <- mclapply(
    1:1000,
    function(x) {
        null_isoforms <- sapply(neighbor_sets, sample, size = 1)
        tibble(ensembl_transcript_id = null_isoforms) %>%
            left_join(impact_mx, by = "ensembl_transcript_id") %>%
            mutate(Permutation = x) %>%
            return()
    }
) %>%
    bind_rows()
# saveRDS(expected_microexons, "data/expected_microexons.rds")
# expected_microexons <- readRDS("data/expected_microexons.rds")

microexon_permutation_test <- sapply(
    1:1000,
    function(x) {
        # i <- observed_microexons$Impact
        m <- sum(
            filter(expected_microexons, Permutation == x)$Microexon, 
            na.rm = TRUE
        )
        # sum(i == m) / nrow(filter(expected_microexons, Permutation == x))
    }
)

ggplot(
    data = data.frame(x = microexon_permutation_test),
    mapping = aes(x = x)
) +
    geom_histogram() +
    geom_vline(
        xintercept = sum(observed_microexons$Microexon[which(!is.na(neighbor_sets))])
    )

1 - (sum((sum(observed_microexons$Microexon) / nrow(observed_microexons)) > (microexon_permutation_test / nrow(observed_microexons))) / 1000)

rate_data <- tibble(
    Class = c("Impact", "Nonimpact"),
    MicroexonOccurrence = c(
        sum(observed_microexons$Microexon) / nrow(observed_microexons),
        mean(microexon_permutation_test) / nrow(observed_microexons)
    ),
    SDmin = c(
        NA, 
        mean(microexon_permutation_test / nrow(observed_microexons)) - sd(microexon_permutation_test / nrow(observed_microexons))
    ),
    SDmax = c(
        NA, 
        mean(microexon_permutation_test / nrow(observed_microexons)) + sd(microexon_permutation_test / nrow(observed_microexons))
    )
)

ggplot(
    data = rate_data,
    mapping = aes(
        x = Class, y = MicroexonOccurrence,
        ymin = SDmin, ymax = SDmax,
        fill = Class
    )
) +
    geom_bar(stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(width = 0.25) +
    scale_fill_manual(
        values = c("Impact" = "steelblue", "Nonimpact" = "firebrick")
    ) +
    theme_minimal()




# Permutation test impact v background
neighbor_sets_bg <- lapply(
    setNames(nm = observed_microexons$ensembl_transcript_id),
    function(tx) {
        observed_microexons %>%
            filter(ensembl_transcript_id == tx) %>% {
                tx_min <- pull(., transcript_length_min)
                tx_max <- pull(., transcript_length_max)
                gc_min <- pull(., gc_min)
                gc_max <- pull(., gc_max)
                sample_set <- impact_mx %>%
                    filter(! Impact & ! GeneImpact) %>%
                    filter(transcript_length <= tx_max & transcript_length >= tx_min) %>%
                    filter(percentage_gc_content <= gc_max & percentage_gc_content >= gc_min) %>%
                    filter(Background) %>%
                    pull(ensembl_transcript_id)
                if (length(sample_set) == 0) { return(NA) }
                else { return(sample_set) }
            }
    }
)
expected_microexons_bg <- mclapply(
    1:1000,
    function(x) {
        null_isoforms <- sapply(neighbor_sets_bg, sample, size = 1)
        tibble(ensembl_transcript_id = null_isoforms) %>%
            left_join(impact_mx, by = "ensembl_transcript_id") %>%
            mutate(Permutation = x) %>%
            return()
    }
) %>%
    bind_rows()
metadata <- read_tsv("data/metadata.tsv")
iexpr <- readRDS("data/iso_tpm_filter.rds")
test_expression <- lapply(
    unique(metadata$Period),
    function(period) {
        samples <- metadata$Sample[metadata$Period == period]
        impact_expression <- impact_mx %>%
            filter(Impact) %>%
            pull(ensembl_transcript_id) %>% {
                mean(iexpr[., samples])
            }
        impact_mx_expression <- impact_mx %>%
            filter(Impact & Microexon) %>%
            pull(ensembl_transcript_id) %>% {
                mean(iexpr[., samples])
            }
        nonimpact_expressions <- sapply(
            1:1000,
            function(x) {
                expected_microexons %>%
                    pull(ensembl_transcript_id) %>% {
                        mean(iexpr[., samples])
                    }
            }
        )
        nonimpact_mx_expressions <- sapply(
            1:1000,
            function(x) {
                expected_microexons %>%
                    filter(Microexon) %>%
                    pull(ensembl_transcript_id) %>% {
                        mean(iexpr[., samples])
                    }
            }
        )
        result <- tibble(
            impact_expression = impact_expression,
            impact_mx_expression = impact_mx_expression,
            nonimpact_expressions = nonimpact_expressions,
            nonimpact_mx_expressions = nonimpact_mx_expressions
        ) %>%
            nest(nonimpact_expressions, nonimpact_mx_expressions, .key = "nonimpact")
    }
)

microexon_permutation_test_bg <- sapply(
    1:1000,
    function(x) {
        # i <- observed_microexons$Impact
        m <- sum(
            filter(expected_microexons_bg, Permutation == x)$Microexon, 
            na.rm = TRUE
        )
        # sum(i == m) / nrow(filter(expected_microexons, Permutation == x))
    }
)

ggplot(
    data = data.frame(x = microexon_permutation_test_bg),
    mapping = aes(x = x)
) +
    geom_histogram() +
    geom_vline(
        xintercept = sum(observed_microexons$Microexon[which(!is.na(neighbor_sets_bg))])
    )

1 - (sum((sum(observed_microexons$Microexon) / nrow(observed_microexons)) > (microexon_permutation_test_bg / nrow(observed_microexons))) / 1000)

rate_data_bg <- tibble(
    Class = c("Impact", "Nonimpact"),
    MicroexonOccurrence = c(
        sum(observed_microexons$Microexon) / nrow(observed_microexons),
        mean(microexon_permutation_test_bg) / nrow(observed_microexons)
    ),
    SDmin = c(
        NA, 
        mean(microexon_permutation_test / nrow(observed_microexons)) - sd(microexon_permutation_test_bg / nrow(observed_microexons))
    ),
    SDmax = c(
        NA, 
        mean(microexon_permutation_test / nrow(observed_microexons)) + sd(microexon_permutation_test_bg / nrow(observed_microexons))
    )
)

(
    ggplot(
        data = rate_data_bg,
        mapping = aes(
            x = Class, y = MicroexonOccurrence,
            ymin = SDmin, ymax = SDmax,
            fill = Class
        )
    ) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_errorbar(width = 0.25) +
        scale_fill_manual(
            values = c("Impact" = "steelblue", "Nonimpact" = "firebrick")
        ) +
        theme_bw() +
        theme(
            text = element_text(size = 5),
            legend.position = "none"
        )
) %>% {
    ggsave(
        filename = "data/figures/microexon_bg.pdf",
        plot = .,
        device = "pdf", width = 2.5, height = 3, units = "cm"
    )
}




# Alt reg genes
gene_alt_mx_impact <- impact_mx %>% {
        df <- mutate(., counter = 1) %>%
            group_by(ensembl_gene_id) %>%
            mutate(num_tx = sum(counter)) %>%
            mutate(num_impact_tx = sum(Impact)) %>%
            mutate(num_microexon_tx = sum(Microexon)) %>%
            distinct(
                ensembl_gene_id, ASD, 
                GeneImpact, num_tx, num_microexon_tx,
                percentage_gc_content
            ) %>%
            mutate(
                AltRegMicroexon = num_tx != num_microexon_tx &
                    num_microexon_tx > 0
            )
        return(df)
    } %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "start_position", "end_position"),
            filters = "ensembl_gene_id",
            values = .$ensembl_gene_id,
            mart = mart
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(gene_length = abs(start_position - end_position)) %>%
    mutate(
        ln_min = gene_length * 0.9,
        ln_max = gene_length * 1.1,
        gc_min = percentage_gc_content * 0.9,
        gc_max = percentage_gc_content * 1.1
    )

gene_mask <- gene_alt_mx_impact %>%
    dplyr::select(ensembl_gene_id) %>%
    left_join(
        getBM(
            attributes = c(
                "ensembl_gene_id", "exon_chrom_start", "exon_chrom_end"
            ),
            filters = "ensembl_gene_id",
            values = .$ensembl_gene_id,
            mart = mart
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(exon_size = abs(exon_chrom_start - exon_chrom_end))
average_genemask_exon_size_nonimpact <- gene_mask %>% {
    genes <- filter(gene_alt_mx_impact, !GeneImpact)$ensembl_gene_id
    filter(., ensembl_gene_id %in% genes)
} %>%
    pull(exon_size) %>%
    mean()
average_genemask_exon_size_impact <- gene_mask %>% {
    genes <- filter(gene_alt_mx_impact, GeneImpact)$ensembl_gene_id
    filter(., ensembl_gene_id %in% genes)
} %>%
    pull(exon_size) %>%
    mean()
    
gene_neighbors <- mcmapply(
    function(this_ln_min, this_ln_max, this_gc_min, this_gc_max, ensg) {
        print(ensg)
        gene_alt_mx_impact %>%
            filter(gene_length >= this_ln_min & gene_length <= this_ln_max) %>%
            filter(percentage_gc_content >= this_gc_min & percentage_gc_content <= this_gc_max) %>%
            pull(ensembl_gene_id) %>% {
                if (length(.) < 1) {
                    return(NA)
                } else {
                    return(.)
                }
            }
    },
    gene_alt_mx_impact$ln_min, gene_alt_mx_impact$ln_max,
    gene_alt_mx_impact$gc_min, gene_alt_mx_impact$gc_max,
    gene_alt_mx_impact$ensembl_gene_id
) %>%
    setNames(nm = gene_alt_mx_impact$ensembl_gene_id)

# impact and alt mx vs nonimpact and alt mx
null_permutations <- gene_alt_mx_impact %>% {
    impact_and_alt_mx <- filter(., GeneImpact) %>%
        pull(ensembl_gene_id)
    print(length(impact_and_alt_mx))
    potential_neighbors <- gene_neighbors[impact_and_alt_mx]
    permutations <- mclapply(
        1:1000,
        function(x) {
            n <- sapply(potential_neighbors, sample, size = 1)
            sum(pull(filter(gene_alt_mx_impact, ensembl_gene_id %in% n), AltRegMicroexon))
        }
    )
}
actual_impact <- filter(gene_alt_mx_impact, GeneImpact & AltRegMicroexon) %>%
    nrow()

1 - (sum(actual_impact > null_permutations) / 1000)

data.frame(
    x = c("ASD LoF", "Non-impacted by ASD LoF"),
    y = c(actual_impact / 824, mean(as.numeric(null_permutations)) / 824),
    ymin = c(
        NA,
        mean(as.numeric(null_permutations) / 824) - sd(as.numeric(null_permutations) / 824)
    ),
    ymax = c(
        NA,
        mean(as.numeric(null_permutations) / 824) + sd(as.numeric(null_permutations) / 824)
    )
) %>%
    {
        ggplot(
            data = .,
            mapping = aes(
                x = x, y = y, ymin = ymin, ymax = ymax, fill = x
            )
        ) +
            geom_bar(stat = "identity", position = "dodge", colour = "black") +
            geom_errorbar(width = 0.25) +
            scale_fill_manual(
                values = c(
                    "ASD LoF" = "firebrick",
                    "Non-impacted by ASD LoF" = "steelblue"
                )
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 5),
                axis.text.x = element_text(angle = 30, hjust = 1),
                legend.position = "none"
            )
    } %>% {
        ggsave(
            filename = "data/figures/microexon_gene.pdf",
            plot = .,
            device = "pdf", width = 2.5, height = 3, units = "cm"
        )
    }
