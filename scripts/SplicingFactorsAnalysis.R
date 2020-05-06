library(tidyverse)

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
    filter(external_gene_id %in% rbp_mele)

rbp_results_files <- list.files("data/rbp", pattern = "*Threshold0.9", full.names = TRUE)
metadata <- read_tsv("data/metadata.tsv")

rbp_results_table <- expand.grid(
    Period = c(sort(unique(metadata$Period)), "Prenatal", "Postnatal"),
    Permutation = 1:100
) %>%
    mutate(
        Filename = paste0(
            "data/rbp/rbp_expected_Period", Period, "_", Permutation, ".rds"
        )
    )

all(sapply(rbp_results_table$Filename, file.exists))

rbp_results_combined <- bind_rows(lapply(rbp_results_files, readRDS))

permutation_analysis <- expand.grid(
    Period = unique(rbp_results_combined$period),
    rbp_ensembl_gene_id = unique(rbp_results_combined$rbp_ensembl_gene_id)
) %>%
    mutate(
        P = mapply(
            function(p, g) {
                res <- filter(rbp_results_combined, period == p & rbp_ensembl_gene_id == g)
                exp_n <- res$expected_n
                obs_n <- res$observed_n
                return(((sum(exp_n >= c(obs_n)) + 1) / 101))
            }, Period, rbp_ensembl_gene_id
        )
    ) %>%
    mutate(fdr = p.adjust(P, method = "BH")) %>%
    left_join(
        distinct(dplyr::select(rbp_ensembl, ensembl_gene_id, external_gene_id)),
        by = c("rbp_ensembl_gene_id" = "ensembl_gene_id")
    ) %>%
    mutate(Period = factor(Period, levels = rev(c(as.character(seq(2, 13)), "Prenatal", "Postnatal")))) %>%
    mutate(LFDR = -log10(fdr))

rbp_correlation_significance <- ggplot() +
    geom_tile(
        data = filter(permutation_analysis, fdr <= 0.05),
        mapping = aes(
            x = external_gene_id, y = Period, fill = fdr
        ),
        colour = "black"
    ) +
    geom_tile(
        data = filter(permutation_analysis, fdr > 0.05),
        mapping = aes(
            x = external_gene_id, y = Period
        ),
        colour = "black",
        fill = "white"
    ) +
    # scale_fill_gradientn(
    #     values = c(0, 0.01, 0.05, 1), colours = c("red", "white", "white"),
    #     limits = c(0, 1)
    # ) +
    scale_fill_gradient(low = "red", high = "white") +
    theme_bw() +
    theme(
        text = element_text(size = 7),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title = element_blank(),
        legend.position = "top"
    )
rbp_correlation_significance

ggsave(
    filename = "data/figures/RBPSignificance_Threshold0.9.pdf",
    plot = rbp_correlation_significance,
    width = 12, height = 5.5, units = "cm", device = "pdf",
    useDingbats = FALSE
)


tmp <- data.frame(A = c(1,2,3,4,5), B = c(0,0,4,5,6), C = c("A", "A", "A", "B", "B")) %>% 
    group_by(C) %>%
    summarise(gt = sum(A > B)) %>%
    print()





test_rbp_permutations <- list.files("data/rbp/", pattern = "*Threshold*", full.names = TRUE) %>%
    lapply(
        ., function(fn) {
            readRDS(fn)
        }
    ) %>%
    bind_rows() %>% 
    group_by(period, rbp_ensembl_gene_id) %>%
    summarize(o_gt_e = sum(observed_n > expected_n)) %>%
    mutate(pvalue = 1 - ((o_gt_e + 1) / 111)) %>%
    mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
    left_join(
        distinct(dplyr::select(rbp_ensembl, ensembl_gene_id, external_gene_id)),
        by = c("rbp_ensembl_gene_id" = "ensembl_gene_id")
    ) %>%
    mutate(selector = paste(period, external_gene_id)) %>% {
        s <- .$selector
        bind_rows(
            .,
            expand.grid(period = .$period, external_gene_id = .$external_gene_id) %>%
                mutate(fdr = 1) %>%
                mutate(selector = paste(period, external_gene_id)) %>%
                filter(! selector %in% s)
        )
    } %>% 
    mutate(Period = factor(period, levels = rev(c(as.character(seq(2, 13)), "Prenatal", "Postnatal"))))
test_rbp_permutations_plot <- test_rbp_permutations %>% {
        ggplot() +
            geom_tile(
                data = filter(., fdr > 0.05),
                mapping = aes(
                    x = external_gene_id, y = Period
                ),
                colour = "black",
                fill = "white"
            ) +
            geom_tile(
                data = .,
                mapping = aes(
                    x = external_gene_id, y = Period, fill = fdr
                ),
                colour = "black"
            ) +
            scale_fill_gradient(low = "red", high = "white") +
            theme_bw() +
            theme(
                text = element_text(size = 7),
                axis.text.x = element_text(angle = 60, hjust = 1),
                axis.title = element_blank(),
                legend.position = "top"
            )
    }
ggsave(
    filename = "data/figures/RBPSignificance_Threshold09.pdf",
    plot = test_rbp_permutations,
    width = 12, height = 5.5, units = "cm", device = "pdf",
    useDingbats = FALSE
)
