library(tidyverse)

asd <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 6.xlsx", sheet = 2
) %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(`Affected status` == 2) %>% 
    filter(
        Consequence %in% c(
            "frameshift_variant", "start_lost", "stop_gained", 
            "splice_donor_variant", "splice_acceptor_variant"
        )
    )
impacted <- variants$`Ensembl Transcript ID`

deg <- lapply(
    2:12,
    function(x) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 2.xlsx", sheet = x
        )
    }
) %>%
    bind_rows() %>%
    mutate(DE = abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
    mutate(id = paste0(`Ensembl Gene ID`, Contrast))

dei <- lapply(
    2:12,
    function(x) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 3.xlsx", sheet = x
        )
    }
) %>%
    bind_rows() %>%
    mutate(ASD = `Gene Symbol` %in% asd) %>%
    mutate(DE = abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
    mutate(id = paste0(`Ensembl Gene ID`, Contrast)) %>%
    mutate(Specific = id %in% pull(filter(deg, DE), id))

de_impact_spec <- lapply(
    unique(dei$Contrast),
    function(contr) {
        all_isoforms <- dei %>%
            filter(Contrast == contr & Specific) %>%
            pull(`Ensembl Transcript ID`)
        de_isoforms <- dei %>%
            filter(DE & Contrast == contr & Specific) %>%
            pull(`Ensembl Transcript ID`)
        this_impacted <- intersect(impacted, all_isoforms)
        fisher.test(
            matrix(c(
                length(intersect(de_isoforms, this_impacted)),
                length(this_impacted[!this_impacted %in% de_isoforms]),
                length(de_isoforms[!de_isoforms %in% this_impacted]),
                length(all_isoforms[!all_isoforms %in% c(de_isoforms, this_impacted)])
            ), ncol = 2, nrow = 2),
            alternative = "greater"
        ) %>%
            broom::tidy() %>%
            mutate(Contrast = contr) %>%
            mutate(SpecificDECount = length(de_isoforms))
    }
) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))
ggplot(
    data = de_impact_spec, mapping = aes(x = SpecificDECount, y = -log10(fdr))
) +
    geom_line()

asd_de_impact_spec <- lapply(
    unique(dei$Contrast),
    function(contr) {
        all_isoforms <- dei %>%
            filter(Contrast == contr & Specific & ASD) %>%
            pull(`Ensembl Transcript ID`)
        de_isoforms <- dei %>%
            filter(DE & Contrast == contr & Specific & ASD) %>%
            pull(`Ensembl Transcript ID`)
        fisher.test(
            matrix(c(
                length(intersect(de_isoforms, impacted)),
                length(impacted[!impacted %in% de_isoforms]),
                length(de_isoforms[!de_isoforms %in% impacted]),
                length(all_isoforms[!all_isoforms %in% c(de_isoforms, impacted)])
            ), ncol = 2, nrow = 2),
            alternative = "greater"
        ) %>%
            broom::tidy() %>%
            mutate(Contrast = contr) %>%
            mutate(SpecificASDDECount = length(de_isoforms))
    }
) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))

asd_de_impact <- lapply(
    unique(dei$Contrast),
    function(contr) {
        all_isoforms <- dei %>%
            filter(Contrast == contr & ASD) %>%
            pull(`Ensembl Transcript ID`)
        de_isoforms <- dei %>%
            filter(DE & Contrast == contr & ASD) %>%
            pull(`Ensembl Transcript ID`)
        fisher.test(
            matrix(c(
                length(intersect(de_isoforms, impacted)),
                length(impacted[!impacted %in% de_isoforms]),
                length(de_isoforms[!de_isoforms %in% impacted]),
                length(all_isoforms[!all_isoforms %in% c(de_isoforms, impacted)])
            ), ncol = 2, nrow = 2),
            alternative = "greater"
        ) %>%
            broom::tidy() %>%
            mutate(Contrast = contr) %>%
            mutate(SpecificASDDECount = length(de_isoforms))
    }
) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))
