library(tidyverse)
library(gprofiler2)

feature_anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]

module_assigns <- bind_rows(
    read_tsv("data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        mutate(Data = "Gene") %>%
        gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
        filter(kME_Key == paste0("kME", trimws(module_label))),
    read_tsv("data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        mutate(Data = "Isoform") %>%
        gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
        filter(kME_Key == paste0("kME", trimws(module_label)))
)

module_enrs <- module_assigns %>%
    distinct(Data, module_label, module_colour) %>%
    mutate(
        GOST = mapply(
            function(dat, ml, mc) {
                enr <- module_assigns %>%
                    filter(
                        Data == dat & module_label == ml & module_colour == mc
                    ) %>%
                    arrange(desc(kME)) %>%
                    pull(Feature) %>%
                    gost(
                        organism = "hsapiens",
                        ordered_query = TRUE,
                        correction_method = "fdr",
                        custom_bg = module_assigns$Feature,
                        source = c("GO:BP", "GO:MF")
                    )
                return(enr$result)
            }, Data, module_label, module_colour
        )
    )
module_enrs <- filter(module_enrs, !sapply(GOST, function(g) is.null(g))) %>%
    arrange(Data, module_label)

saveRDS(
    setNames(
        pull(filter(module_enrs, Data == "Gene"), GOST),
        paste(
            pull(filter(module_enrs, Data == "Gene"), module_label),
            pull(filter(module_enrs, Data == "Gene"), module_colour)
        )
    ),
    "data/WGCNA/GeneNetwork_DS2_MM20_GOST.rds"
)
writexl::write_xlsx(
    setNames(
        pull(filter(module_enrs, Data == "Gene"), GOST),
        paste(
            pull(filter(module_enrs, Data == "Gene"), module_label),
            pull(filter(module_enrs, Data == "Gene"), module_colour)
        )
    ),
    "data/WGCNA/GeneNetwork_DS2_MM20_GOST.xlsx"
)
saveRDS(
    setNames(
        pull(filter(module_enrs, Data == "Isoform"), GOST),
        paste(
            pull(filter(module_enrs, Data == "Isoform"), module_label),
            pull(filter(module_enrs, Data == "Isoform"), module_colour)
        )
    ),
    "data/WGCNA/IsoformNetwork_DS2_MM20_GOST.rds"
)
writexl::write_xlsx(
    setNames(
        pull(filter(module_enrs, Data == "Isoform"), GOST),
        paste(
            pull(filter(module_enrs, Data == "Isoform"), module_label),
            pull(filter(module_enrs, Data == "Isoform"), module_colour)
        )
    ),
    "data/WGCNA/IsoformNetwork_DS2_MM20_GOST.xlsx"
)


gene_module_enrs <- readRDS(
    "data/genes/Networks/Enrichments/Network_DS2_MM20_GO.rds"
)
isoform_module_enrs <- readRDS(
    "data/isoforms/Networks/Enrichments/Network_DS2_MM20_GO.rds"
)

gene_enr_plots <- lapply(
    setNames(nm = names(gene_module_enrs)),
    function(ml) {
        if (nrow(gene_module_enrs[[ml]]) < 1) {
            return(NULL)
        }
        module_colour <- filter(
            module_assigns$gene, 
            module_label == as.numeric(ml)
        ) %>%
            pull(module_colour) %>%
            unique()
        gene_module_enrs[[ml]] %>%
            filter(term.size <= 1000) %>%
            arrange(p.value) %>%
            mutate(Order = seq(nrow(.), 1)) %>%
            top_n(10, Order) %>%
            ggplot(
                data = .,
                mapping = aes(
                    x = -log10(p.value), 
                    y = reorder(str_wrap(term.name, 40), Order), 
                    size = overlap.size, shape = domain
                )
            ) +
                geom_point(fill = module_colour) +
                labs(
                    title = paste("Gene", ml, module_colour),
                    x = expression("-log"["10"]*"FDR")
                ) +
                scale_size_continuous(range = c(5, 15), name = "Overlap") +
                scale_shape_manual(
                    name = "Source",
                    values = c(
                        "BP" = 21,
                        "MF" = 22
                    )
                ) +
                guides(
                    shape = guide_legend(override.aes = list(size = 10))
                ) +
                theme_bw() +
                theme(
                    text = element_text(size = 24),
                    axis.title.y = element_blank()
                )
    }
)

isoform_enr_plots <- lapply(
    setNames(nm = names(isoform_module_enrs)),
    function(ml) {
        if (nrow(isoform_module_enrs[[ml]]) < 1) {
            return(NULL)
        }
        module_colour <- filter(
            module_assigns$isoform, 
            module_label == as.numeric(ml)
        ) %>%
            pull(module_colour) %>%
            unique()
        isoform_module_enrs[[ml]] %>%
            filter(term.size <= 1000) %>%
            arrange(p.value) %>%
            mutate(Order = seq(nrow(.), 1)) %>%
            top_n(10, Order) %>%
            ggplot(
                data = .,
                mapping = aes(
                    x = -log10(p.value), 
                    y = reorder(str_wrap(term.name, 40), Order), 
                    size = overlap.size, shape = domain
                )
            ) +
            geom_point(fill = module_colour) +
            labs(
                title = paste("Isoform", ml, module_colour),
                x = expression("-log"["10"]*"FDR")
            ) +
            scale_size_continuous(range = c(5, 15), name = "Overlap") +
            scale_shape_manual(
                name = "Source",
                values = c(
                    "BP" = 21,
                    "MF" = 22
                )
            ) +
            guides(
                shape = guide_legend(override.aes = list(size = 10))
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 24),
                axis.title.y = element_blank()
            )
    }
)

# Paired plots

gene_isoform_module_pairs <- data.frame(
    GeneModule =    c(1, 2, 3,  4,  5, 6,  7,  8),
    IsoformModule = c(1, 2, 4, 35, 11, 9, 23, 41)
)

gi_enr_plots <- apply(
    gene_isoform_module_pairs, 1, function(gi) {
        gml <- as.numeric(gi[["GeneModule"]])
        iml <- as.numeric(gi[["IsoformModule"]])
        gcol <- module_assigns$gene %>%
            filter(module_label == gml) %>%
            pull(module_colour) %>%
            unique()
        icol <- module_assigns$isoform %>%
            filter(module_label == iml) %>%
            pull(module_colour) %>%
            unique()
        if (nrow(gene_module_enrs[[as.character(gml)]]) == 0) {
            genr <- data.frame(p.value = Inf, term.name = "")
        } else {
            genr <- gene_module_enrs[[as.character(gml)]] %>%
                filter(term.size <= 1000) %>%
                mutate(Data = "Gene") %>%
                arrange(p.value) %>%
                mutate(Order = seq(nrow(.), 1)) %>%
                top_n(10, Order)
        }
        if (nrow(isoform_module_enrs[[as.character(iml)]]) == 0) {
            ienr <- data.frame(p.value = Inf, term.name = "")
        } else {
            ienr <- isoform_module_enrs[[as.character(iml)]] %>%
                filter(term.size <= 1000) %>%
                mutate(Data = "Isoform") %>%
                arrange(p.value) %>%
                mutate(Order = seq(nrow(.), 1)) %>%
                top_n(10, Order)
        }
        term_orders <- data.frame(
            term.name = unique(c(genr$term.name, ienr$term.name))
        ) %>%
            left_join(
                dplyr::select(genr, term.name, p.value), by = "term.name"
            ) %>%
            rename(gval = p.value) %>%
            left_join(
                dplyr::select(ienr, term.name, p.value), by = "term.name"
            ) %>%
                rename(ival = p.value) %>%
            mutate(
                val = mapply(function(x, y) min(x, y, na.rm = TRUE), gval, ival)
            ) %>%
            arrange(-val) %>%
            mutate(Order = seq(1, nrow(.)))
        ggplot(
            mapping = aes(
                x = -log10(p.value), 
                y = reorder(str_wrap(term.name, 40), Order.y), 
                size = overlap.size, shape = domain, colour = Data
            )
        ) +
            geom_point(
                data = genr %>%
                    full_join(term_orders, by = "term.name"), 
                fill = gcol, stroke = 3
            ) +
            geom_point(
                data = ienr %>%
                    full_join(term_orders, by = "term.name"), 
                fill = icol, stroke = 3
            ) +
            labs(
                title = paste(
                    paste("Gene", gml, gcol),
                    paste("Isoform", iml, icol),
                    sep = "\n"
                ),
                x = expression("-log"["10"]*"FDR")
            ) +
            scale_size_continuous(range = c(5, 15), name = "Overlap") +
            scale_shape_manual(
                name = "Source",
                values = c(
                    "BP" = 21,
                    "MF" = 22
                )
            ) +
            guides(
                shape = guide_legend(override.aes = list(size = 10))
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 24),
                axis.title.y = element_blank()
            )
    }
)
