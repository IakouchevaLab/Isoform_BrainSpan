library(tidyverse)
library(WGCNA)

gene_cts <- readRDS("data/RegressGeneCounts.rds")
iso_cts <- readRDS("data/RegressIsoformCounts.rds")

iso_cts

gene_network <- readRDS("data/genes/Networks/Network_DS2_MM20.rds")
isoform_network <- readRDS("data/isoforms/Networks/Network_DS2_MM20.rds")

module_assignments <- read_tsv("data/ModuleInfo.tsv")

gene_module_assignments <- module_assignments %>%
    filter(Data == "Gene") %>%
    distinct(module_label, ensembl_gene_id)
gene_colors <- setNames(
    pull(gene_module_assignments, module_label),
    pull(gene_module_assignments, ensembl_gene_id)
) %>% as.list()

iso_gene_module_assignments <- module_assignments %>%
    filter(Data == "Isoform") %>%
    distinct(module_label, ensembl_gene_id)
iso_gene_colors <- setNames(
    pull(iso_gene_module_assignments, module_label),
    pull(iso_gene_module_assignments, ensembl_gene_id)
) %>% as.list()


gene_input_mp <- fixDataStructure(t(gene_cts))
isoform_input_mp <- fixDataStructure(t(iso_cts))
set_labels <- c("gene", "isoform")
multi_expr <- list(gene = gene_input_mp[[1]], isoform = gene_input_mp[[1]])
multi_color <- list(gene = iso_gene_colors)

mp <- modulePreservation(
    multiData = multi_expr, multiColor = multi_color,
    dataIsExpr = TRUE, networkType = "signed",
    referenceNetworks = 1,
    testNetworks = 1,
    nPermutations = 1,
    quickCor = 1,
    verbose = 5
)

# Jaccard analysis

isoform_modules <- module_assignments %>%
    filter(module_label != 0) %>%
    filter(Data == "Isoform") %>%
    pull(module_label) %>%
    unique()
gene_modules <- module_assignments %>%
    filter(module_label != 0) %>%
    filter(Data == "Gene") %>%
    pull(module_label) %>%
    unique()

module_jaccard <- expand.grid(
    isoform_module = isoform_modules,
    gene_module = gene_modules
) %>%
    mutate(
        intersection = mapply(
            function(i, g) {
                gl <- pull(filter(module_assignments, Data == "Gene", module_label == g), 'ensembl_gene_id')
                il <- pull(filter(module_assignments, Data == "Isoform", module_label == i), 'ensembl_gene_id')
                length(intersect(gl, il))
            },
            isoform_module, gene_module
        ),
        union = mapply(
            function(i, g) {
                gl <- pull(filter(module_assignments, Data == "Gene", module_label == g), 'ensembl_gene_id')
                il <- pull(filter(module_assignments, Data == "Isoform", module_label == i), 'ensembl_gene_id')
                length(union(gl, il))
            },
            isoform_module, gene_module
        )
    ) %>%
    mutate(jaccard = intersection / union) %>%
    mutate(jaccard_dissimilarity = 1 - jaccard)
