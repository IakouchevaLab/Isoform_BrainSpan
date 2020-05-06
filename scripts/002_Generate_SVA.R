library(tidyverse)
library(edgeR)
library(sva)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

# i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
dir.create("data/isoforms/sva_results")
dir.create("data/genes/sva_results")
metadata <- read_tsv("data/metadata.tsv")

iso_cts <- readRDS("data/iso_cts_filter.rds")
iso_cpm <- cpm(calcNormFactors(DGEList(iso_cts), method = "TMM"))
m1 <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
m0 <- model.matrix(
    ~ Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
foreach(i = 1:20) %dopar% {
    fn <- paste0("data/isoforms/sva_results/", i, ".rds")
    sv <- svaseq(iso_cpm, m1, m0, n.sv = i)
    rownames(sv$sv) <- metadata$Sample
    colnames(sv$sv) <- paste0("SV", seq(1, ncol(sv$sv)))
    saveRDS(sv, fn)
}


gn_cts <- readRDS("data/gene_cts_filter.rds")
gn_cpm <- cpm(calcNormFactors(DGEList(gn_cts), method = "TMM"))
m1 <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
m0 <- model.matrix(
    ~ Regioncode + Sex + Ethnicity + Site,
    data = metadata
)
foreach(i = 1:20) %dopar% {
    fn <- paste0("data/genes/sva_results/", i, ".rds")
    sv <- svaseq(gn_cpm, m1, m0, n.sv = i)
    rownames(sv$sv) <- metadata$Sample
    colnames(sv$sv) <- paste0("SV", seq(1, ncol(sv$sv)))
    saveRDS(sv, fn)
}
