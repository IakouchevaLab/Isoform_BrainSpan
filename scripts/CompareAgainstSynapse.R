library(tidyverse)

load("data/source/BrainSpan.RSEM_Quant.isoform.tpm.RData")
iso_tpm <- tpm
load("data/source/BrainSpan.RSEM_Quant.isoform.counts.RData")
iso_cts <- counts
load("data/source/RSEM_Quant.genes.counts.RData")
gene_cts <- counts[, grepl("BrainSpan", colnames(counts))]
load("data/source/RSEM_Quant.genes.tpm.RData")
gene_tpm <- tpm[, grepl("BrainSpan", colnames(counts))]
metadata <- read_tsv("data/source/brainSpan.phenotype.meta.final.tsv") %>%
    mutate(Sample = paste("BrainSpan", Braincode, Regioncode, sep = "_"))

HSB103_AMY_TPM <- read_tsv("Synapse/HSB103_AMY.RSEM_Quant.isoforms.results") %>%
    dplyr::select(
        transcript_id, gene_id, TPM
    ) %>%
    separate(transcript_id, into = c("ensembl_transcript_id", "version"), sep = "\\.")
all(iso_tpm[, "BrainSpan_HSB103_AMY"][pull(HSB103_AMY_TPM, ensembl_transcript_id)] == pull(HSB103_AMY_TPM, TPM))

HSB113_A1C_TPM <- read_tsv("Synapse/HSB113_A1C.RSEM_Quant.isoforms.results") %>%
    dplyr::select(
        transcript_id, gene_id, TPM
    ) %>%
    separate(transcript_id, into = c("ensembl_transcript_id", "version"), sep = "\\.")
all(iso_tpm[, "BrainSpan_HSB113_A1C"][pull(HSB113_A1C_TPM, ensembl_transcript_id)] == pull(HSB113_A1C_TPM, TPM))
