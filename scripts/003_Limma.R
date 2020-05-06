library(tidyverse)
library(edgeR)
library(limma)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

dir.create("data/isoforms/limma_intermediates")
dir.create("data/genes/limma_intermediates")

metadata <- read_tsv("data/metadata.tsv")

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# iso_cts <- readRDS("data/iso_cts_filter.rds")
# sva_results <- readRDS(paste0("data/isoforms/sva_results/", i, ".rds"))
# mm <- model.matrix(
#     ~ Period + Regioncode + Sex + Ethnicity + Site,
#     data = metadata[match(colnames(iso_cts), metadata$Sample), ] %>%
#         mutate(Period = as.factor(Period))
# )
# mm_p2p <- model.matrix(
#     ~ Prenatal + Period + Regioncode + Sex + Ethnicity + Site,
#     data = metadata[match(colnames(iso_cts), metadata$Sample), ] %>%
#         mutate(Period = as.factor(Period))
# )
# rownames(mm) <- metadata$Sample
# rownames(mm_p2p) <- metadata$Sample
# # all(rownames(mm) == rownames(sva_results$sv))
# # colnames(sva_results$sv) <- paste0("SV", seq(1, ncol(sva_results$sv)))
# mm <- cbind(mm, sva_results$sv)
# mm_p2p <- cbind(mm_p2p, sva_results$sv)
# block <- metadata$Braincode[match(colnames(iso_cts), metadata$Sample)]
# iso_norm <- calcNormFactors(DGEList(iso_cts), method = "TMM")
# # ADJ
# iso_voom <- voom(iso_norm, design = mm)
# iso_dc <- duplicateCorrelation(
#     iso_voom,
#     design = mm,
#     block = block
# )
# iso_voom <- voom(
#     iso_norm, design = mm, block = block, correlation = iso_dc$consensus
# )
# iso_dc <- duplicateCorrelation(
#     iso_voom,
#     design = mm,
#     block = block
# )
# lm_fit <- lmFit(
#     iso_voom, design = mm, block = block, correlation = iso_dc$consensus
# )
# saveRDS(
#     lm_fit,
#     paste0("data/isoforms/limma_intermediates/lm_fit_full_", i, ".rds")
# )
# cont_mat <- makeContrasts(
#     P02P03 = Period3,
#     P03P04 = Period4 - Period3,
#     P04P05 = Period5 - Period4,
#     P05P06 = Period6 - Period5,
#     P06P07 = Period7 - Period6,
#     P07P08 = Period8 - Period7,
#     P08P09 = Period9 - Period8,
#     P09P10 = Period10 - Period9,
#     P10P11 = Period11 - Period10,
#     P11P12 = Period12 - Period11,
#     P12P13 = Period13 - Period12,
#     levels = colnames(mm)
# )
# ctr_fit <- contrasts.fit(lm_fit, cont_mat)
# ebayes <- eBayes(ctr_fit)
# tt_full <- lapply(
#     colnames(cont_mat),
#     function(coef) {
#         topTable(ebayes, coef = coef, number = Inf) %>%
#             mutate(
#                 ensembl_transcript_id = rownames(.),
#                 Contrast = coef
#             )
#     }
# ) %>% setNames(., colnames(cont_mat)) %>%
#     bind_rows()
# # P2P
# iso_voom_p2p <- voom(iso_norm, design = mm_p2p)
# iso_dc_p2p <- duplicateCorrelation(
#     iso_voom_p2p,
#     design = mm_p2p,
#     block = block
# )
# iso_voom_p2p <- voom(
#     iso_norm, design = mm_p2p,
#     block = block, correlation = iso_dc_p2p$consensus
# )
# iso_dc_p2p <- duplicateCorrelation(
#     iso_voom_p2p,
#     design = mm_p2p,
#     block = block
# )
# lm_fit_p2p <- lmFit(
#     iso_voom_p2p, design = mm_p2p,
#     block = block, correlation = iso_dc_p2p$consensus
# )
# saveRDS(
#     lm_fit_p2p,
#     paste0("data/isoforms/limma_intermediates/lm_fit_p2p_", i, ".rds")
# )
# cont_mat_p2p <- makeContrasts(
#     PrePost = -PrenatalTRUE,
#     levels = colnames(mm_p2p)
# )
# ctr_fit_p2p <- contrasts.fit(lm_fit_p2p, cont_mat_p2p)
# ebayes_p2p <- eBayes(ctr_fit_p2p)
# tt_p2p <- lapply(
#     colnames(cont_mat_p2p),
#     function(coef) {
#         topTable(ebayes_p2p, coef = coef, number = Inf) %>%
#             mutate(
#                 ensembl_transcript_id = rownames(.),
#                 Contrast = coef
#             )
#     }
# ) %>% setNames(., colnames(cont_mat_p2p)) %>%
#     bind_rows()
# tt <- bind_rows(
#     tt_full, tt_p2p
# )
# saveRDS(
#     tt,
#     paste0("data/isoforms/limma_intermediates/tt_SV", i, ".rds")
# )

gene_cts <- readRDS("data/gene_cts_filter.rds")
sva_results <- readRDS(paste0("data/genes/sva_results/", i, ".rds"))
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
mm_p2p <- model.matrix(
    ~ Prenatal + Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
rownames(mm_p2p) <- metadata$Sample
# all(rownames(mm) == rownames(sva_results$sv))
# colnames(sva_results$sv) <- paste0("SV", seq(1, ncol(sva_results$sv)))
mm <- cbind(mm, sva_results$sv)
mm_p2p <- cbind(mm_p2p, sva_results$sv)
block <- metadata$Braincode[match(colnames(gene_cts), metadata$Sample)]
gene_norm <- calcNormFactors(DGEList(gene_cts), method = "TMM")
# ADJ
gene_voom <- voom(gene_norm, design = mm)
gene_dc <- duplicateCorrelation(
    gene_voom, 
    design = mm, 
    block = block
)
gene_voom <- voom(
    gene_norm, design = mm, 
    block = block, correlation = gene_dc$consensus
)
gene_dc <- duplicateCorrelation(
    gene_voom, 
    design = mm, 
    block = block
)
lm_fit <- lmFit(
    gene_voom, design = mm, 
    block = block, correlation = gene_dc$consensus
)
saveRDS(
    lm_fit, 
    paste0("data/genes/limma_intermediates/lm_fit_full_", i, ".rds")
)
cont_mat <- makeContrasts(
    P02P03 = Period3,
    P03P04 = Period4 - Period3,
    P04P05 = Period5 - Period4,
    P05P06 = Period6 - Period5,
    P06P07 = Period7 - Period6,
    P07P08 = Period8 - Period7,
    P08P09 = Period9 - Period8,
    P09P10 = Period10 - Period9,
    P10P11 = Period11 - Period10,
    P11P12 = Period12 - Period11,
    P12P13 = Period13 - Period12,
    levels = colnames(mm)
)
ctr_fit <- contrasts.fit(lm_fit, cont_mat)
ebayes <- eBayes(ctr_fit)
tt_full <- lapply(
    colnames(cont_mat),
    function(coef) {
        topTable(ebayes, coef = coef, number = Inf) %>%
            mutate(
                ensembl_gene_id = rownames(.),
                Contrast = coef
            )
    }
) %>% setNames(., colnames(cont_mat)) %>%
    bind_rows()
# P2P
gene_voom_p2p <- voom(gene_norm, design = mm_p2p)
gene_dc_p2p <- duplicateCorrelation(
    gene_voom_p2p, 
    design = mm_p2p, 
    block = block
)
gene_voom_p2p <- voom(
    gene_norm, design = mm_p2p, 
    block = block, correlation = gene_dc_p2p$consensus
)
gene_dc_p2p <- duplicateCorrelation(
    gene_voom_p2p, 
    design = mm_p2p, 
    block = block
)
lm_fit_p2p <- lmFit(
    gene_voom_p2p, design = mm_p2p, 
    block = block, correlation = gene_dc_p2p$consensus
)
saveRDS(
    lm_fit_p2p, 
    paste0("data/genes/limma_intermediates/lm_fit_p2p_", i, ".rds")
)
cont_mat_p2p <- makeContrasts(
    PrePost = -PrenatalTRUE,
    levels = colnames(mm_p2p)
)
ctr_fit_p2p <- contrasts.fit(lm_fit_p2p, cont_mat_p2p)
ebayes_p2p <- eBayes(ctr_fit_p2p)
tt_p2p <- lapply(
    colnames(cont_mat_p2p),
    function(coef) {
        topTable(ebayes_p2p, coef = coef, number = Inf) %>%
            mutate(
                ensembl_gene_id = rownames(.),
                Contrast = coef
            )
    }
) %>% setNames(., colnames(cont_mat_p2p)) %>%
    bind_rows()
tt <- bind_rows(
    tt_full, tt_p2p
)
saveRDS(
    tt, 
    paste0("data/genes/limma_intermediates/tt_SV", i, ".rds")
)
