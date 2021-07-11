library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(doParallel)
registerDoParallel(cores = detectCores() - 1)

metadata <- read_tsv("data/metadata.tsv")
cell_types_zhong_raw <- read_csv("data/estimated_props_Zhong_originalData.csv") %>%
    rename(Sample = X1)
cell_types_zhong_regress <- read_csv("data/estimated_props_Zhong_regressedData.csv") %>%
    rename(Sample = X1)
# cell_types_lake <- read_csv("data/estimated_props_Bisque_lake_level1.csv") %>%
    # rename(Sample = X1)
gene_cts <- readRDS("data/gene_cts_filter.rds")
sva_results <- readRDS(paste0("data/isoforms/sva_results/", 13, ".rds"))


pdf("data/figures/CellTypeZhong_VariancePartition.pdf", width = 8, height = 6)
# Raw counts data
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(mm, cell_types_zhong_raw[which(colnames(cell_types_zhong_raw) != "Sample")])
meta <- left_join(metadata, cell_types_zhong_raw, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + Astrocytes + `GABAergic neurons` + Microglia + Neurons + OPC
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeZhongRaw_Period_Voom_VariancePartition")
plotvp
# Raw counts data - full model
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(mm, cell_types_zhong_raw[which(colnames(cell_types_zhong_raw) != "Sample")])
meta <- left_join(metadata, cell_types_zhong_raw, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + (1|Regioncode) + (1|Sex) + (1|Ethnicity) + (1|Site) + Astrocytes + `GABAergic neurons` + Microglia + Neurons + OPC
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeZhongRaw_FullModel_Voom_VariancePartition")
plotvp

# Regressed counts data
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(
    mm, 
    cell_types_zhong_regress[which(colnames(cell_types_zhong_regress) != "Sample")]
)
meta <- left_join(metadata, cell_types_zhong_regress, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + Astrocytes + `GABAergic neurons` + Microglia + Neurons + OPC
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeZhongRegress_Period_Voom_VariancePartition")
plotvp

# Regressed counts data - Full model
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(
    mm, 
    cell_types_zhong_regress[which(colnames(cell_types_zhong_regress) != "Sample")]
)
meta <- left_join(metadata, cell_types_zhong_regress, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + (1|Regioncode) + (1|Sex) + (1|Ethnicity) + (1|Site) + Astrocytes + `GABAergic neurons` + Microglia + Neurons + OPC
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeZhongRegress_Period_Voom_VariancePartition")
plotvp

dev.off()

# sv <- colnames(dplyr::select(mm, starts_with("SV")))
# cell_types <- c("Ast", "End", "Ex", "In", "Mic", "Oli", "OPC", "Per")
# expand.grid(sv, cell_types) %>%
#     rename(SV = Var1, CellType = Var2) %>%
#     mutate(Correlation = apply(., 1, function(r) {
#         print(r[["SV"]])
#         cor(sva_results$sv[, r[["SV"]]], cell_types_lake[, r[["CellType"]]])
#     }))
# 
# # Linear combinations
# sv_ct <- cbind(sva_results$sv, cell_types_lake[which(colnames(cell_types_lake) != "Sample")])
# caret::findLinearCombos(mm[, 2:ncol(mm)])
