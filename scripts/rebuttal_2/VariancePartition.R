library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(doParallel)
registerDoParallel(cores = detectCores() - 1)

metadata <- read_tsv("data/metadata.tsv")
cell_types_damon_ar <- read_csv("data/estimated_props_Damon_afterRegression.csv") %>%
    rename(Sample = X1)
cell_types_damon_br <- read_csv("data/estimated_props_Damon_beforeRegression.csv") %>%
    rename(Sample = X1)
gene_cts <- readRDS("data/gene_cts_filter.rds")
sva_results <- readRDS(paste0("data/genes/sva_results/", 16, ".rds"))


pdf("data/figures/rebuttal_2/CellTypeDamon_VariancePartition.pdf", paper = 'a4r')
# Raw counts data
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(mm, cell_types_damon_br[which(colnames(cell_types_damon_br) != "Sample")])
meta <- left_join(metadata, cell_types_damon_br, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + ExDp1 + ExDp2 + ExM + `ExM-U` + ExN + InCGE + InMGE + IP + OPC + oRG + Per + PgG2M + PgS + vRG
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeDamon_beforeRegression_CTOnly_Period_Voom_VariancePartition")
plotvp

# Regressed counts data
mm <- model.matrix(
    ~ Period + Regioncode + Sex + Ethnicity + Site,
    data = metadata[match(colnames(gene_cts), metadata$Sample), ] %>%
        mutate(Period = as.factor(Period))
)
rownames(mm) <- metadata$Sample
mm <- cbind(mm, sva_results$sv)
mm_ct <- cbind(mm, cell_types_damon_ar[which(colnames(cell_types_damon_ar) != "Sample")])
meta <- left_join(metadata, cell_types_damon_ar, by = "Sample")
ge <- DGEList(counts=gene_cts) %>%
    calcNormFactors()
vobj <- voom(ge, design = mm_ct)
gene_expression <- vobj$E
form <- ~ Period + ExDp1 + ExDp2 + ExM + `ExM-U` + ExN + InCGE + InMGE + IP + OPC + oRG + Per + PgG2M + PgS + vRG
varPart <- fitExtractVarPartModel(gene_expression, form, meta)
vp <- sortCols(varPart)
plotvp <- plotVarPart(vp, main="CellTypeDamon_afterRegression_CTOnly_Period_Voom_VariancePartition", paper = 'a4r')
plotvp
dev.off()

