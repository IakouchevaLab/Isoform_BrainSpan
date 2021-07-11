setwd("/Users/pramod/Desktop/Proteomics/16p_ipsc_VP/")
library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(doParallel)
registerDoParallel()
data<-get(load("1month_Organoids_RSEM_Quant.genes.counts.RData"))
meta<-readxl::read_xlsx("1m_metadata.xlsx")
# data<- data[, -c(37:38)]
data<-as.data.frame(data)
data1<- data[,match(colnames(data), meta$Sample)]
data1<- data[,meta$Sample]
cbind(colnames(data1), meta$Sample)

Expressed_genes<-readxl::read_xlsx("16p11.2_list_of_expressed_genes.xlsx", sheet = 1)
data2<- data1[Expressed_genes$`1m_expressed_genes`,]
data2<-na.omit(data2)


# Variance Partition
gExpr <- DGEList(counts=data2)
gExpr <- calcNormFactors(gExpr)
design <- model.matrix( ~ Genotype + Run +Clone, meta)
vobjGenes <- voom(gExpr, design )
geneExpres <-vobjGenes$E

form <- ~ (1|Genotype) + (1|Individual) + (1|Run) + (1|Clone) 
varPart <- fitExtractVarPartModel( vobjGenes, form, meta)
vp <- sortCols(varPart )
plotvp<-plotVarPart(vp)
ggsave("Vp_1m_Clone.pdf",plotvp)
plotPercentBars(vp[1:10,])
vp<-as.data.frame( vp)
vp$Ensembl_ID<-rownames(vp)
writexl::write_xlsx(vp,"1m_Org_VariancePartition_variability_Clone.xlsx", col_names = T)

# Top most variable genes
sortedLab <- (varPart[order(varPart$Clone, decreasing=TRUE),])
i <- sortedLab$Clone > 0.5
topLab <- sortedLab[i, ]
topLab$Ensembl_ID<-rownames(topLab)
writexl::write_xlsx(topLab, "1m_Genes_mostVariable_byClone_0.5.xlsx", col_names = T)

# compare<-readxl::read_xlsx("Genes_mostVariable_byLab_0.5_patricia.xlsx", sheet = 1)
# head(topLab)
# dim(compare)
# dim(topLab)
# common<-intersect(compare$Ensembl_ID, rownames(topLab))
# writexl::write_xlsx(as.data.frame(common), "common_1m.xlsx", col_names = TRUE)
