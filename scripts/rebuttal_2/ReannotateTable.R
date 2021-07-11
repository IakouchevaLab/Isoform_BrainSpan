library(tidyverse)
library(biomaRt)
mart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)
tbl <- readxl::read_xlsx("data/Supplementary Table 1_NEW_Validation_JU.xlsx", sheet=2)
anno <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "external_transcript_name"),
    filters = c("ensembl_transcript_id"),
    mart = mart,
    values = tbl$ensembl_transcript_id
)
tbl_new <- left_join(tbl, anno, by = c("ensembl_transcript_id" = "ensembl_transcript_id"))
write_csv(tbl_new, "data/Supplementary Table 1_NEW_Validation_JU_KCReannotated.csv")
