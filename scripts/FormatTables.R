library(tidyverse)
library(readxl)
library(writexl)

annotations <- read.csv(
    "data/source/annotation.transcript.ensg75.txt", header = TRUE, row.names = 1
)

# DE Tables
contrasts <- setNames(nm = c(
    "P02P03", "P03P04", "P04P05", "P05P06", "P06P07", "P07P08",
    "P08P09", "P09P10", "P11P12", "P12P13", "PrePost"
))
gde <- readRDS("data/genes/limma_intermediates/tt_SV16.rds") %>%
    left_join(
        annotations %>% 
            dplyr::select(
                ensembl_gene_id, external_gene_id
            ) %>%
            distinct(),
        by = "ensembl_gene_id"
    ) %>%
    rename(
        `Gene Symbol` = external_gene_id, `Ensembl Gene ID` = ensembl_gene_id,
        `P-Value` = P.Value, FDR = adj.P.Val
    ) %>%
    dplyr::select(
        `Gene Symbol`, `Ensembl Gene ID`, logFC, AveExpr, t, B, `P-Value`, FDR,
        Contrast
    )
gde_list <- lapply(contrasts, function(ctr) filter(gde, Contrast == ctr))

write_xlsx(gde_list, "data/SupplementaryTables/DEGenes.xlsx")

ide <- readRDS("data/isoforms/limma_intermediates/tt_SV16.rds") %>%
    left_join(
        annotations %>% 
            dplyr::select(
                ensembl_gene_id, external_gene_id, 
                ensembl_transcript_id, external_transcript_id
            ) %>%
            distinct(),
        by = "ensembl_transcript_id"
    ) %>%
    rename(
        `Gene Symbol` = external_gene_id, 
        `Transcript Symbol` = external_transcript_id,
        `Ensembl Gene ID` = ensembl_gene_id,
        `Ensembl Transcript ID` = ensembl_transcript_id,
        `P-Value` = P.Value, FDR = adj.P.Val
    ) %>%
    dplyr::select(
        `Gene Symbol`, `Transcript Symbol`, 
        `Ensembl Gene ID`, `Ensembl Transcript ID`,
        logFC, AveExpr, t, B, `P-Value`, FDR,
        Contrast
    )
ide_list <- lapply(contrasts, function(ctr) filter(ide, Contrast == ctr))
write_xlsx(ide_list, "data/SupplementaryTables/DEIsoforms.xlsx")

# DE Enrichments
gene_de_go <- readRDS("data/genes/Enrichments/DE_GOST_SPEC.rds") %>%
    lapply(
        ., function(tb) {
            if (is.null(tb)) {
                return(data.frame())
            }
            tb %>%
                filter(term_size <= 1000) %>%
                dplyr::select(
                    source, term_id, term_name, p_value, intersection_size
                ) %>%
                rename(
                    Source = source, 
                    `Term ID` = term_id,
                    `Term Name` = term_name,
                    `FDR` = p_value,
                    `Intersection Size` = intersection_size
                ) %>%
                arrange(FDR)
        }
    )
write_xlsx(gene_de_go, "data/SupplementaryTables/GeneDEGOST.xlsx")

isoform_de_go <- readRDS("data/isoforms/Enrichments/DE_GOST_SPEC.rds") %>%
    lapply(
        ., function(tb) {
            if (is.null(tb)) {
                return(data.frame())
            }
            tb %>%
                filter(term_size <= 1000) %>%
                dplyr::select(
                    source, term_id, term_name, p_value, intersection_size
                ) %>%
                rename(
                    Source = source, 
                    `Term ID` = term_id,
                    `Term Name` = term_name,
                    `FDR` = p_value,
                    `Intersection Size` = intersection_size
                ) %>%
                arrange(FDR)
        }
    )
write_xlsx(isoform_de_go, "data/SupplementaryTables/IsoformDEGOST.xlsx")

# Module assignments
gene_modules <- read_tsv(
    "data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv"
) %>%
    rename(
        `Ensembl Gene ID` = Feature, 
        `Module Label` = module_label, 
        `Module Colour` = module_colour
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(ensembl_gene_id, external_gene_id),
        by = c("Ensembl Gene ID" = "ensembl_gene_id")
    ) %>%
    rename(`Gene Symbol` = external_gene_id) %>%
    gather(key = "kMEKey", value = "kME", starts_with("kME")) %>%
    filter(as.numeric(gsub("kME", "", kMEKey)) == `Module Label`) %>%
    dplyr::select(
        `Gene Symbol`, `Ensembl Gene ID`, `Module Label`, `Module Colour`, kME
    ) %>%
    arrange(`Module Label`)
gene_modules_list <- lapply(
    setNames(nm = sort(unique(gene_modules$`Module Label`))),
    function(ml) {
        gene_modules %>%
            filter(`Module Label` == ml) %>%
            arrange(desc(kME)) %>%
            distinct()
    }
)
write_xlsx(
    gene_modules_list, "data/SupplementaryTables/GeneModuleAssignment.xlsx"
)

isoform_modules <- read_tsv(
    "data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv"
) %>%
    rename(
        `Ensembl Transcript ID` = Feature, 
        `Module Label` = module_label, 
        `Module Colour` = module_colour
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(
                ensembl_gene_id, external_gene_id,
                ensembl_transcript_id, external_transcript_id
            ),
        by = c(
            "Ensembl Transcript ID" = "ensembl_transcript_id"
        )
    ) %>%
    gather(key = "kMEKey", value = "kME", starts_with("kME")) %>%
    filter(as.numeric(gsub("kME", "", kMEKey)) == `Module Label`) %>%
    rename(
        `Gene Symbol` = external_gene_id,
        `Transcript Symbol` = external_transcript_id,
        `Ensembl Gene ID` = ensembl_gene_id
    ) %>%
    dplyr::select(
        `Gene Symbol`, `Transcript Symbol`, 
        `Ensembl Gene ID`, `Ensembl Transcript ID`,
        `Module Label`, `Module Colour`, kME
    ) %>%
    arrange(`Module Label`)
isoform_modules_list <- lapply(
    setNames(nm = sort(unique(isoform_modules$`Module Label`))),
    function(ml) {
        isoform_modules %>%
            filter(`Module Label` == ml) %>%
            arrange(desc(kME))
    }
)
write_xlsx(
    isoform_modules_list, 
    "data/SupplementaryTables/IsoformModuleAssignment.xlsx"
)

# Module enrichments
gene_module_go <- readRDS("data/WGCNA/GeneNetwork_DS2_MM20_GOST.rds") %>%
    lapply(
        ., function(tb) {
            if (is.null(tb)) {
                return(data.frame())
            }
            tb %>%
                filter(term_size <= 1000) %>%
                dplyr::select(
                    source, term_id, term_name, p_value, intersection_size
                ) %>%
                rename(
                    Source = source, 
                    `Term ID` = term_id,
                    `Term Name` = term_name,
                    `FDR` = p_value,
                    `Intersection Size` = intersection_size
                ) %>%
                arrange(FDR)
        }
    )
write_xlsx(gene_module_go, "data/SupplementaryTables/GeneModuleGOST.xlsx")
isoform_module_go <- readRDS("data/WGCNA/IsoformNetwork_DS2_MM20_GOST.rds") %>%
    lapply(
        ., function(tb) {
            if (is.null(tb)) {
                return(data.frame())
            }
            tb %>%
                filter(term_size <= 1000) %>%
                dplyr::select(
                    source, term_id, term_name, p_value, intersection_size
                ) %>%
                rename(
                    Source = source, 
                    `Term ID` = term_id,
                    `Term Name` = term_name,
                    `FDR` = p_value,
                    `Intersection Size` = intersection_size
                ) %>%
                arrange(FDR)
        }
    )
write_xlsx(isoform_module_go, "data/SupplementaryTables/IsoformModuleGOST.xlsx")
