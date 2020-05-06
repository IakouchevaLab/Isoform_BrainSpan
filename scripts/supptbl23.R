library(biomaRt)
mart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)



deg_raw <- readRDS("data/genes/limma_intermediates/tt_SV16.rds") %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "external_gene_name"),
            filters = "ensembl_gene_id",
            values = unique(.$ensembl_gene_id),
            mart = mart
        ),
        by = "ensembl_gene_id"
    ) %>%
    setNames(
        nm = c("logFC", "AveExpr", "t", "P Value", "FDR", "B", "Ensembl Gene ID", "Contrast", "Gene Symbol")
    ) %>%
    dplyr::select(
        c("Gene Symbol", "Ensembl Gene ID", "logFC", "AveExpr", "t", "P Value", "FDR", "B", "Contrast")
    )
deg_list <- lapply(
    setNames(nm = sort(unique(deg_raw$Contrast))),
    function(ctr) {
        filter(deg_raw, Contrast == ctr) %>%
            arrange(FDR)
    }
) %>%
    writexl::write_xlsx(
        "data/SupplementaryTables/Supplementary Table 2.xlsx"
    )
dei_raw <- readRDS("data/isoforms/limma_intermediates/tt_SV16.rds") %>%
    left_join(
        getBM(
            attributes = c("ensembl_transcript_id", "external_gene_name", "external_transcript_name", "ensembl_gene_id"),
            filters = "ensembl_transcript_id",
            values = unique(.$ensembl_transcript_id),
            mart = mart
        ),
        by = "ensembl_transcript_id"
    ) %>%
    setNames(
        nm = c("logFC", "AveExpr", "t", "P Value", "FDR", "B", "Ensembl Transcript ID", "Contrast", "Gene Symbol", "Transcript Symbol", "Ensembl Gene ID")
    ) %>%
    dplyr::select(
        c("Gene Symbol", "Transcript Symbol", "Ensembl Gene ID", "Ensembl Transcript ID", "logFC", "AveExpr", "t", "P Value", "FDR", "B", "Contrast")
    )
dei_list <- lapply(
    setNames(nm = sort(unique(dei_raw$Contrast))),
    function(ctr) {
        filter(dei_raw, Contrast == ctr) %>%
            arrange(FDR)
    }
) %>%
    writexl::write_xlsx(
        "data/SupplementaryTables/Supplementary Table 3.xlsx"
    )


deg <- lapply(
    readxl::excel_sheets(
        "data/SupplementaryTables/Supplementary Table 2.xlsx"
    )[-1],
    function(sheet) {
        print(sheet)
        readxl::read_xlsx(
            path = "data/SupplementaryTables/Supplementary Table 2.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()
dei <- lapply(
    readxl::excel_sheets(
        "data/SupplementaryTables/Supplementary Table 3.xlsx"
    )[-1],
    function(sheet) {
        print(sheet)
        readxl::read_xlsx(
            path = "data/SupplementaryTables/Supplementary Table 3.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()


de_pct <- expand.grid(
    comparison = c(
        "P02P03", "P03P04", "P04P05", "P05P06", "P06P07", "P07P08", 
        "P08P09", "P09P10", "P10P11", "P11P12", "P12P13",
        "PrePost", "AmongPre", "AmongPost"
    ),
    data = c("gene", "isoform", "isoform as gene"),
    stringsAsFactors = FALSE
) %>%
    as_tibble() %>%
    bind_cols(
        apply(., 1, function(x) {
            de <- if (x[["data"]] == "gene") {
                mutate(deg, feature = `Ensembl Gene ID`)
            } else if (x[["data"]] == "isoform") {
                mutate(dei, feature = `Ensembl Transcript ID`)
            } else {
                mutate(dei, feature = `Ensembl Gene ID`)
            }
            num_de <- de %>%
                filter(
                    if (x[["comparison"]] == "AmongPre") {
                        Contrast %in% c("P02P03", "P03P04", "P04P05", "P05P06", "P06P07")
                    } else if (x[["comparison"]] == "AmongPost") {
                        Contrast %in% c("P08P09", "P09P10", "P10P11", "P10P11", "P11P12", "P12P13")
                    } else {
                        Contrast == x[["comparison"]]
                    }
                ) %>%
                filter(abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
                distinct(feature) %>%
                nrow()
            tibble(
                num_de = num_de,
                total_analyzed = nrow(distinct(de, feature)),
            )
        }) %>%
            bind_rows()
    ) %>%
    mutate(pct_de = num_de * 100 / total_analyzed)

write_csv(de_pct, "data/de_pct.csv")

