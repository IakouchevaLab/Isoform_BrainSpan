library(tidyverse)
library(biomaRt)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt",
    header = TRUE, sep = ",", row.names = 1,
    stringsAsFactors = FALSE
) %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "gene_biotype"),
            filters = c("ensembl_gene_id"),
            values = unique(.$ensembl_gene_id),
            mart = useMart(
                biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl",
                host = "GRCh37.ensembl.org"
            )
        ),
        by = "ensembl_gene_id"
    )

deg <- lapply(
    readxl::excel_sheets("data/SupplementaryTables/Supplementary Table 2.xlsx")[-1],
    function(sheet) {
        readxl::read_xlsx("data/SupplementaryTables/Supplementary Table 2.xlsx", sheet = sheet)
    }
) %>%
    bind_rows() %>%
    left_join(
        annotations[, c("ensembl_gene_id", "gene_biotype")],
        by = c("Ensembl Gene ID" = "ensembl_gene_id")
    )
dei <- lapply(
    readxl::excel_sheets("data/SupplementaryTables/Supplementary Table 3.xlsx")[-1],
    function(sheet) {
        readxl::read_xlsx("data/SupplementaryTables/Supplementary Table 3.xlsx", sheet = sheet)
    }
) %>%
    bind_rows() %>%
    left_join(
        annotations[, c("ensembl_transcript_id", "transcript_biotype")],
        by = c("Ensembl Transcript ID" = "ensembl_transcript_id")
    )
gmod <- lapply(
    readxl::excel_sheets("data/SupplementaryTables/Supplementary Table 7.xlsx")[-1],
    function(sheet) {
        readxl::read_xlsx("data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = sheet)
    }
) %>%
    bind_rows() %>%
    left_join(
        annotations[, c("ensembl_gene_id", "gene_biotype")],
        by = c("Ensembl Gene ID" = "ensembl_gene_id")
    )
imod <- lapply(
    readxl::excel_sheets("data/SupplementaryTables/Supplementary Table 8.xlsx")[-1],
    function(sheet) {
        readxl::read_xlsx("data/SupplementaryTables/Supplementary Table 8.xlsx", sheet = sheet)
    }
) %>%
    bind_rows() %>%
    left_join(
        annotations[, c("ensembl_transcript_id", "transcript_biotype")],
        by = c("Ensembl Transcript ID" = "ensembl_transcript_id")
    )

deg_pc <- deg %>% {
    df <- .
    tibble(
        Contrast = unique(.$Contrast)
    ) %>%
        mutate(
            NumDE = sapply(Contrast, function(ctr) {
                filter(df, Contrast == ctr) %>%
                    filter(abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
                    nrow()
            }),
            NumDE_proteincoding = sapply(Contrast, function(ctr) {
                filter(df, Contrast == ctr) %>%
                    filter(abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
                    filter(gene_biotype == "protein_coding") %>%
                    nrow()
            })
        )
} %>%
    mutate(PC_proportion = NumDE_proteincoding / NumDE)

dei_pc <- dei %>% {
    df <- .
    tibble(
        Contrast = unique(.$Contrast)
    ) %>%
        mutate(
            NumDE = sapply(Contrast, function(ctr) {
                filter(df, Contrast == ctr) %>%
                    filter(abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
                    nrow()
            }),
            NumDE_proteincoding = sapply(Contrast, function(ctr) {
                filter(df, Contrast == ctr) %>%
                    filter(abs(logFC) >= log2(1.5) & FDR <= 0.05) %>%
                    filter(transcript_biotype == "protein_coding") %>%
                    nrow()
            })
        )
}

gmod_pc <- gmod %>% {
    df <- .
    tibble(
        module_label = .$`Module Label`,
        module_colour = .$`Module Colour`,
    ) %>%
        distinct() %>%
        mutate(
            NumTotal = sapply(module_label, function(ml) {
                nrow(filter(df, `Module Label` == ml))
            }),
            NumPC = sapply(module_label, function(ml) {
                nrow(filter(df, `Module Label` == ml & gene_biotype == "protein_coding"))
            })
        )
} %>%
    mutate(PC_proportion = NumPC / NumTotal)

imod_pc <- imod %>% {
    df <- .
    tibble(
        module_label = .$`Module Label`,
        module_colour = .$`Module Colour`,
    ) %>%
        distinct() %>%
        mutate(
            NumTotal = sapply(module_label, function(ml) {
                nrow(filter(df, `Module Label` == ml))
            }),
            NumPC = sapply(module_label, function(ml) {
                nrow(filter(df, `Module Label` == ml & transcript_biotype == "protein_coding"))
            })
        )
} %>%
    mutate(PC_proportion = NumPC / NumTotal)

writexl::write_xlsx(
    list(
        DEG_PC = deg_pc,
        DEI_PC = dei_pc,
        GENE_MOD_PC = gmod_pc,
        ISOFORM_MOD_PC = imod_pc
    ),
    "data/results_protein_coding_proportions.xlsx"
)
