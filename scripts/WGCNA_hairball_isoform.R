library(tidyverse)
library(igraph)
library(ggnetwork)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", 
    header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE
)

gexpr <- readRDS("data/RegressGeneCounts.rds")
iexpr <- readRDS("data/RegressIsoformCounts.rds")

gene_modules <- lapply(
    setNames(
        nm = readxl::excel_sheets(
            "data/SupplementaryTables/Supplementary Table 7.xlsx"
        )[-1]
    ),
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 7.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()
isoform_modules <- lapply(
    setNames(
        nm = readxl::excel_sheets(
            "data/SupplementaryTables/Supplementary Table 8.xlsx"
        )[-1]
    ),
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 8.xlsx",
            sheet = sheet
        )
    }
) %>%
    bind_rows()

variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 6.xlsx",
    sheet = 2
)
control_targets <- c(
    filter(variants, `Affected status` == 1)$`Ensembl Gene ID`,
    filter(variants, `Affected status` == 1)$`Ensembl Transcript ID`
)
asd_targets <- c(
    filter(variants, `Affected status` == 2)$`Ensembl Gene ID`,
    filter(variants, `Affected status` == 2)$`Ensembl Transcript ID`
)
asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]
asd_genes_ensembl <- as.character(annotations$ensembl_gene_id[match(
    asd_genes, annotations$external_gene_id
)])

vertices <- bind_rows(
    gene_modules %>%
        mutate(name = `Ensembl Gene ID`) %>%
        mutate(Data = "Gene"),
    isoform_modules %>%
        mutate(name = `Ensembl Transcript ID`) %>%
        mutate(Data = "Isoform")
) %>%
    dplyr::select(name, everything()) %>%
    mutate(ASD = `Ensembl Gene ID` %in% asd_genes_ensembl) %>%
    mutate(
        Target = case_when(
            name %in% control_targets & name %in% asd_targets ~ "Both",
            ! name %in% control_targets & name %in% asd_targets ~ "ASD",
            name %in% control_targets & ! name %in% asd_targets ~ "Control",
            ! name %in% control_targets & ! name %in% asd_targets ~ "Neither"
        )
    )

yang_ppi <- readLines(
    "data/source/CuratedLists/Yang-Vidal_Cell_2016.psi"
) %>%
    str_split(" ", n = 5) %>%
    unlist() %>%
    matrix(ncol = 5, byrow = TRUE) %>%
    as.data.frame() %>%
    setNames(nm = c("Uniprot1", "Uniprot2", "Ensembl1", "Ensembl2", "Info")) %>%
    group_by(Uniprot1, Uniprot2, Ensembl1, Ensembl2) %>%
    summarise(Info = toString(Info)) %>%
    ungroup() %>%
    mutate(
        Ensembl1 = str_replace_all(Ensembl1, "ensembl:", ""),
        Ensembl2 = str_replace_all(Ensembl2, "ensembl:", "")
    ) %>%
    mutate(
        Ensembl1 = str_replace_all(Ensembl1, "\\.[0-9]+", ""),
        Ensembl2 = str_replace_all(Ensembl2, "\\.[0-9]+", "")
    ) %>%
    mutate(
        ENSG1 = as.character(str_match(Ensembl1, "ENSG[0-9]+")),
        ENSG2 = as.character(str_match(Ensembl2, "ENSG[0-9]+")),
        ENST1 = as.character(str_match(Ensembl1, "ENST[0-9]+")),
        ENST2 = as.character(str_match(Ensembl2, "ENST[0-9]+"))
    ) %>%
    dplyr::select(starts_with("ENS", ignore.case = FALSE)) %>%
    mutate(
        G1inData = ENSG1 %in% rownames(gexpr), 
        G2inData = ENSG2 %in% rownames(gexpr), 
        I1inData = ENST1 %in% rownames(iexpr),
        I2inData = ENST2 %in% rownames(iexpr)
    ) %>%
    mutate(
        GenePCC = mapply(
            function(g1indata, g2indata, g1, g2) {
                if (!g1indata | !g2indata) {
                    return(0)
                } else {
                    cor(gexpr[g1, ], gexpr[g2, ])
                }
            }, G1inData, G2inData, ENSG1, ENSG2
        ),
        IsoformPCC = mapply(
            function(i1indata, i2indata, i1, i2) {
                if (!i1indata | !i2indata) {
                    return(0)
                } else {
                    cor(iexpr[i1, ], iexpr[i2, ])
                }
            }, I1inData, I2inData, ENST1, ENST2
        )
    )
isoform_ppi <- yang_ppi

# Filter for significant modules and get module networks

gM1_edges <- isoform_ppi %>%
    dplyr::select(contains("ENSG"), GenePCC) %>%
    rename(F1 = ENSG1, F2 = ENSG2, weight = GenePCC) %>%
    filter(
        F1 %in% pull(filter(vertices, Data == "Gene" & `Module Label` == 1), name)
    ) %>%
    filter(
        F2 %in% pull(filter(vertices, Data == "Gene" & `Module Label` == 1), name)
    ) %>%
    filter(
        (F1 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name)) | 
            (F2 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name))
    )
iM1_edges <- isoform_ppi %>%
    dplyr::select(contains("ENST"), IsoformPCC) %>%
    rename(F1 = ENST1, F2 = ENST2, weight = IsoformPCC) %>%
    filter(
        F1 %in% pull(filter(vertices, Data == "Isoform" & `Module Label` == 1), name)
    ) %>%
    filter(
        F2 %in% pull(filter(vertices, Data == "Isoform" & `Module Label` == 1), name)
    ) %>%
    filter(
        (F1 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name)) | 
            (F2 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name))
    )
iM30_edges <- isoform_ppi %>%
    dplyr::select(contains("ENST"), IsoformPCC) %>%
    rename(F1 = ENST1, F2 = ENST2, weight = IsoformPCC) %>%
    filter(
        F1 %in% pull(filter(vertices, Data == "Isoform" & `Module Label` == 30), name)
    ) %>%
    filter(
        F2 %in% pull(filter(vertices, Data == "Isoform" & `Module Label` == 30), name)
    ) %>%
    filter(
        (F1 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name)) | 
            (F2 %in% pull(filter(vertices, Target %in% c("ASD", "Both")), name))
    )
writexl::write_xlsx(
    list(
        gM1 = gM1_edges, iM1 = iM1_edges, iM30 = iM30_edges, Vertices = vertices
    ),
    "data/ppi/IsoformPPIAnnotated.xlsx"
)

build_graph <- function(in_edges, vertices, module_colour) {
    this_graph <- graph_from_data_frame(
        setNames(in_edges, c("F1", "F2", "weight")), 
        vertices = vertices %>%
            filter(name %in% unlist(c(in_edges[, 1], in_edges[, 2])))
    )
    vertex_layout <- layout_with_fr(this_graph) %>%
        as.data.frame() %>%
        setNames(nm = c("x", "y")) %>%
        mutate(name = names(V(this_graph))) %>%
        left_join(vertices, by = "name") %>%
        mutate(Target = Target %in% c("ASD", "Both"))
    edge_layout <- setNames(in_edges, c("F1", "F2", "weight")) %>%
        left_join(
            vertex_layout %>%
                dplyr::select(name, x, y),
            by = c("F1" = "name")
        ) %>%
        rename(x1 = x, y1 = y) %>%
        left_join(
            vertex_layout %>%
                dplyr::select(name, x, y),
            by = c("F2" = "name")
        ) %>%
        rename(x2 = x, y2 = y)
    ggplot() +
        geom_segment(
            data = edge_layout,
            mapping = aes(
                x = x1, y = y1, xend = x2, yend = y2, size = weight
            ),
            colour = "lightgrey"
        ) +
        geom_point(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, shape = ASD, colour = Target
            ),
            size = 5, stroke = 1.5, fill = module_colour
        ) +
        ggrepel::geom_text_repel(
            data = vertex_layout,
            mapping = aes(
                x = x, y = y, 
                label = ifelse(
                    vertex_layout$Data == "Gene",
                    `Gene Symbol`, `Transcript Symbol`
                )
            ),
            size = 5
        ) +
        scale_shape_manual(
            values = c("TRUE" = 24, "FALSE" = 21),
            limits = c("TRUE", "FALSE"),
            breaks = c("TRUE", "FALSE")
        ) +
        scale_colour_manual(
            values = c("TRUE" = "red", "FALSE" = "black"),
            limits = c("TRUE", "FALSE")
        ) +
        scale_size_continuous(
            limits = c(
                floor(min(edge_layout$weight) * 10) / 10,
                ceiling(max(edge_layout$weight) * 10) / 10
            ), 
            range = c(0.5, 3)
        ) +
        theme_blank()
}





pdf("data/ppi/IsoformPPI.pdf", width = 16, height = 12)
tryCatch({build_graph(gM1_edges, vertices, "turquoise")}, error = function(e) { return(NULL) })
tryCatch({build_graph(iM1_edges, vertices, "turquoise")}, error = function(e) { return(NULL) })
tryCatch({build_graph(iM30_edges, vertices, "turquoise")}, error = function(e) { return(NULL) })
dev.off()
