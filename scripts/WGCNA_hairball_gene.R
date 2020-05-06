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
gM1 <- filter(vertices, Data == "Gene" & `Module Label` == 1)
iM1 <- filter(vertices, Data == "Isoform" & `Module Label` == 1)
iM30 <- filter(vertices, Data == "Isoform" & `Module Label` == 30)

gene_ppi <- read.table(
    "data/source/CuratedLists/allVidal.lit17.mvp.hprdInvivo.hintBinary.txt",
    header = TRUE, stringsAsFactors = FALSE
)[, 1:2] %>%
    setNames(nm = c("ENSG1", "ENSG2"))

gM1_edges <- gene_ppi %>%
    filter(ENSG1 %in% gM1$name & ENSG2 %in% gM1$name) %>%
    distinct() %>%
    filter(ENSG1 != ENSG2) %>%
    mutate(
        weight = mapply(
            function(f1, f2) {
                cor(gexpr[f1, ], gexpr[f2, ])
            }, 
            ENSG1, ENSG2
        )
    )
iM1_edges <- gene_ppi %>%
    left_join(
        annotations %>%
            dplyr::select(ensembl_gene_id, ensembl_transcript_id),
        by = c("ENSG1" = "ensembl_gene_id")
    ) %>%
    rename(ENST1 = ensembl_transcript_id) %>%
    left_join(
        annotations %>%
            dplyr::select(ensembl_gene_id, ensembl_transcript_id),
        by = c("ENSG2" = "ensembl_gene_id")
    ) %>%
    rename(ENST2 = ensembl_transcript_id) %>%
    dplyr::select(ENST1, ENST2) %>%
    distinct() %>%
    filter(ENST1 %in% iM1$name & ENST2 %in% iM1$name) %>%
    filter(ENST1 != ENST2) %>%
    mutate(
        weight = mapply(
            function(f1, f2) {
                cor(iexpr[f1, ], iexpr[f2, ])
            }, 
            ENST1, ENST2
        )
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(
                ensembl_gene_id, ensembl_transcript_id
            ),
        by = c("ENST1" = "ensembl_transcript_id")
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(
                ensembl_gene_id, ensembl_transcript_id
            ),
        by = c("ENST2" = "ensembl_transcript_id")
    ) %>%
    filter(ensembl_gene_id.x != ensembl_gene_id.y) %>%
    dplyr::select(-contains("gene"))
iM30_edges <- gene_ppi %>%
    left_join(
        annotations %>%
            dplyr::select(ensembl_gene_id, ensembl_transcript_id),
        by = c("ENSG1" = "ensembl_gene_id")
    ) %>%
    rename(ENST1 = ensembl_transcript_id) %>%
    left_join(
        annotations %>%
            dplyr::select(ensembl_gene_id, ensembl_transcript_id),
        by = c("ENSG2" = "ensembl_gene_id")
    ) %>%
    rename(ENST2 = ensembl_transcript_id) %>%
    filter(ENST1 %in% iM30$name & ENST2 %in% iM30$name) %>%
    dplyr::select(ENST1, ENST2) %>%
    distinct() %>%
    filter(ENST1 != ENST2) %>%
    mutate(
        weight = mapply(
            function(f1, f2) {
                cor(iexpr[f1, ], iexpr[f2, ])
            }, 
            ENST1, ENST2
        )
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(
                ensembl_gene_id, ensembl_transcript_id
            ),
        by = c("ENST1" = "ensembl_transcript_id")
    ) %>%
    left_join(
        annotations %>%
            dplyr::select(
                ensembl_gene_id, ensembl_transcript_id
            ),
        by = c("ENST2" = "ensembl_transcript_id")
    ) %>%
    filter(ensembl_gene_id.x != ensembl_gene_id.y) %>%
    dplyr::select(-contains("gene"))

gM1_edges <- gM1_edges %>%
    filter(ENSG1 %in% asd_targets | ENSG2 %in% asd_targets)
iM1_edges <- iM1_edges %>%
    filter(ENST1 %in% asd_targets | ENST2 %in% asd_targets)
iM30_edges <- iM30_edges %>%
    filter(ENST1 %in% asd_targets | ENST2 %in% asd_targets)
    
gM1_vertices <- vertices %>%
    filter(`Module Label` == 1 & Data == "Gene") %>%
    filter(
        name %in% unlist(c(gM1_edges$ENSG1, gM1_edges$ENSG2)) |
            Target %in% c("ASD", "Both")
    )
iM1_vertices <- vertices %>%
    filter(`Module Label` == 1 & Data == "Isoform") %>%
    filter(
        name %in% unlist(c(iM1_edges$ENST1, iM1_edges$ENST2)) |
            Target %in% c("ASD", "Both")
    )
iM30_vertices <- vertices %>%
    filter(`Module Label` == 30 & Data == "Isoform") %>%
    filter(
        name %in% unlist(c(iM30_edges$ENST1, iM30_edges$ENST2)) |
            Target %in% c("ASD", "Both")
    )
write_tsv(gM1_edges, "data/ppi/genePPI_gM1_edges.tsv")
write_tsv(gM1_vertices, "data/ppi/genePPI_gM1_vertices.tsv")
write_tsv(iM1_edges, "data/ppi/genePPI_iM1_edges.tsv")
write_tsv(iM1_vertices, "data/ppi/genePPI_iM1_vertices.tsv")
write_tsv(iM30_edges, "data/ppi/genePPI_iM30_edges.tsv")
write_tsv(iM30_vertices, "data/ppi/genePPI_iM30_vertices.tsv")

# writexl::write_xlsx(
#     list(
#         gM1 = gM1_edges, iM1 = iM1_edges, iM30 = iM30_edges, Vertices = vertices
#     ),
#     "data/ppi/GenePPIAnnotated.xlsx"
# )

build_graph <- function(in_edges, vertices, module_colour) {
    this_graph <- graph_from_data_frame(
        setNames(in_edges, c("F1", "F2", "weight")), 
        vertices = vertices %>%
            filter(name %in% c(in_edges[, 1], in_edges[, 2]))
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

pdf("data/ppi/GenePPI.pdf", width = 16, height = 12)
print(lapply(
    seq(0.9, 0.99, 0.01),
    function(x) {
        tryCatch({
            build_graph(
                in_edges = gM1_edges %>%
                    filter(
                        (ENSG1 %in% asd_targets | ENSG2 %in% asd_targets) & 
                            abs(weight) >= x
                    ),
                vertices = vertices,
                module_colour = "turquoise"
            ) +
                labs(
                    title = paste("gM1, PCC >= ", x)
                )
        }, error = function(e) {
            return(NULL)
        })
    }
))
print(lapply(
    seq(0.9, 0.99, 0.01),
    function(x) {
        tryCatch({
            build_graph(
                in_edges = iM1_edges %>%
                    filter(
                        (ENST1 %in% asd_targets | ENST2 %in% asd_targets) & 
                            abs(weight) >= x
                    ),
                vertices = vertices,
                module_colour = "turquoise"
            ) +
                labs(
                    title = paste("iM1, PCC >= ", x)
                )
        }, error = function(e) {
            return(NULL)
        })
    }
))
print(
    build_graph(
        in_edges = iM30_edges %>%
            filter(ENST1 %in% asd_targets | ENST2 %in% asd_targets),
        vertices = vertices,
        module_colour = "steelblue"
    ) +
        labs(
            title = paste("iM30")
        )
)

dev.off()
