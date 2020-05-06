library(tidyverse)
library(cowplot)
library(igraph)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", sep = ",",
    header = TRUE, row.names = 1
)

yang_ppi <- readxl::read_xlsx(
    "data/source/CuratedLists/Yang et al. - 2016 - Widespread Expansion of Protein Interaction Capabilities by Alternative Splicing(2).xlsx",
    sheet = 2
) %>%
    filter(Interaction_Found == "positive") %>%
    mutate(
        Interaction = mapply(
            function(gs, is) {
                paste0(gs, is, is, gs)
            },
            Gene_Symbol, Interactor_Symbol
        )
    )
corominas_ppi

gene_module_assignments <- lapply(
    readxl::excel_sheets(
        "data/SupplementaryTables/Supplementary Table 7.xlsx"
    )[-1],
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = sheet
        )
    }
) %>%
    bind_rows()

isoform_module_assignments <- lapply(
    readxl::excel_sheets(
        "data/SupplementaryTables/Supplementary Table 8.xlsx"
    )[-1],
    function(sheet) {
        readxl::read_xlsx(
            "data/SupplementaryTables/Supplementary Table 8.xlsx", sheet = sheet
        )
    }
) %>%
    bind_rows()

geneExpr <- readRDS("data/RegressGeneCounts.rds")
isoformExpr <- readRDS("data/RegressIsoformCounts.rds")

# iM1
iM1_Members <- isoform_module_assignments %>%
    filter(`Module Label` == 1) %>%
    pull(`Ensembl Transcript ID`)
iM1_Hubs <- isoform_module_assignments %>%
    filter(`Module Label` == 1) %>%
    arrange(desc(kME)) %>%
    top_n(20, kME) %>%
    pull(`Ensembl Transcript ID`)
iM1_Coex <- expand.grid(
    Hub = as.character(iM1_Hubs),
    Member = as.character(iM1_Members)
) %>%
    mutate(
        Hub = as.character(Hub),
        Member = as.character(Member)
    ) %>%
    filter(Hub != Member) %>%
    mutate(
        Correlation = mapply(
            function(hub, member) {
                cor(
                    isoformExpr[hub, ], isoformExpr[member, ]
                )
            }, Hub, Member
        )
    ) %>%
    group_by(Hub) %>%
    filter(Correlation == max(Correlation)) %>%
    ungroup() %>%

iM1_Coex_annotate <- iM1_Coex %>%
    mutate(
        Hub_Symbol = sapply(Hub, function(hub) {
            annotations %>%
                distinct(ensembl_transcript_id, external_gene_id) %>%
                filter(ensembl_transcript_id == hub) %>%
                pull(external_gene_id)
        }),
        Member_Symbol = sapply(Member, function(member) {
            annotations %>%
                distinct(ensembl_transcript_id, external_gene_id) %>%
                filter(ensembl_transcript_id == member) %>%
                pull(external_gene_id)
        })
    )

iM1_Hub_Graph <- graph_from_edgelist(
    as.matrix(iM1_Coex[, c(1, 2)]), directed = FALSE
)
iM1_Layout <- data.frame(layout_with_fr(iM1_Hub_Graph))
iM1_Edgeset <- get.edgelist(iM1_Hub_Graph, names = FALSE)
iM1_Edges <- data.frame(
    iM1_Layout[iM1_Edgeset[, 1], ], iM1_Layout[iM1_Edgeset[, 2], ]
)
colnames(iM1_Edges) <- c("X1", "Y1", "X2", "Y2")
iM1_Hub_Graph_Table <- cbind(
    iM1_Layout, 
    data.frame(Gene = iM1_Hubs[["Gene Symbol"]])
)
iM1_Interactions <- left_join(
    iM1_Edges, 
    iM1_Hub_Graph_Table,
    by = c("X1" = "X1", "Y1" = "X2")
) %>%
    rename(Gene1 = Gene) %>%
    left_join(
        ., 
        iM1_Hub_Graph_Table,
        by = c("X2" = "X1", "Y2" = "X2")
    ) %>%
    rename(Gene2 = Gene) %>%
    mutate(PPI = paste0(Gene1, Gene2) %in% isoform_ppi$Interaction)

ggplot() +
    geom_segment(
        data = iM1_Interactions,
        mapping = aes(
            x = X1, y = Y1, xend = X2, yend = Y2
        ), 
        size = 0.5, colour = "lightgrey"
    ) +
    geom_point(
        data = iM1_Hub_Graph_Table,
        mapping = aes(
            x = X1, y = X2
        ),
        shape = 21, size = 5, fill = "turquoise"
    ) +
    theme_nothing()



iM1_Hubs_ADJ <- WGCNA::adjacency(
    t(isoformExpr[iM1_Hubs[["Ensembl Transcript ID"]], ]), 
    type = "signed", power = 3
)
iM1_Hub_Graph <- graph.adjacency(
    as.matrix(iM1_Hubs_ADJ), mode = "undirected", weighted = TRUE, diag = FALSE
)
iM1_Layout <- data.frame(layout_with_fr(iM1_Hub_Graph))
iM1_Edgeset <- get.edgelist(iM1_Hub_Graph, names = FALSE)
iM1_Edges <- data.frame(
    iM1_Layout[iM1_Edgeset[, 1], ], iM1_Layout[iM1_Edgeset[, 2], ]
)
colnames(iM1_Edges) <- c("X1", "Y1", "X2", "Y2")
iM1_Hub_Graph_Table <- cbind(iM1_Layout, data.frame(
    Gene = iM1_Hubs[["Gene Symbol"]])
)
iM1_Interactions <- left_join(
    iM1_Edges, 
    iM1_Hub_Graph_Table,
    by = c("X1" = "X1", "Y1" = "X2")
) %>%
    rename(Gene1 = Gene) %>%
    left_join(
        ., 
        iM1_Hub_Graph_Table,
        by = c("X2" = "X1", "Y2" = "X2")
    ) %>%
    rename(Gene2 = Gene) %>%
    mutate(PPI = paste0(Gene1, Gene2) %in% isoform_ppi$Interaction)
    
ggplot() +
    geom_segment(
        data = iM1_Interactions,
        mapping = aes(
            x = X1, y = Y1, xend = X2, yend = Y2
        ), 
        size = 0.5, colour = "lightgrey"
    ) +
    geom_point(
        data = iM1_Hub_Graph_Table,
        mapping = aes(
            x = X1, y = X2
        ),
        shape = 21, size = 5, fill = "turquoise"
    ) +
    theme_nothing()
