library(tidyverse)

zhong_celltypes <- readxl::read_xlsx(
    file.path(
        "data/source/CuratedLists",
        "ZhongNature2018_humanPreFrontalCortexCellTypes.xlsx"
    ),
    sheet = 1, skip = 4
)[, c("cluster", "gene")]
wang_celltypes <- readxl::read_xlsx(
    file.path(
        "data/source/CuratedLists",
        "WangScience2018CellTypesPsychencode.xlsx"
    ),
    sheet = 1, skip = 1
) %>%
    mutate(
        cluster = sapply(
            Cluster,
            function(clst) {
                switch(
                    str_match(
                        clst,
                        "Ex|In|Astro|OPC|Oligo|Endo|Per|Microglia"
                    )[, 1],
                    "Ex" = "Excitatory neurons",
                    "In" = "Interneurons",
                    "Astro" = "Astrocytes",
                    "Oligo" = "Oligodendrocytes",
                    "OPC" = "OPC", 
                    "Endo" = "Endothelial",
                    "Per" = "Pericytes",
                    "Microglia" = "Microglia"
                )
            }
        )
    ) %>%
    rename(gene = Gene) %>%
    dplyr::select(gene, cluster)

intersect_gene_clusters <- inner_join(
    zhong_celltypes,
    wang_celltypes,
    by = "gene", 
    suffix = c("_zhong", "_wang")
) %>% 
    distinct() %>%
    dplyr::select(gene, cluster_zhong, cluster_wang)
