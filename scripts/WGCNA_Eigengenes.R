library(tidyverse)
source("scripts/utility/plot_expression.R")

cell_types <- stack(c(
    "Oligodendrocytes" = "IsoformME6",
    "NPCs" = "IsoformME10",
    "Microglia" = "IsoformME35",
    "Interneurons" = "IsoformME17",
    "Excitatory neurons" = "IsoformME2",
    "Astrocytes" = "IsoformME25"
))
colnames(cell_types) <- c("ME", "CellType")
metadata <- read_tsv("data/metadata.tsv")
module_map <- bind_rows(
    read_tsv("data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        distinct(module_label, module_colour) %>%
        mutate(Data = "Gene"),
    read_tsv("data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        distinct(module_label, module_colour) %>%
        mutate(Data = "Isoform")
) %>%
    mutate(Module = paste0(Data, "ME", module_label))
MEs <- bind_rows(
    readRDS("data/genes/Networks/Network_DS2_MM20.rds")$MEs %>%
        as.data.frame() %>%
        mutate(Sample = rownames(.), Data = "Gene") %>%
        gather(starts_with("ME"), key = "ME", value = "Expression"),
    readRDS("data/isoforms/Networks/Network_DS2_MM20.rds")$MEs %>%
        as.data.frame() %>%
        mutate(Sample = rownames(.), Data = "Isoform") %>%
        gather(starts_with("ME"), key = "ME", value = "Expression")
) %>%
    mutate(Module = paste0(Data, ME)) %>%
    left_join(module_map, by = c("Data", "Module")) %>%
    left_join(cell_types, by = c("Module" = "ME")) %>%
    left_join(metadata, by = "Sample")

celltype_me_plot <- (ggplot(
    data = filter(MEs, !is.na(CellType)),
    mapping = aes(
        x = Days, y = Expression, group = CellType, colour = CellType
    )
) +
    geom_line(stat = "smooth", method = "loess", size = 2) +
    labs(
        x = "Days", y = "Module Eigengene Expression"
    ) +
    scale_colour_brewer(
        palette = "Set1"
    ) +
    theme(
        text = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.justification = c(1, 0),
        legend.position = c(0.95, 0.05)
    )) %>%
    expression_plot_fill(metadata)
saveRDS(
    celltype_me_plot,
    "data/figures/CellTypeMEs.rds"
)
