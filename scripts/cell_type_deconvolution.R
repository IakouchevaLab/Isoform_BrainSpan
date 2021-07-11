library(tidyverse)

cell_type_props <- read_csv("data/estimated_props_Bisque_lake_level1.csv") %>%
    rename(Sample = X1)

cell_type_long <- cell_type_props %>%
    pivot_longer(
        cols = c(Ast, End, Ex, In, Mic, Oli, OPC, Per),
        names_to = "CellType",
        values_to = "Fraction"
    )

plt <- ggplot(
    data = cell_type_long,
    mapping = aes(
        x = Sample, y = Fraction, fill = CellType
    )
) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +
    theme_bw() +
    theme(
        text = element_text(size = 12),
        axis.text.y = element_text(size = 5, hjust = 1),
        legend.position = "top"
    )
ggsave(
    "data/figures/cell_type_deconvolution.pdf",
    plt,
    device = "pdf", width = 10, height = 40, units = "in"
)
