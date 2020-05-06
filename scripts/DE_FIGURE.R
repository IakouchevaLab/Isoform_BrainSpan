# Combine and lay out various DE figures

library(tidyverse)
library(cowplot)
library(gtable)
source("scripts/utility/plotting.R")

# DE_PostAnalysis.R
de_frequency_combined_plot_untrans <- readRDS(
    "data/figures/DEFrequencyCombinedPlot_untrans.rds"
)
de_frequency_combined_plot_untrans_noPrePost <- readRDS(
    "data/figures/DEFrequencyCombinedPlot_untrans_noPrePost.rds"
)
# DE_Enrichment_Analysis_Fisher.R
ct_list_spec_combined <- readRDS(
    "data/figures/DE_FisherEnrich_CTandList_Combined.rds"
)
ggsave(
    filename = "PaperFigures/de_frequency_combined_plot_untrans.pdf",
    plot = de_frequency_combined_plot_untrans,
    device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "PaperFigures/de_frequency_combined_plot_untrans_noPrePost.pdf",
    plot = de_frequency_combined_plot_untrans_noPrePost,
    device = "pdf", width = 16, height = 12
)
ggsave(
    filename = "data/figures/ct_list_spec_combined.pdf",
    plot = ct_list_spec_combined + theme(axis.text.x = element_text(angle = 90, hjust = 1)),
    device = "pdf", width = 16, height = 12
)

de_figure <- plot_grid(
    de_frequency_combined_plot +
        scale_x_discrete(position = "top") +
        scale_y_continuous(
            trans = "log1p", breaks = c(100, 1000, 10000)
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 0),
            axis.title.y = element_blank()
        ),
    ct_list_spec_combined +
        guides(
            fill = guide_colourbar(
                title = expression("-log"["10"]*"(adj. P-Value)"),
                direction = "vertical"
            )
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 24),
            panel.spacing = unit(0, units = "lines"),
            legend.position = "right",
            legend.key.height = unit(3, "lines"),
            legend.key.width = unit(2, "lines")
        ),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(0.5, 0.8), 
    labels = "AUTO", label_size = 40, label_fontface = "bold"
)
# de_figure

ggsave(
    filename = "data/figures/DE_FIGURE.pdf",
    plot = de_figure,
    device = "pdf", width = 16, height = 16
)

de_figure_untrans <- plot_grid(
    de_frequency_combined_plot_untrans +
        scale_x_discrete(position = "top") +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 0),
            axis.title.y = element_blank()
        ),
    ct_list_spec_combined +
        guides(
            fill = guide_colourbar(
                title = expression("-log"["10"]*"(adj. P-Value)"),
                direction = "vertical"
            )
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 24),
            panel.spacing = unit(0, units = "lines"),
            legend.position = "right",
            legend.key.height = unit(3, "lines"),
            legend.key.width = unit(2, "lines")
        ),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(0.5, 0.8), 
    labels = "AUTO", label_size = 40, label_fontface = "bold"
)
# de_figure

ggsave(
    filename = "data/figures/DE_FIGURE_untrans.pdf",
    plot = de_figure_untrans,
    device = "pdf", width = 16, height = 16
)



de_figure_untrans_noPrePost <- plot_grid(
    de_frequency_combined_plot_untrans_noPrePost +
        scale_y_continuous(position = "left") +
        scale_x_discrete(position = "top") +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 0),
            axis.title.y = element_blank()
        ),
    ct_list_spec_combined +
        guides(
            fill = guide_colourbar(
                title = expression("-log"["10"]*"(adj. P-Value)"),
                direction = "vertical"
            )
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 24),
            panel.spacing = unit(0, units = "lines"),
            legend.position = "right",
            legend.key.height = unit(3, "lines"),
            legend.key.width = unit(2, "lines")
        ),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(0.5, 0.8), 
    labels = "AUTO", label_size = 40, label_fontface = "bold"
)
# de_figure

ggsave(
    filename = "data/figures/DE_FIGURE_untrans_noPrePost.pdf",
    plot = de_figure_untrans_noPrePost,
    device = "pdf", width = 16, height = 16
)
