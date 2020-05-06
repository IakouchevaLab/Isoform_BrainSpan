library(tidyverse)
library(cowplot)

dend_plot <- readRDS("data/figures/GeneIsoformMEDendrogram.rds")
module_assoc_plot <- readRDS("data/figures/ModuleAssociationTiles.rds")
list_enrichment_plot <- readRDS("data/figures/ModuleListEnrichmentFisher.rds")
# hub_networks <- readRDS("data/figures/WGCNA_networks.rds")
# actual_v_expect_vp <- readRDS("data/figures/WGCNA_PositionalImpact_ASDLOF.rds")
celltype_me_plot <- readRDS("data/figures/CellTypeMEs.rds")
# WGCNA_SNVs.R
vpm_diff <- readRDS("data/figures/VPM_byStatus_lollipop.rds")
m1_gost <- readRDS("data/figures/GIM1_ENR.rds")
im30_gost <- readRDS("data/figures/IM30_ENR.rds")

sigModule_funEnrich_combine <- plot_grid(
    m1_gost + 
        theme(
            text = element_text(size = 20),
            legend.position = "none"
        ), 
    im30_gost + 
        theme(
            text = element_text(size = 20),
            legend.position = "none"
        ),
    get_legend(
        m1_gost +
            guides(size = guide_legend(title = "Feature\nIntersection")) +
            theme(
                legend.title = element_text(size = 18)
            )
    ),
    nrow = 1, rel_widths = c(1, 1, 0.3)
)
sigModule_funEnrich_combine

Modules_Figures <- plot_grid(
    plot_grid(
        dend_plot +
            coord_cartesian(clip = "off") +
            theme_nothing() +
            theme(
                axis.text = element_blank(),
                plot.margin = margin()
            ),
        module_assoc_plot +
            theme(
                text = element_text(size = 20),
                strip.text = element_text(size = 20),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_blank(),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                legend.key.height = unit(1, "lines"),
                plot.margin = margin()
            ),
        list_enrichment_plot +
            theme(
                strip.text = element_text(size = 20),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_blank(),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                legend.key.height = unit(1, "lines"),
                panel.spacing = unit(0, units = "lines"),
                plot.margin = margin()
            ),
        celltype_me_plot +
            theme(
                axis.title = element_text(size = 24),
                axis.text = element_text(size = 24),
                legend.background = element_rect(
                    fill = "white", colour = "white", size = 1
                ),
                legend.direction = "horizontal",
                legend.text = element_text(size = 20),
                legend.key.size = unit(2, "lines")
            ),
        nrow = 1, align = "h", axis = "tb",
        rel_widths = c(0.15, 0.4, 0.45, 0.75)
    ),
    plot_grid(
        vpm_diff +
            theme(
                plot.margin = margin(l = 2, unit = "lines")
            ),
        sigModule_funEnrich_combine,
        nrow = 1, rel_widths = c(0.75, 1)
    ),
    nrow = 2, rel_heights = c(1, 0.8), axis = "lr", align = "v"
)
# Modules_Figures

ggsave(
    filename = "data/figures/WGCNA_FIGURE.pdf",
    plot = Modules_Figures,
    device = "pdf", width = 28, height = 20, useDingbats = FALSE
)
