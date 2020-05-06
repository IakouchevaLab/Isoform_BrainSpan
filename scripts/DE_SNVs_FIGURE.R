library(tidyverse)
library(cowplot)

# DE_SNVs_2.R
de_iso_gene_lof_spec_enrichment_barplot <- readRDS(
    "data/figures/DE_Enr_LoFTargets_Bar.rds"
)

# DE_SNVs.R
diff_asd_iso_targets_breakdown <- readRDS(
    "data/isoforms/figures/DiffASDIsoLoFTargetFreq.rds"
)
heat_and_dend <- readRDS(
    "data/isoforms/figures/DiffASDIsoLoFTargetHeat.rds"
)
# DE_SNVs_4.R
# consequence_breakdown_plot <- readRDS("data/figures/LOF_DE_Breakdown.rds")
# CompareSNVTargetExpressions.R
compare_means <- readRDS(
    "data/isoforms/figures/SNVTargets_tpmWilcoxCompare.rds"
)
# compare_means_cortex <- readRDS(
#     "data/isoforms/figures/SNVTargets_tpmWilcoxCompareCortex.rds"
# )
# SelectedExpressionProfiles.R
expr_plot_combine <- readRDS("data/figures/DiffASDIsoTargets_ExprPltsCombine.rds")

de_snv_figure <- plot_grid(
    compare_means +
        theme(
            axis.text.x = element_text(size = 20, angle = 30, hjust = 1),
            plot.margin = margin(b = 5),
            rect = element_rect(fill = NA, colour = NA),
            legend.position = "top"
        ),
    plot_grid(
        plot_grid(
            diff_asd_iso_targets_breakdown +
                labs(y = "Proportion of Isoforms") +
                scale_y_continuous(expand = c(0, 0)) +
                theme(
                    text = element_text(size = 20),
                    axis.title.x = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1),
                    legend.position = "top"
                ),
            heat_and_dend,
            ncol = 2
        ),
        expr_plot_combine +
            theme(
                text = element_text(size = 20),
                legend.position = "none",
                plot.margin = margin()
            ),
        nrow = 2, rel_heights = c(1, 0.8), align = "v", axis = "lr"
    ),
    ncol = 2, nrow = 1, rel_widths = c(0.5, 1),
    labels = "AUTO", label_size = 40, label_fontface = "bold"
)

ggsave(
    filename = "data/figures/DE_SNV_FIGURE.pdf",
    plot = de_snv_figure,
    device = "pdf", width = 30, height = 18
)

