library(tidyverse)
library(cowplot)

read_sheets <- function(path) {
    lapply(
        setNames(nm = readxl::excel_sheets(path)),
        function(sheet) {
            readxl::read_xlsx(path, sheet = sheet)
        }
    )
}

gene_tables <- list(
    de = read_sheets("data/genes/Enrichments/DE_GOST.xlsx"),
    de_spec = read_sheets("data/genes/Enrichments/DE_GOST_SPEC.xlsx"),
    de_spec_maxSize1000 = lapply(
        read_sheets("data/genes/Enrichments/DE_GOST_SPEC.xlsx"),
        function(enr) {
            filter(enr, term_size <= 1000)
        }
    ),
    de_spec_maxSize400 = lapply(
        read_sheets("data/genes/Enrichments/DE_GOST_SPEC.xlsx"),
        function(enr) {
            filter(enr, term_size <= 400)
        }
    ),
    lof_de = read_sheets("data/genes/Enrichments/LOF_DE_GOST.xlsx"),
    lof_de_spec = read_sheets("data/genes/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
    lof_de_spec_maxSize1000 = lapply(
        read_sheets("data/genes/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
        function(enr) {
            filter(enr, term_size <= 1000)
        }
    ),
    lof_de_spec_maxSize400 = lapply(
        read_sheets("data/genes/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
        function(enr) {
            filter(enr, term_size <= 400)
        }
    )
)

isoform_tables <- list(
    de = read_sheets("data/isoforms/Enrichments/DE_GOST.xlsx"),
    de_spec = read_sheets("data/isoforms/Enrichments/DE_GOST_SPEC.xlsx"),
    de_spec_maxSize1000 = lapply(
        read_sheets("data/isoforms/Enrichments/DE_GOST_SPEC.xlsx"),
        function(enr) {
            filter(enr, term_size <= 1000)
        }
    ),
    de_spec_maxSize400 = lapply(
        read_sheets("data/isoforms/Enrichments/DE_GOST_SPEC.xlsx"),
        function(enr) {
            filter(enr, term_size <= 400)
        }
    ),
    lof_de = read_sheets("data/isoforms/Enrichments/LOF_DE_GOST.xlsx"),
    lof_de_spec = read_sheets("data/isoforms/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
    lof_de_spec_maxSize1000 = lapply(
        read_sheets("data/isoforms/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
        function(enr) {
            filter(enr, term_size <= 1000)
        }
    ),
    lof_de_spec_maxSize400 = lapply(
        read_sheets("data/isoforms/Enrichments/LOF_DE_SPEC_GOST.xlsx"),
        function(enr) {
            filter(enr, term_size <= 400)
        }
    )
)

plot_gost <- function(enr, n, title = "") {
    ggplot(
        data = enr %>%
            arrange(p_value) %>%
            mutate(Order = seq(nrow(.), 1)) %>%
            top_n(n, Order),
        mapping = aes(
            x = -log10(p_value), 
            y = reorder(str_wrap(term_name, 40), Order),
            size = intersection_size, fill = -log10(p_value)
        )
    ) +
        geom_point(shape = 21) +
        geom_vline(xintercept = -log10(0.05), colour = "red") +
        labs(
            x = expression("-log"["10"]*"FDR"),
            title = title
        ) +
        guides(
            fill = guide_colourbar(title = expression("-log"["10"]*"FDR")),
            size = guide_legend(title = "Overlap")
        ) +
        scale_fill_gradient(low = "white", high = "red") +
        scale_size_continuous(range = c(5, 15)) +
        theme_bw() +
        theme(
            text = element_text(size = 24),
            axis.title.y = element_blank()
        )
}
# plot_gost(isoform_tables$de$P02P03 %>% top_n(5, -p_value))

gene_tables_plots <- lapply(
    names(gene_tables),
    function(gt) {
        lapply(
            names(gene_tables[[gt]]),
            function(period) {
                plot_gost(
                    gene_tables[[gt]][[period]],
                    10,
                    title = paste0(gt, ", ", period)
                )
            }
        )
    }
)
pdf("data/genes/Enrichments/DE_GOST_Figures.pdf", width = 12, height = 12)
invisible(lapply(gene_tables_plots, print))
dev.off()

isoform_tables_plots <- lapply(
    names(isoform_tables),
    function(it) {
        lapply(
            names(isoform_tables[[it]]),
            function(period) {
                plot_gost(
                    isoform_tables[[it]][[period]],
                    10,
                    title = paste0(it, ", ", period)
                )
            }
        )
    }
)
pdf("data/isoforms/Enrichments/DE_GOST_Figures.pdf", width = 12, height = 12)
invisible(lapply(isoform_tables_plots, print))
dev.off()

################################################################################
# Plot all tables together                                                     #
################################################################################

all_names <- intersect(names(gene_tables), names(isoform_tables))

combined_plots <- lapply(
    setNames(nm = all_names),
    function(table_name) {
        contrasts <- intersect(
            names(gene_tables[[table_name]]), 
            names(isoform_tables[[table_name]])
        )
        lapply(
            setNames(nm = contrasts),
            function(contrast) {
                enr <- bind_rows(
                    gene_tables[[table_name]][[contrast]] %>%
                        mutate(Data = "Gene") %>%
                        arrange(p_value) %>%
                        mutate(Order = seq(nrow(.), 1)),
                    isoform_tables[[table_name]][[contrast]] %>%
                        mutate(Data = "Isoform") %>%
                        arrange(p_value) %>%
                        mutate(Order = seq(nrow(.), 1))
                ) %>%
                    group_by(Data) %>%
                    top_n(5, Order) %>%
                    ungroup()
                ggplot(
                    data = enr,
                    mapping = aes(
                        x = -log10(p_value), 
                        y = reorder(str_wrap(term_name, 20), Order),
                        size = intersection_size, fill = -log10(p_value)
                    )
                ) +
                    facet_grid(Data ~ ., scales = "free_y") +
                    geom_point(shape = 21) +
                    # geom_vline(xintercept = -log10(0.05), colour = "red") +
                    labs(
                        x = expression("-log"["10"]*"FDR"),
                        title = paste(table_name, contrast)
                    ) +
                    guides(
                        fill = guide_colourbar(
                            title = expression("-log"["10"]*"FDR")
                        ),
                        size = guide_legend(title = "Overlap")
                    ) +
                    scale_fill_gradient(low = "white", high = "red") +
                    scale_size_continuous(range = c(5, 15)) +
                    theme_bw() +
                    theme(
                        text = element_text(size = 24),
                        axis.title.y = element_blank()
                    )
            }
        )
    }
)
# combined_plots[[1]][[1]]
# pdf("data/figures/DE_GOST_Figures.pdf", width = 16, height = 16)
# invisible(lapply(combined_plots, function(plt) invisible(lapply(plt, print))))
# dev.off()

################################################################################
# Final Plots                                                                  #
################################################################################

# DE Spec

de_spec_finalPlots <- lapply(
    list(
        combined_plots$de_spec_maxSize1000$P04P05 +
            labs(title = "P04P05") +
            scale_x_continuous(
                breaks = c(5, 10)
            ),
        combined_plots$de_spec_maxSize1000$P07P08 +
            labs(title = "P07P08"),
        combined_plots$de_spec_maxSize1000$P08P09 +
            labs(title = "P08P09")
    ),
    function(plt) {
        plt + 
            scale_fill_gradient(
                low = "white", high = "red", limits = c(0, 25)
            ) +
            scale_size_continuous(
                limits = c(0, 200)
            ) +
            labs(
                x = expression("-log"["10"]*"(FDR)")
            ) +
            theme(
                text = element_text(size = 24),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20)
            )
    }
)
final_plot_deSpec <- plot_grid(
    plot_grid(
        plotlist = lapply(
            de_spec_finalPlots,
            function(plt) plt + theme(legend.position = "none")
        ),
        nrow = 1, 
        labels = "AUTO", label_size = 40, label_fontface = "bold"
    ),
    get_legend(de_spec_finalPlots[[1]]),
    ncol = 2, rel_widths = c(1, 0.12)
)
ggsave(
    filename = "data/figures/DE_GOST_DESpec_Selected.pdf",
    final_plot_deSpec,
    device = "pdf", width = 16, height = 12,
    useDingbats = FALSE
)

# LOF DE SPEC

lof_de_spec_finalPlots <- lapply(
    list(
        combined_plots$lof_de_spec_maxSize1000$P04P05 +
            labs(title = "P04P05"),
        combined_plots$lof_de_spec_maxSize1000$P07P08 +
            labs(title = "P07P08"),
        combined_plots$lof_de_spec_maxSize1000$P08P09 +
            labs(title = "P08P09")
    ),
    function(plt) {
        plt + 
            scale_fill_gradient(
                low = "white", high = "red", limits = c(0, 5)
            ) +
            scale_size_continuous(
                limits = c(0, 15)
            ) +
            labs(
                x = expression("-log"["10"]*"(FDR)")
            ) +
            theme(
                text = element_text(size = 24),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 10)
            )
    }
)
final_plot_LofDeSpec <- plot_grid(
    plot_grid(
        plotlist = lapply(
            lof_de_spec_finalPlots,
            function(plt) plt + theme(legend.position = "none")
        ),
        nrow = 1, 
        labels = "AUTO", label_size = 40, label_fontface = "bold"
    ),
    get_legend(lof_de_spec_finalPlots[[1]]),
    ncol = 2, rel_widths = c(1, 0.12)
)
ggsave(
    filename = "data/figures/LOF_DE_GOST_DESpec_Selected.pdf",
    final_plot_LofDeSpec,
    device = "pdf", width = 16, height = 8
)

# PREPOST

prepost_finalPlots <- lapply(
    list(
        combined_plots$de_spec_maxSize1000$PrePost +
            labs(title = "PrePost, All DE"),
        combined_plots$lof_de_spec_maxSize1000$PrePost +
            labs(title = "PrePost, LoF targets")
    ),
    function(plt) {
        plt + 
            scale_fill_gradient(
                low = "white", high = "red", limits = c(0, 15)
            ) +
            scale_size_continuous(
                limits = c(0, 155)
            ) +
            labs(
                x = expression("-log"["10"]*"(FDR)")
            ) +
            theme(
                text = element_text(size = 24),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 10)
            )
    }
)
final_plot_PrePost <- plot_grid(
    plot_grid(
        plotlist = lapply(
            prepost_finalPlots,
            function(plt) plt + theme(legend.position = "none")
        ),
        nrow = 1, 
        labels = "AUTO", label_size = 40, label_fontface = "bold"
    ),
    get_legend(prepost_finalPlots[[1]]),
    ncol = 2, rel_widths = c(1, 0.12)
)
final_plot_PrePost
ggsave(
    filename = "data/figures/LOF_DE_PREPOST_GOST_Selected.pdf",
    final_plot_PrePost,
    device = "pdf", width = 16, height = 8
)
