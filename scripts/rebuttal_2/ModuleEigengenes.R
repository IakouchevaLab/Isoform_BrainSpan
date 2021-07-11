library(tidyverse)

source("scripts/utility/plot_expression.R")

#' Add standard geoms/parts to an expression profile plot
#'
#' @param plt Input plot
#' @param metadata Metadata for period definitions 
#'                 (should include "Days" and "Period" fields)
#' 
#' @return Filled expression profile ggplot object
expression_plot_fill <- function(plt, metadata) {
    plt +
        geom_vline(
            xintercept = get_period_breaks(metadata)[[6]], colour = "red"
        )  +
        scale_x_continuous(
            trans = "log1p",
            breaks = c(50, 100, 1000, 10000),
            expand = c(0, 0)
        )
}

metadata <- read_tsv("data/metadata.tsv")

module_info <- read_tsv("data/ModuleInfo.tsv") %>%
    select(
        module_label, module_colour, Data
    ) %>%
    distinct() %>%
    rename(Module = module_label, Colour = module_colour, group = Data)

isoform_modules <- readRDS("data/isoforms/Networks/Network_DS2_MM20.rds")
isoform_me <- isoform_modules$MEs
gene_modules <- readRDS("data/genes/Networks/Network_DS2_MM20.rds")
gene_me <- gene_modules$MEs

module_annotations <- bind_rows(
    tibble(
        group = "Gene",
        ME = c("ME1", "ME2", "ME3", "ME4", "ME5", "ME6", "ME7", "ME8"),
        Annotation = c(
            "RNA processing/splicing", 
            "Trans-synaptic signaling",
            "Mitochondrion/translation",
            "Angiogenesis",
            "Cilium/cytoskeleton",
            "Embryonic morphogenesis",
            "Transmembrane transport",
            "NA"
        )
    ),
    tibble(
        group = "Isoform",
        ME = c("ME1" , "ME2" , "ME3" , "ME4" , "ME5" , "ME6" , "ME7" , "ME8" , 
               "ME9" , "ME10", "ME11", "ME12", "ME13", "ME14", "ME15", "ME16", 
               "ME17", "ME18", "ME19", "ME20", "ME21", "ME22", "ME23", "ME24", 
               "ME25", "ME26", "ME27", "ME28", "ME29", "ME30", "ME31", "ME32", 
               "ME33", "ME34", "ME35", "ME36", "ME37", "ME38", "ME39", "ME40", 
               "ME41", "ME42", "ME43", "ME44", "ME45", "ME46", "ME47", "ME48", 
               "ME49", "ME50", "ME51", "ME52", "ME53", "ME54", "ME55"),
        Annotation = c(
            "RNA processing/splicing",
            "Trans-synaptic signaling",
            "Protein catabolism/GTPase",
            "Mitochondrion/translation",
            "RNA processing/splicing",
            "Neurons ensheathment/actin/glia",
            "Angiogenesis",
            "Exocytosis",
            "Axonogenesis/brain development",
            "Mitosis/cell cycle",
            "Cilium/microtubule",
            "ATP metabolism",
            "Behavior/cognition/learning and memory",
            "Histone modification/carbohydrate",
            "Adhesion/locomotion/Wnt signaling",
            "Neuron projection development/actin",
            "Neuron/synapse",
            "Chemical synaptic transmission/neurogenesis",
            "Chromatin/histone",
            "Cytoskeleton/microtubule transport",
            "Translation",
            "Tricarboxylic acid transport",
            "Transmembrane/anion transport",
            "Synaptic signaling/neurotransmitter transport",
            "Catabolism",
            "Transcription/RNA processing/catabolism",
            "Neuron projection morphogenesis/cytoskeleton/Ras GTPase binding",
            "Interferon/cytokine signaling",
            "Protein localization/targeting",
            "Synapse organization/dendrite development",
            "Chromatin organization/binding/neuron migration",
            "Protein folding/chaperone",
            "Cation/ion transport",
            "Response to hormone/retinoic acid",
            "Angiogenesis/tube morphogenesis",
            "Immune response/cytokine/leukocyte",
            "Phospholipid binding/neurotransmitter",
            "Extracellular matrix organization",
            "Mesenchyme/embryonic morphogenesis",
            "Cell adhesion/neuron fate/cerebral cortex regionalization",
            "NA",
            "Synaptic signaling/exocytosis/neurotransmitter secretion",
            "Synaptic vesicle exocytosis/sodium:potassium exchange",
            "Chemical synaptic transmission/voltage-gated channel activity",
            "ATP-gated ion channel activity",
            "Microtubules/cilium",
            "Protein/histone demethylation",
            "Carbohydrate binding",
            "Cell cycle/nuclear transport/clathrin",
            "Amino acid/glutamate transport",
            "Kinase activity",
            "Protein targeting/localization",
            "Transmembrane transport",
            "Beta-tubulin/EGF/cadherin binding",
            "Innate immune response"
        )
    )
)

isoform_me_gene_me_correlation <- lapply(
    setNames(nm = colnames(isoform_me)),
    function(i) {
        sapply(
            setNames(nm = colnames(gene_me)),
            function(g) {
                cor(isoform_me[[i]], gene_me[[g]])
            }
        )
    }
)

isoform_max_correlations <- sort(sapply(isoform_me_gene_me_correlation, function(x) max(abs(x))))
i_plt_dat <- isoform_me %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols = starts_with("ME"), names_to = "ME", values_to = "Expression") %>%
    mutate(group = "Isoform") %>%
    mutate(MaxGeneCorrelation = isoform_max_correlations[ME]) %>%
    mutate(Distinct = MaxGeneCorrelation < 0.25) %>%
    mutate(Module = as.numeric(str_replace(ME, 'ME', ''))) %>%
    filter(Module != 0) %>%
    mutate(MaxCorrelation = isoform_max_correlations[ME]) %>%
    left_join(metadata, by = c("Sample" = "Sample")) %>%
    left_join(module_annotations, by = c("group", "ME")) %>%
    left_join(module_info, by = c("group", "Module")) %>%
    mutate(
        label_colour = mapply(
            function(x, m) {
                paste0(
                    m, " ( ",
                    paste(col2rgb(x)[, 1], collapse = ', '),
                    ")"
                )
            },
            Colour, Module
        )
    )
g_plt_dat <- gene_me %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols = starts_with("ME"), names_to = "ME", values_to = "Expression") %>%
    mutate(group = "Gene") %>%
    mutate(Distinct = FALSE) %>%
    mutate(Module = as.numeric(str_replace(ME, 'ME', ''))) %>%
    filter(Module != 0) %>%
    left_join(metadata, by = c("Sample" = "Sample")) %>%
    left_join(module_annotations, by = c("group", "ME")) %>%
    left_join(module_info, by = c("group", "Module")) %>%
    mutate(
        label_colour = mapply(
            function(x, m) {
                paste0(
                    m, " ( ",
                    paste(col2rgb(x)[, 1], collapse = ', '),
                    ")"
                )
            },
            Colour, Module
        )
    )

i_plts <- lapply(
    seq(1, 55, by = 8),
    function(i) {
        i_i_plt_dat <- i_plt_dat %>%
            filter(Module %in% seq(i, i + 7))
        expr_plt <- ggplot() +
            facet_wrap(~label_colour, ncol = 1) +
            geom_line(
                data = i_i_plt_dat,
                mapping = aes(
                    x = Days, y = Expression, group = ME
                ),
                stat = "smooth", method = "loess", size = 1
            ) +
            geom_text(
                data = i_i_plt_dat %>%
                    distinct(ME, Module, MaxGeneCorrelation),
                mapping = aes(
                    label = paste0("PCC:", sprintf("%.2f", MaxGeneCorrelation))
                ),
                x = log1p(100), y = 0.075, hjust = 0, size = 2, colour = "black"
            ) +
            scale_colour_manual(values = setNames(nm = module_info$Colour)) +
            #scale_y_continuous(limits = c(-0.08, 0.08)) +
            coord_cartesian(ylim = c(-0.1, 0.1)) +
            guides(colour = FALSE) +
            theme_bw() +
            theme(
                text = element_text(size = 20),
                plot.title = element_text(size = 20),
                strip.text.x = element_text(size = 3),
                strip.background = element_rect(fill = NA, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black"),
                panel.grid = element_blank(),
                strip.text = element_text(size = 10),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank()
            )
        annot_plt <- ggplot() +
            facet_wrap(~Module, ncol = 1) +
            geom_text(
                data = i_i_plt_dat,
                mapping = aes(
                    x = 0, y = 0, label = str_wrap(Annotation, width = 15)
                ),
                size = 2.5
            ) +
            coord_cartesian(ylim = c(-0.1, 0.1)) +
            theme_bw() +
            theme(
                text = element_text(size = 20),
                plot.title = element_blank(),
                strip.text = element_blank(),
                strip.background = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                plot.margin = margin(0, unit = "pt")
            )
        # expr_plt_tbl <- ggplot_gtable(ggplot_build(expr_plt))
        # strip_t <- which(grepl('strip-t', expr_plt_tbl$layout$name))
        # fills <- pull(arrange(distinct(i_i_plt_dat, Module, Colour), desc(Module)), Colour)
        # k <- 1
        # for (i in strip_t) {
        #     j <- which(grepl('rect', expr_plt_tbl$grobs[[i]]$grobs[[1]]$childrenOrder))
        #     expr_plt_tbl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
        #     k <- k+1
        # }
        combined_plt <- cowplot::plot_grid(
            plotlist = list(expr_plt, annot_plt), ncol = 2, align = "h",
            rel_widths = c(0.75, 0.8)
        )
        return(combined_plt)
    }
)
g_expr_plt <- ggplot() +
    facet_wrap(~label_colour, ncol = 1) +
    geom_line(
        data = g_plt_dat,
        mapping = aes(
            x = Days, y = Expression, group = ME
        ),
        stat = "smooth", method = "loess", size = 1
    ) +
    #scale_y_continuous(limits = c(-0.08, 0.08)) +
    coord_cartesian(ylim = c(-0.1, 0.1)) +
    scale_colour_manual(values = setNames(nm = module_info$Colour)) +
    guides(colour = FALSE) +
    # labs(title = "Gene Eigengenes") +
    theme_bw() +
    theme(
        text = element_text(size = 20),
        strip.text.x = element_text(size = 3),
        strip.background = element_rect(fill = NA, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        strip.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    )
# g_expr_plt_tbl <- ggplot_gtable(ggplot_build(g_expr_plt))
# strip_t <- which(grepl('strip-t', g_expr_plt_tbl$layout$name))
# fills <- pull(arrange(distinct(g_plt_dat, Module, Colour), desc(Module)), Colour)
# k <- 1
# for (i in strip_t) {
#     j <- which(grepl('rect', g_expr_plt_tbl$grobs[[i]]$grobs[[1]]$childrenOrder))
#     g_expr_plt_tbl$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#     k <- k+1
# }
g_plt_anno <- cowplot::plot_grid(
    plotlist = list(
        g_expr_plt, 
        ggplot() +
            facet_wrap(~Module, ncol = 1) +
            geom_text(
                data = g_plt_dat,
                mapping = aes(
                    x = 0, y = 0, label = str_wrap(Annotation, width = 20)
                ),
                size = 2
            ) +
            coord_cartesian(ylim = c(-0.1, 0.1)) +
            theme_bw() +
            theme(
                text = element_text(size = 20),
                plot.title = element_blank(),
                strip.text = element_blank(),
                strip.background = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                plot.margin = margin(0, unit = "pt")
            )
    ), 
    ncol = 2, align = "h", rel_widths = c(0.75, 0.8)
)
plots <- cowplot::plot_grid(
    plotlist = c(list(g_plt_anno), i_plts), nrow = 1, align = "h"
)
# paginated_plots <- lapply(
#     i_plts,
#     function(i) {
#         cowplot::plot_grid(
#             plotlist = list(g_plt_anno, i), ncol = 2, align = "h"
#         )
#     }
# )
ggsave(file="data/figures/rebuttal_2/ModuleEigengenes.pdf", device = 'pdf', plot = plots, width = 297, height = 210, units = "mm")
ggsave(file="data/figures/rebuttal_2/ModuleEigengenes.png", device = 'png',plot = plots, width = 297, height = 210, units = "mm")
ggsave(file="data/figures/rebuttal_2/ModuleEigengenes.svg", device = 'svg',plot = plots, width = 297, height = 210, units = "mm")

counter=1
for (plt in c(list(g_plt_anno), i_plts)) {
    # pdf(paste0("data/figures/rebuttal_2/ModuleEigengenes_", counter, ".pdf"), paper="a4r", units)
    # print(plt)
    # dev.off()
    ggsave(
        filename = paste0("data/figures/rebuttal_2/ModuleEigengenes_ggsave_", counter, ".pdf"),
        plot = plt,
        width = 297 / 8,
        height = 210,
        units = 'mm'
    )
    counter=counter+1
}
