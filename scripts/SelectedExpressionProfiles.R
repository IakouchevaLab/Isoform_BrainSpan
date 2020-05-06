library(tidyverse)
library(biomaRt)
library(doParallel)
source("scripts/utility/plot_expression.R")

lof <- c(
    "splice_donor_variant" = "Splice Donor",
    "splice_acceptor_variant" = "Splice Acceptor",
    "frameshift_variant" = "Frameshift",
    "stop_gained" = "Stop Gain",
    "start_lost" = "Start Loss"
)
mutations <- read_tsv("data/SNVs/SatterstromProcessedVEP.txt")
lof_mutations <- mutations %>%
    filter(LoF)

anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]
# # Add extra information to annotation table
# ensembl <- useMart(
#     biomart = "ensembl",
#     dataset = "hsapiens_gene_ensembl",
#     host = "GRCh37.ensembl.org"
# )
# anno_expand <- anno %>%
#     left_join(
#         getBM(
#             attributes = c(
#                 "ensembl_gene_id", "gene_biotype",
#                 "start_position", "end_position"
#             ),
#             filters = "ensembl_gene_id",
#             values = unique(.$ensembl_gene_id),
#             mart = ensembl
#         ),
#         by = "ensembl_gene_id"
#     ) %>%
#     mutate(gene_length = abs(start_position - end_position))
# 
# sv_gn <- 13
# sv_tx <- 16
# 
# tt_genes <- readRDS(
#     paste0("data/genes/limma_intermediates/tt_SV", sv_gn, ".rds")
# ) %>%
#     left_join(anno_expand, by = "ensembl_gene_id") %>%
#     dplyr::select(-contains("transcript")) %>%
#     distinct() %>%
#     rename(biotype = gene_biotype) %>%
#     mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
#     mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
#     left_join(
#         lof_mutations %>%
#             dplyr::select(
#                 Chromosome, Variant_start, Variant_end, Ref, Alt,
#                 Affected_status, ensembl_gene_id, ensembl_transcript_id
#             ),
#         by = "ensembl_gene_id"
#     ) %>%
#     mutate(Data = "Gene")
# tt_iso <- readRDS(
#     paste0("data/isoforms/limma_intermediates/tt_SV", sv_tx, ".rds")
# ) %>%
#     left_join(anno_expand, by = "ensembl_transcript_id") %>%
#     rename(biotype = transcript_biotype) %>%
#     mutate(selector = paste0(ensembl_gene_id, Contrast)) %>%
#     mutate(Significant = adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5)) %>%
#     left_join(
#         lof_mutations %>%
#             dplyr::select(
#                 Chromosome, Variant_start, Variant_end, Ref, Alt,
#                 Affected_status, ensembl_gene_id, ensembl_transcript_id
#             ),
#         by = c("ensembl_transcript_id", "ensembl_gene_id")
#     ) %>%
#     mutate(Data = "Isoform")
# 
# de_genes_list <- tt_genes %>%
#     filter(Significant) %>%
#     pull(selector)
# 
# de_iso_list <- tt_iso %>%
#     filter(Significant) %>%
#     pull(selector)
# 
# tt_combined <- bind_rows(
#     tt_genes %>%
#         mutate(Specific = ! selector %in% de_iso_list),
#     tt_iso %>%
#         mutate(Specific = ! selector %in% de_genes_list)
# )
# 
# write_tsv(tt_combined, "data/limmaResultsAll_LoFTargets.tsv")

tt_combined <- read_tsv(
    "data/limmaResultsAll_LoFTargets.tsv",
    col_types = cols(.default = "c")
)

metadata <- read_tsv("data/metadata.tsv")
# iexpr <- readRDS("data/RegressIsoformCounts.rds")
# gexpr <- readRDS("data/RegressGeneCounts.rds")
iexpr <- readRDS("data/iso_tpm_filter.rds")
gexpr <- readRDS("data/gene_tpm_filter.rds")

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

# genes <- c("SCN2A", "DYRK1A", "BTRC", "CELF2", "DLG2")
# genes <- c("BTRC", "DYRK1A")
genes <- c(
    "ADNP", "ANK2", "ANKRD11", "ARID1B", "ASH1L", "ASXL3", "BCL11A", "CACNA2D3",
    "CELF4", "CHD2", "CTNNB1", "DNMT3A", "FOXP1", "GABRB3", "GRIA2", "KDM5B", 
    "KMT2C", "KMT2E", "LDB1", "MBD5", "NCOA1", "NRXN1", "PHF12", "PPP2R5D", 
    "PTK7", "RFX3", "SATB1", "SCN2A", "SETD5", "SHANK2", "SIN3A", "SKI",
    "SMARCC2", "TAOK1", "TBL1XR1", "TBR1", "TCF4", "TCF7L2", "TRIP12", "WAC",
    "ZMYND8", "SCN2A", "DYRK1A", "BTRC", "CELF2", "DLG2"
)
diff_asd_genes <- genes <- c(
    "ADNP", "ANK2", "ANKRD11", "ASXL3", "BCL11A", "CELF4", "CHD2", "CTNNB1", 
    "GABRB3", "GRIA2", "KMT2C", "KMT2E", "MBD5", "PHF12", "PPP2R5D", "PTK7", 
    "RFX3", "SATB1", "SETD5", "SHANK2", "SIN3A", "TBL1XR1", "TCF4", "TRIP12", 
    "WAC", "ZMYND8"
    # , "SCN2A", "DLG2", "CELF2", "DYRK1A", "BTRC"
)

expressions <- lapply(
    setNames(nm = diff_asd_genes),
    function(gene) {
        ensg <- unique(
            anno$ensembl_gene_id[anno$external_gene_id == gene]
        )
        enst <- unique(
            anno$ensembl_transcript_id[anno$external_gene_id == gene]
        )
        # new_rownames <- c(ensg, intersect(rownames(iexpr), enst))
        # expr <- bind_rows(
        #     as.data.frame(gexpr[new_rownames[1], , drop = FALSE]),
        #     as.data.frame(iexpr[new_rownames[-1], , drop = FALSE])
        # ) %>%
        #     as.matrix()
        # rownames(expr) <- new_rownames
        # expr %>%
        iexpr[intersect(enst, rownames(iexpr)), , drop = FALSE] %>%
            expr_to_z() %>%
            as.data.frame() %>%
            mutate(Feature = rownames(.)) %>%
            mutate(
                Data = ifelse(str_detect(Feature, "ENSG"), "Gene", "Isoform")
            ) %>%
            mutate(
                Target = sapply(Feature, function(feature) {
                    if (str_detect(feature, "ENSG")) { return(TRUE) }
                    else {
                        return(
                            feature %in% (
                                lof_mutations %>% 
                                    filter(Affected_status == 2) %>%
                                    pull(ensembl_transcript_id)
                            )
                        )
                    }
                })
            ) %>%
            mutate(
                Class = paste0(Data, Target) %>%
                    str_replace(., "GeneTRUE", "Gene") %>%
                    str_replace(., "IsoformTRUE", "Isoform Target") %>%
                    str_replace(., "IsoformFALSE", "Isoform NonTarget")
            ) %>%
            reshape2::melt(
                id.vars = c(
                    "Feature", "Data", "Target", "Class"
                )
            ) %>%
            rename(
                Sample = variable, 
                `Expression Z-Score` = value
            ) %>%
            left_join(metadata, by = "Sample") %>%
            left_join(
                bind_rows(
                    distinct(data.frame(
                        Feature = anno$ensembl_gene_id,
                        Gene = anno$external_gene_id
                    )),
                    distinct(data.frame(
                        Feature = anno$ensembl_transcript_id,
                        Gene = anno$external_gene_id
                    ))
                ),
                by = "Feature"
            ) %>%
            mutate(external_gene_id = gene)
    }
) %>%
    bind_rows() %>%
    left_join(
        dplyr::select(anno, ensembl_transcript_id, transcript_biotype), 
        by = c("Feature" = "ensembl_transcript_id")
    ) %>%
    filter(transcript_biotype == "protein_coding")

expr_plots <- lapply(
    setNames(nm = sort(genes)),
    function(gene) {
        p <- ggplot() +
            facet_grid(Gene ~ .) +
            geom_line(
                data = expressions %>%
                    filter(external_gene_id == gene) %>%
                    filter(!Target),
                mapping = aes(
                    x = Days, y = `Expression Z-Score`,
                    group = Feature, colour = Class
                ),
                stat = "smooth", method = "loess", size = 2
            ) +
            geom_line(
                data = expressions %>%
                    filter(external_gene_id == gene) %>%
                    filter(Target),
                mapping = aes(
                    x = Days, y = `Expression Z-Score`,
                    group = Feature, colour = Class
                ),
                stat = "smooth", method = "loess", size = 2
            ) +
            labs(
                x = "Days",
                y = "Expression Z-Score"
            ) +
            scale_colour_manual(
                breaks = c(
                    "Gene",
                    "Isoform Target",
                    "Isoform NonTarget"
                ),
                values = c(
                    "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
                    "Isoform Target" = rgb(210 / 255, 53 / 255, 43 / 255),
                    "Isoform NonTarget" = rgb(72 / 255, 126 / 255, 179 / 255)
                ),
                drop = FALSE
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 20)
            )
        expression_plot_fill(p, metadata)
    }
)
cowplot::plot_grid(
    plotlist = lapply(
        expr_plots,
        function(plt) plt + theme(legend.position = "none")
    )
) %>% {
    ggsave(
        filename = "data/figures/26ASD_expression.pdf",
        plot = .,
        device = "pdf", width = 24, height = 16
    )
}
# saveRDS(expr_plots, "data/figures/DiffASDIsoTargets_ExprPlts.rds")
# pdf("data/figures/DiffASDIsoTargets_ExprPlts.pdf", width = 16, height = 9)
# invisible(lapply(expr_plots, print))
# dev.off()
# 
# # Combine a few profiles
# combine_genes <- c("GABRB3", "KMT2C", "MBD5", "PTK7")
# expr_plot_combine <- ggplot() +
#     facet_grid(. ~ Gene) +
#     geom_line(
#         data = expressions %>%
#             filter(external_gene_id %in% combine_genes) %>%
#             filter(!Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 1
#     ) +
#     geom_line(
#         data = expressions %>%
#             filter(external_gene_id %in% combine_genes) %>%
#             filter(Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 1.5
#     ) +
#     labs(
#         x = "Days",
#         y = "Expression Z-Score"
#     ) +
#     scale_colour_manual(
#         breaks = c(
#             "Gene",
#             "Isoform Target",
#             "Isoform NonTarget"
#         ),
#         values = c(
#             "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
#             "Isoform Target" = rgb(210 / 255, 53 / 255, 43 / 255),
#             "Isoform NonTarget" = rgb(72 / 255, 126 / 255, 179 / 255)
#         ),
#         drop = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 20)
#     )
# expr_plot_combine <- expression_plot_fill(expr_plot_combine, metadata)
# saveRDS(expr_plot_combine, "data/figures/DiffASDIsoTargets_ExprPltsCombine.rds")
# 
# # Overlay all expressions
# 
# expressions_overlay <- iexpr[
#     intersect(unique(
#         anno$ensembl_transcript_id[anno$external_gene_id %in% diff_asd_genes]
#     ), rownames(iexpr)), , drop = FALSE
# ] %>%
#     expr_to_z() %>%
#     as.data.frame() %>%
#     mutate(Feature = rownames(.)) %>%
#     mutate(
#         Data = ifelse(str_detect(Feature, "ENSG"), "Gene", "Isoform")
#     ) %>%
#     mutate(
#         Target = sapply(Feature, function(feature) {
#             if (str_detect(feature, "ENSG")) { return(TRUE) }
#             else {
#                 return(
#                     feature %in% (
#                         lof_mutations %>% 
#                             filter(Affected_status == 2) %>%
#                             pull(ensembl_transcript_id)
#                     )
#                 )
#             }
#         })
#     ) %>%
#     mutate(
#         Class = paste0(Data, Target) %>%
#             str_replace(., "GeneTRUE", "Gene") %>%
#             str_replace(., "IsoformTRUE", "Isoform Target") %>%
#             str_replace(., "IsoformFALSE", "Isoform NonTarget")
#     ) %>%
#     reshape2::melt(
#         id.vars = c(
#             "Feature", "Data", "Target", "Class"
#         )
#     ) %>%
#     rename(
#         Sample = variable, 
#         `Expression Z-Score` = value
#     ) %>%
#     left_join(metadata, by = "Sample") %>%
#     left_join(
#         bind_rows(
#             distinct(data.frame(
#                 Feature = anno$ensembl_gene_id,
#                 Gene = anno$external_gene_id
#             )),
#             distinct(data.frame(
#                 Feature = anno$ensembl_transcript_id,
#                 Gene = anno$external_gene_id
#             ))
#         ),
#         by = "Feature"
#     ) %>%
#     # mutate(external_gene_id = gene) %>%
#     left_join(
#         dplyr::select(
#             anno, ensembl_transcript_id, transcript_biotype, external_gene_id
#         ), 
#         by = c("Feature" = "ensembl_transcript_id")
#     ) %>%
#     filter(transcript_biotype == "protein_coding")
# p_overlay <- ggplot() +
#     # facet_grid(Gene ~ .) +
#     geom_line(
#         data = expressions_overlay %>%
#             filter(!Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 2
#     ) +
#     geom_line(
#         data = expressions_overlay %>%
#             filter(Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 2
#     ) +
#     labs(
#         x = "Days",
#         y = "Expression Z-Score"
#     ) +
#     scale_colour_manual(
#         breaks = c(
#             "Gene",
#             "Isoform Target",
#             "Isoform NonTarget"
#         ),
#         values = c(
#             "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
#             "Isoform Target" = rgb(84 / 255, 188 / 255, 194 / 255),
#             "Isoform NonTarget" = "grey"
#         ),
#         drop = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 20)
#     )
# p_overlay <- expression_plot_fill(p_overlay, metadata)
# p_overlay
# 
# 
# # Experimental genes
# 
# exp_genes <- c(
#     "SCN2A", "DLG2", "CELF2", "DYRK1A", "BTRC"
# )
# experimental_expressions <- expressions %>%
#     filter(external_gene_id %in% exp_genes)
# ggplot() +
#     facet_grid(. ~ Gene) +
#     geom_line(
#         data = expressions %>%
#             filter(external_gene_id %in% exp_genes) %>%
#             filter(!Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 1
#     ) +
#     geom_line(
#         data = expressions %>%
#             filter(external_gene_id %in% exp_genes) %>%
#             filter(Target),
#         mapping = aes(
#             x = Days, y = `Expression Z-Score`,
#             group = Feature, colour = Class
#         ),
#         stat = "smooth", method = "loess", size = 1.5
#     ) +
#     labs(
#         x = "Days",
#         y = "Expression Z-Score"
#     ) +
#     scale_colour_manual(
#         breaks = c(
#             "Gene",
#             "Isoform Target",
#             "Isoform NonTarget"
#         ),
#         values = c(
#             "Gene" = rgb(232 / 255, 125 / 255, 114 / 255),
#             "Isoform Target" = rgb(210 / 255, 53 / 255, 43 / 255),
#             "Isoform NonTarget" = rgb(72 / 255, 126 / 255, 179 / 255)
#         ),
#         drop = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         text = element_text(size = 20)
#     )
# expr_plot_combine <- expression_plot_fill(expr_plot_combine, metadata)