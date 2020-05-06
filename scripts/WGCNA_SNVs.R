library(tidyverse)
library(biomaRt)
library(cowplot)
library(doParallel)

################################################################################
# Globals                                                                      #
################################################################################

registerDoParallel(cores = detectCores() - 1)
n_perm <- 1000

lof <- c(
    "splice_donor_variant" = "Splice Donor",
    "splice_acceptor_variant" = "Splice Acceptor",
    "frameshift_variant" = "Frameshift",
    "stop_gained" = "Stop Gain",
    "start_lost" = "Start Loss"
)

#' Combine lengths such that overlapping coordinates are not counted more than 
#' once and gaps are removed
#' 
#' @param feature_starts Coordinate start positions
#' @param feature_ends Coordeinate end positions
#' 
#' @return Total combined length
combine_lengths <- function(feature_starts, feature_ends) {
    # Create data frame for easy sorting/mapping
    feature_coordinates <- arrange(data.frame(
        start = feature_starts,
        end = feature_ends
    ), feature_starts)
    # Iterate through coordinate sets while updating a total length
    cur_start <- NA
    cur_end <- NA
    total_length <- 0
    for (i in 1:nrow(feature_coordinates)) {
        if (i == 1 & is.na(cur_start) & is.na(cur_end)) {
            cur_start <- feature_coordinates[i, "start"]
            cur_end <- feature_coordinates[i, "end"]
            next
        }
        next_start <- feature_coordinates[i, "start"]
        next_end <- feature_coordinates[i, "end"]
        if (next_start <= cur_end) {
            cur_end <- next_end
        } else {
            total_length <- total_length + (cur_end - cur_start) + 1
            cur_start <- next_start
            cur_end <- next_end
        }
    }
    # There's one last chunk to process
    total_length <- total_length + (cur_end - cur_start) + 1
    return(total_length)
}

#' Select neighbor randomly based on GC content and length
#' 
#' @param query Vector of feature queries
#' @param query_type Gene or isoform query
#' @param sel List of possible candidate neighbors, should contain query
#' @param gc Vector of GC content per candidate (must be same order as sel)
#' @param gc_range Range of possible GC content
#' @param len Vector of lengths per candidate (must be same order as sel)
#' @param len_range Range of possible lengths
#' @param n Number of selections to return
#' @param self Should returning query be possible? (FALSE)
#' 
#' @return Random selection(s) that constitutes a neighbor to the query
select_neighbor <- function(
    query, sel, gc, gc_range = 0.1, len, len_range = 0.1, n, self = FALSE
) {
    require(dplyr)
    if (length(sel) != length(len)) {
        stop("len vector is not same length as sel vector")
    }
    if (length(sel) != length(gc)) {
        stop("len vector is not same length as sel vector")
    }
    selection_vector <- foreach (i = 1:length(query)) %dopar% {
        q <- query[i]
        query_idx <- which(sel == q)
        query_gc <- gc[query_idx]
        query_len <- len[query_idx]
        if (!self) {
            sel <- sel[-query_idx]
            gc <- gc[-query_idx]
            len <- len[-query_idx]
        }
        sel_df <- data.frame(
            sel = sel, gc = gc, len = len
        )
        sel_df_filter <- sel_df %>%
            filter(
                gc <= gc * (1 + gc_range) & gc >= gc * (1 - gc_range)
            ) %>%
            filter(
                len <= len * (1 + len_range) & len >= len * (1 - len_range)
            ) %>%
            pull(sel) %>%
            as.character()
        if (length(sel_df_filter) == 0) {selection <- q}
        else {selection <- sample(sel_df_filter, n, replace = TRUE)}
        return(selection)
    }
    return(selection_vector)
}

################################################################################
# Data                                                                         #
################################################################################

asd_genes <- readxl::read_xlsx(
    "data/source/CuratedLists/ASDRelevantGeneListsFromLiterature.xlsx",
    sheet = "SatterstromASD"
)[[1]]

mutations <- read_tsv("data/SNVs/SatterstromProcessedVEP.txt")
lof_mutations <- mutations %>%
    filter(LoF)

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", sep = ",",
    header = TRUE, row.names = 1
)

module_assigns <- bind_rows(
    read_tsv("data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        mutate(ensembl_gene_id = Feature) %>%
        gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
        filter(module_label == as.numeric(gsub("kME", "", kME_Key))) %>%
        mutate(Data = "Gene") %>%
        left_join(
            biomaRt::getBM(
                attributes = c(
                    "ensembl_gene_id", 
                    "start_position", "end_position",
                    "percentage_gene_gc_content"
                ),
                filters = c("ensembl_gene_id"),
                values = .$Feature,
                mart = biomaRt::useMart(
                    biomart = "ensembl",
                    dataset = "hsapiens_gene_ensembl",
                    host = "GRCh37.ensembl.org"
                )
            ),
            by = c("Feature" = "ensembl_gene_id")
        ) %>%
        rename(
            Feature_start = start_position,
            Feature_end = end_position
        ),
    read_tsv("data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
        dplyr::select(
            Feature, module_label, module_colour, starts_with("kME")
        ) %>%
        gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
        filter(module_label == as.numeric(gsub("kME", "", kME_Key))) %>%
        mutate(Data = "Isoform") %>%
        left_join(
            dplyr::select(
                annotations, ensembl_transcript_id, ensembl_gene_id
            ),
            by = c("Feature" = "ensembl_transcript_id")
        ) %>%
        left_join(
            biomaRt::getBM(
                attributes = c(
                    "ensembl_transcript_id", 
                    "transcript_start", "transcript_end",
                    "percentage_gene_gc_content"
                ),
                filters = c("ensembl_transcript_id"),
                values = .$Feature,
                mart = biomaRt::useMart(
                    biomart = "ensembl",
                    dataset = "hsapiens_gene_ensembl",
                    host = "GRCh37.ensembl.org"
                )
            ),
            by = c("Feature" = "ensembl_transcript_id")
        ) %>%
        rename(
            Feature_start = transcript_start,
            Feature_end = transcript_end
        )
) %>%
    mutate(Feature_length = Feature_end - Feature_start)

module_info <- module_assigns %>%
    distinct(Data, module_label, module_colour) %>%
    bind_cols(
        apply(
            ., 1, function(x) {
                dat <- as.character(x[["Data"]])
                ml <- as.numeric(x[["module_label"]])
                mc <- as.character(x[["module_colour"]])
                mm <- filter(
                        module_assigns,
                        Data == dat & module_label == ml & module_colour == mc
                    )
                data.frame(
                    module_length = combine_lengths(
                        mm$Feature_start, mm$Feature_end
                    )
                )
            }
        ) %>%
            bind_rows()
    )
module_assigns <- left_join(
    module_assigns,
    module_info,
    by = c("Data", "module_label", "module_colour")
)
################################################################################
# Analysis                                                                     #
################################################################################

vpm <- module_info %>%
    mutate(
        Affected_status = rep(list(c(1, 2)), nrow(.))
    ) %>%
    unnest(Affected_status) %>%
    mutate(
        Variants = mapply(
            function(dat, ml, mc, mlen, af) {
                join_by = if (dat == "Gene") {
                    c("Feature" = "ensembl_gene_id")
                } else {
                    c("Feature" = "ensembl_transcript_id")
                }
                module_assigns %>%
                    filter(
                        Data == dat & module_label == ml & module_colour == mc
                    ) %>%
                    inner_join(
                        filter(lof_mutations, Affected_status == af),
                        by = join_by
                    ) %>%
                    distinct(`#Uploaded_variation`) %>%
                    nrow()
            },
            Data, module_label, module_colour, module_length, Affected_status
        )
    ) %>%
    mutate(VPK = Variants / (module_length / 1000)) %>%
    mutate(
        ScaleFactor = sapply(Affected_status, function(af) {
            filter(lof_mutations, Affected_status == af) %>%
                distinct(`#Uploaded_variation`) %>%
                nrow()
        }) / 1000000
    ) %>%
    mutate(VPM = VPK / ScaleFactor) %>%
    filter(module_label != 0)
vpm

vpm_plot <- ggplot(
    data = vpm,
    mapping = aes(
        x = module_label, y = VPM, 
        fill = module_colour, colour = Data
    )
) +
    facet_grid(. ~ Affected_status) +
    geom_point(
        shape = 21, size = 3, stroke = 1.1
    ) +
    scale_fill_manual(
        values = setNames(
            nm = vpm$module_colour
        )
    ) +
    guides(fill = FALSE) +
    labs(
        x = "Module", y = "Variants Per Million"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
vpm_plot
saveRDS(vpm_plot, "data/figures/VPM_byStatus.rds")
ggsave(
    filename = "data/figures/VPM_byStatus.pdf",
    plot = vpm_plot,
    device = "pdf", width = 16, height = 12
)

# Do modules have significantly higher VPM than by chance? Permutation test

# message("Begin permutation calculations")
# module_assigns_permute <- module_assigns %>%
#     mutate(
#         random = select_neighbor(
#             query = Feature, sel = Feature, 
#             gc = percentage_gene_gc_content, gc_range = 0.1,
#             len = Feature_length, len_range = 0.1, 
#             n = n_perm, self = TRUE
#         )
#     )
# 
# saveRDS(module_assigns_permute, "data/WGCNA_10000Neighbors.rds")
# module_assigns_permute <- readRDS("data/WGCNA_10000Neighbors.rds")

# Calculate VPM for permuted values
# vpm_permute <- vector(mode = "list", length = n_perm)
# for (n in 1:n_perm) {
#     message("VPM Permutation Calculation: ", n)
#     perm_assign <- module_assigns_permute %>%
#         mutate(
#             Feature = sapply(random, function(r) r[[n]])
#         ) %>%
#         dplyr::select(-random)
#     vpm_permute[[n]] <- module_info %>%
#         mutate(
#             Affected_status = rep(list(c(1, 2)), nrow(.))
#         ) %>%
#         unnest(Affected_status) %>%
#         mutate(
#             Variants = mapply(
#                 function(dat, ml, mc, mlen, af) {
#                     join_by = if (dat == "Gene") {
#                         c("Feature" = "ensembl_gene_id")
#                     } else {
#                         c("Feature" = "ensembl_transcript_id")
#                     }
#                     perm_assign %>%
#                         filter(
#                             Data == dat & 
#                                 module_label == ml & 
#                                 module_colour == mc
#                         ) %>%
#                         inner_join(
#                             filter(lof_mutations, Affected_status == af),
#                             by = join_by
#                         ) %>%
#                         distinct(`#Uploaded_variation`) %>%
#                         nrow()
#                 },
#                 Data, module_label, module_colour, module_length, Affected_status
#             )
#         ) %>%
#         mutate(VPK = Variants / (module_length / 1000)) %>%
#         mutate(
#             ScaleFactor = sapply(Affected_status, function(af) {
#                 filter(lof_mutations, Affected_status == af) %>%
#                     distinct(`#Uploaded_variation`) %>%
#                     nrow()
#             }) / 1000000
#         ) %>%
#         mutate(VPM = VPK / ScaleFactor) %>%
#         filter(module_label != 0)
# }
# saveRDS(vpm_permute, "data/VPM_Permutations.rds")
vpm_permute <- readRDS("data/VPM_Permutations.rds")

# Calculate distribution of VPM differences
null_vpm_difference <- data.frame()
for (i in 1:length(vpm_permute)) {
    null_vpm_difference <- bind_rows(
        null_vpm_difference,
        vpm_permute[[i]] %>%
            dplyr::select(
                Data, module_label, module_colour, module_length, 
                Affected_status, VPM
            ) %>%
            spread(key = Affected_status, value = VPM, ) %>%
            rename(Affected_status_1 = `1`, Affected_status_2 = `2`) %>%
            mutate(VPM_Difference = Affected_status_2 - Affected_status_1) %>%
            mutate(Permutation = i)
    )
        
}

write_tsv(null_vpm_difference, "data/SNVs/nullVPM.tsv")

empirical_vpm_difference <- vpm %>%
    dplyr::select(
        Data, module_label, module_colour, module_length, 
        Affected_status, VPM
    ) %>%
    spread(key = Affected_status, value = VPM, ) %>%
    rename(Affected_status_1 = `1`, Affected_status_2 = `2`) %>%
    mutate(VPM_Difference = Affected_status_2 - Affected_status_1) %>%
    mutate(
        p_value_gt = mapply(
            function(dat, ml, mc, vpm_diff) {
                null_diff <- null_vpm_difference %>%
                    filter(
                        Data == dat &
                            module_label == ml & 
                            module_colour == mc
                    ) %>%
                    pull(VPM_Difference)
                gt_prop = sum(null_diff >= vpm_diff) / (length(null_diff) + 1)
            },
            Data, module_label, module_colour, VPM_Difference
        )
    ) %>%
    gather(
        Affected_status_1, Affected_status_2, 
        key = "Affected_status", value = "VPM"
    ) %>%
    mutate(
        Affected_status = ifelse(
            Affected_status == "Affected_status_1", 
            "Control", "Case"
        )
    ) %>%
    mutate(adj.P.Val = p.adjust(p_value_gt, method = "BH")) %>%
    mutate(Star = ifelse(adj.P.Val <= 0.05, "*", "")) %>%
    group_by(Data, module_label) %>%
    mutate(y_star = max(VPM) * 1.05) %>%
    mutate(y_bar = max(VPM) * 1.01) %>%
    mutate(x_min = module_label - 0.25, x_max = module_label + 0.25)

write_tsv(empirical_vpm_difference, "data/SNVs/empiricalVPM.tsv")

vpm_diff <- ggplot(
    data = empirical_vpm_difference,
    mapping = aes(
        x = module_label, y = VPM, 
        fill = module_colour, 
        colour = Affected_status
    )
) +
    facet_wrap(~ Data, scales = "free_x") +
    labs(
        x = "Module", y = "Variants per Million"
    ) +
    geom_segment(
        mapping = aes(
            x = module_label, xend = module_label, y = 0, yend = VPM
        ),
        colour = "black"
    ) +
    geom_point(
        mapping = aes(
            x = module_label, y = VPM,
            shape = Affected_status
        ),
        size = 5, stroke = 2
    ) +
    # geom_bar(
    #     mapping = aes(
    #         size = Data
    #     ),
    #     stat = "identity", position = "dodge"
    # ) +
    geom_text(
        mapping = aes(
            label = Star, y = y_star
        ), colour = "black", size = 10
    ) +
    # geom_segment(
    #     data = filter(empirical_vpm_difference, adj.P.Val <= 0.05),
    #     mapping = aes(
    #         x = x_min, xend = x_max, y = y_bar, yend = y_bar
    #     ), colour = "black"
    # ) +
    scale_colour_manual(
        values = c(
            "Case" = rgb(210 / 255, 53 / 255, 43 / 255),
            "Control" = rgb(72 / 255, 126 / 255, 179 / 255)
        )
    ) +
    scale_fill_manual(
        values = setNames(
            nm = empirical_vpm_difference$module_colour
        )
    ) +
    scale_size_manual(
        values = c(
            "Gene" = 1.5,
            "Isoform" = 0.75
        )
    ) +
    scale_shape_manual(
        values = c(
            "Case" = 23,
            "Control" = 22
        )
    ) +
    guides(
        colour = guide_legend(override.aes = list(size = 2)),
        fill = FALSE,
        size = FALSE
    ) + 
    theme_bw() +
    theme(
        text = element_text(size = 30),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.justification = c(1, 1),
        legend.position = c(0.5, 0.9),
        panel.spacing = unit(0, "lines")
    )
vpm_diff
ggsave(
    filename = "data/figures/VPM_byStatus_lollipop.pdf",
    plot = vpm_diff,
    device = "pdf", width = 16, height = 12
)
saveRDS(vpm_diff, "data/figures/VPM_byStatus_lollipop.rds")

# gM1, iM1, iM30 are all significantly impacted by case LoF mutations vs ctrl

m1_gost <- bind_rows(
    readRDS("data/WGCNA/GeneNetwork_DS2_MM20_GOST.rds")[["1 turquoise"]] %>%
        mutate(Data = "Gene"), 
    readRDS("data/WGCNA/IsoformNetwork_DS2_MM20_GOST.rds")[["1 turquoise"]] %>%
        mutate(Data = "Isoform")
) %>%
    filter(term_size <= 1000) %>%
    mutate(Order = seq(nrow(.), 1)) %>%
    mutate(Data = factor(Data, levels = c("Gene", "Isoform"))) %>%
    group_by(Data) %>%
    top_n(10, Order) %>%
    ggplot(
        data = .,
        mapping = aes(
            x = -log10(p_value), 
            y = reorder(str_wrap(term_name, 30), Order),
            size = intersection_size, shape = Data
        )
    ) +
    labs(
        x = expression("-log"["10"]*"FDR"),
        title = "g/iM1 turquoise"
    ) +
    geom_point(
        colour = "black", fill = "turquoise"
    ) +
    scale_size_continuous(range = c(5, 15)) +
    scale_shape_manual(
        values = c("Gene" = 21, "Isoform" = 22)
    ) +
    guides(
        fill = FALSE,
        size = guide_legend(title = "Feature Intersection"),
        shape = guide_legend(override.aes = list(size = 10, fill = NA))
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.title.y = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.95, 0.05),
        legend.background = element_rect(
            fill = "white", colour = "black"
        )
    )
m1_gost
saveRDS(m1_gost, "data/figures/GIM1_ENR.rds")

im30_gost <- readRDS(
    "data/WGCNA/IsoformNetwork_DS2_MM20_GOST.rds"
)[["30 steelblue"]] %>%
        mutate(Data = "Isoform") %>%
    filter(term_size <= 1000) %>%
    mutate(Order = seq(nrow(.), 1)) %>%
    mutate(Data = factor(Data, levels = c("Gene", "Isoform"))) %>%
    group_by(Data) %>%
    top_n(10, Order) %>%
    ggplot(
        data = .,
        mapping = aes(
            x = -log10(p_value), 
            y = reorder(str_wrap(term_name, 30), Order),
            size = intersection_size, shape = Data
        )
    ) +
    labs(
        x = expression("-log"["10"]*"FDR"),
        title = "iM30 steelblue"
    ) +
    geom_point(
        colour = "black", fill = "steelblue"
    ) +
    scale_size_continuous(range = c(5, 15)) +
    scale_shape_manual(
        values = c("Gene" = 21, "Isoform" = 22)
    ) +
    guides(
        fill = FALSE,
        size = guide_legend(title = "Feature Intersection"),
        shape = guide_legend(override.aes = list(size = 10))
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30),
        axis.title.y = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.95, 0.05),
        legend.background = element_rect(
            fill = "white", colour = "black"
        )
    )
im30_gost
saveRDS(im30_gost, "data/figures/IM30_ENR.rds")

sigModules <- data.frame(
    Data = c("Gene", "Isoform", "Isoform"),
    module_label = c(1, 1, 30),
    module_colour = c("turquoise", "turquoise", "steelblue")
)

sigModules_funEnrich <- apply(
    sigModules,
    1,
    function(x) {
        dat <- as.character(x[["Data"]])
        ml <- as.numeric(x[["module_label"]])
        mc <- as.character(x[["module_colour"]])
        message(dat, ml, mc)
        ma <- module_assigns %>%
            filter(Data == dat & module_label == ml & module_colour == mc) %>%
            filter(kME_Key == paste0("kME", trimws(ml))) %>%
            arrange(desc(kME))
        bg <- if (dat == "Gene") {
            module_assigns$Feature[module_assigns$Data == "Gene"]
        } else {
            module_assigns$Feature[module_assigns$Data == "Isoform"]
        }
        res <- gprofiler2::gost(
            query = as.character(pull(ma, Feature)),
            organism = "hsapiens", 
            ordered_query = TRUE,
            correction_method = "fdr",
            custom_bg = bg,
            source = c("GO:BP", "GO:MF")
        )$result
        if (is.null(res)) {
            return(res)
        } else {
            res <- res %>%
                mutate(Data = dat, module_label = ml, module_colour = mc)
        }
    }
)

sigModules_funEnrich_sepPlots <- lapply(
    sigModules_funEnrich,
    function(enr) {
        enr_sort <- enr %>%
            filter(term_size <= 1000) %>%
            mutate(Order = seq(nrow(.), 1)) %>%
            top_n(10, Order) %>%
            mutate(Data = factor(Data, levels = c("Gene", "Isoform")))
        ggplot(
            data = enr_sort,
            mapping = aes(
                x = -log10(p_value), 
                y = reorder(str_wrap(term_name, 30), Order),
                fill = module_colour, size = intersection_size
            )
        ) +
            labs(
                x = expression("-log"["10"]*"FDR"),
                title = paste(
                    enr_sort$Data[[1]],
                    enr_sort$module_label[[1]],
                    enr_sort$module_colour[[1]]
                )
            ) +
            geom_point(
                shape = 21, colour = "black"
            ) +
            scale_fill_manual(
                values = setNames(nm = module_assigns$module_colour)
            ) +
            scale_size_continuous(range = c(5, 15)) +
            guides(
                fill = FALSE,
                size = guide_legend(title = "Feature Intersection")
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 30),
                axis.title.y = element_blank(),
                legend.justification = c(1, 0),
                legend.position = c(0.95, 0.05),
                legend.background = element_rect(
                    fill = "white", colour = "black"
                )
            )
    }
)
sigModules_funEnrich_sepPlots[[3]]

saveRDS(
    sigModules_funEnrich_sepPlots,
    "data/figures/SignificantVPM_ModuleGOST_sepPlots.rds"
)

# Which genes/isoforms are actually impacted in the significant cases?
sigModules_targets <- apply(
    sigModules, 1, function(x) {
        dat <- as.character(x[["Data"]])
        ml <- as.numeric(x[["module_label"]])
        mc <- as.character(x[["module_colour"]])
        message(dat, ml, mc)
        merge_col <- if (dat == "Gene") {
            "ensembl_gene_id"
        } else {
            "ensembl_transcript_id"
        }
        module_assigns_lof <- module_assigns %>%
            filter(Data == dat & module_label == ml & module_colour == mc) %>%
            filter(kME_Key == paste0("kME", trimws(ml))) %>%
            arrange(desc(kME)) %>%
            inner_join(
                lof_mutations,
                by = c("Feature" = merge_col)
            ) %>%
            distinct(
                Data, module_label, module_colour, 
                Feature, `#Uploaded_variation`, Affected_status
            )
        external_ids <- getBM(
            attributes = c(
                "external_gene_name", "external_transcript_name", 
                "ensembl_gene_id", "ensembl_transcript_id"
            ),
            filters = merge_col,
            values = unique(as.character(module_assigns_lof$Feature)),
            mart = useMart(
                biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl",
                host = "GRCh37.ensembl.org"
            )
        )
        module_assigns_lof <- module_assigns_lof %>%
            left_join(
                external_ids,
                by = c("Feature" = merge_col)
            )
        if (dat == "Gene") {
            module_assigns_lof %>%
                mutate(ensembl_transcript_id = NA, external_transcript_id = NA)
        }
        return(module_assigns_lof)
    }
) %>%
    bind_rows() %>%
    dplyr::select(-contains("ensembl")) %>%
    mutate(
        external_feature_id = ifelse(
            Data == "Gene", external_gene_name, external_transcript_name
        )
    ) %>%
    dplyr::select(
        -external_transcript_name
    ) %>%
    distinct() %>%
    mutate(ASD = external_gene_name %in% asd_genes)

sigModules_targets_counts <- sigModules_targets %>%
    count(
        Data, module_label, module_colour, Feature, external_feature_id,
        ASD, Affected_status
    ) %>%
    rename(NumberVariants = n)
View(sigModules_targets_counts)

write_tsv(sigModules_targets, "data/SNVs/WGCNA_SigVPMModule_Targets.tsv")

# Extra

sigModules_extra <- data.frame(
    Data = c("Isoform", "Isoform", "Isoform"),
    module_label = c(30, 17, 27),
    module_colour = c("steelblue", "grey60", "white")
)

sigModules_funEnrich_extra <- apply(
    sigModules_extra,
    1,
    function(x) {
        dat <- as.character(x[["Data"]])
        ml <- as.numeric(x[["module_label"]])
        mc <- as.character(x[["module_colour"]])
        message(dat, ml, mc)
        ma <- module_assigns %>%
            filter(Data == dat & module_label == ml & module_colour == mc) %>%
            filter(kME_Key == paste0("kME", trimws(ml))) %>%
            arrange(desc(kME))
        bg <- if (dat == "Gene") {
            module_assigns$Feature[module_assigns$Data == "Gene"]
        } else {
            module_assigns$Feature[module_assigns$Data == "Isoform"]
        }
        res <- gprofiler2::gost(
            query = as.character(pull(ma, Feature)),
            organism = "hsapiens", 
            ordered_query = TRUE,
            correction_method = "fdr",
            custom_bg = bg,
            source = c("GO:BP", "GO:MF")
        )$result
        if (is.null(res)) {
            return(res)
        } else {
            res <- res %>%
                mutate(Data = dat, module_label = ml, module_colour = mc)
        }
    }
)

sigModules_funEnrich_sepPlots_extra <- lapply(
    sigModules_funEnrich_extra,
    function(enr) {
        enr_sort <- enr %>%
            filter(term_size <= 1000) %>%
            mutate(Order = seq(nrow(.), 1)) %>%
            top_n(10, Order) %>%
            mutate(Data = factor(Data, levels = c("Gene", "Isoform")))
        ggplot(
            data = enr_sort,
            mapping = aes(
                x = -log10(p_value), 
                y = reorder(str_wrap(term_name, 30), Order),
                fill = module_colour, size = intersection_size
            )
        ) +
            labs(
                x = expression("-log"["10"]*"FDR"),
                title = paste(
                    enr_sort$Data[[1]],
                    enr_sort$module_label[[1]],
                    enr_sort$module_colour[[1]]
                )
            ) +
            geom_point(
                shape = 21, colour = "black"
            ) +
            scale_fill_manual(
                values = setNames(nm = module_assigns$module_colour)
            ) +
            scale_size_continuous(range = c(5, 15)) +
            guides(
                fill = FALSE,
                size = guide_legend(title = "Feature Intersection")
            ) +
            theme_bw() +
            theme(
                text = element_text(size = 30),
                axis.title.y = element_blank(),
                legend.justification = c(1, 0),
                legend.position = c(0.95, 0.05),
                legend.background = element_rect(
                    fill = "white", colour = "black"
                )
            )
    }
)
sigModule_funEnrich_combine_extra <- plot_grid(
    plotlist = sigModules_funEnrich_sepPlots_extra %>%
        lapply(
            ., function(p) {
                p + 
                    guides(
                        size = guide_legend(title = "Feature\nIntersection")
                    ) +
                    theme(
                        text = element_text(size = 18),
                        plot.title = element_text(size = 18),
                        axis.text.y = element_text(size = 14),
                        legend.title = element_text(size = 18),
                        legend.text = element_text(size = 14)
                    )
            }
        ),
    nrow = 1
)
sigModule_funEnrich_combine_extra
