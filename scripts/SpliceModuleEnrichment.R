library(tidyverse)
library(biomaRt)
library(cowplot)
library(doParallel)

################################################################################
# Globals                                                                      #
################################################################################

registerDoParallel(cores = detectCores() - 1)
n_perm <- 1000

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

annotations <- read.table(
    "data/source/annotation.transcript.ensg75.txt", sep = ",",
    header = TRUE, row.names = 1
)

# module_assigns <- bind_rows(
#     read_tsv("data/genes/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
#         mutate(ensembl_gene_id = Feature) %>%
#         gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
#         filter(module_label == as.numeric(gsub("kME", "", kME_Key))) %>%
#         mutate(Data = "Gene") %>%
#         left_join(
#             biomaRt::getBM(
#                 attributes = c(
#                     "ensembl_gene_id",
#                     "start_position", "end_position",
#                     "percentage_gene_gc_content"
#                 ),
#                 filters = c("ensembl_gene_id"),
#                 values = .$Feature,
#                 mart = biomaRt::useMart(
#                     biomart = "ensembl",
#                     dataset = "hsapiens_gene_ensembl",
#                     host = "GRCh37.ensembl.org"
#                 )
#             ),
#             by = c("Feature" = "ensembl_gene_id")
#         ) %>%
#         rename(
#             Feature_start = start_position,
#             Feature_end = end_position
#         ),
#     read_tsv("data/isoforms/Networks/Network_DS2_MM20_ModuleAssign.tsv") %>%
#         dplyr::select(
#             Feature, module_label, module_colour, starts_with("kME")
#         ) %>%
#         gather(starts_with("kME"), key = "kME_Key", value = "kME") %>%
#         filter(module_label == as.numeric(gsub("kME", "", kME_Key))) %>%
#         mutate(Data = "Isoform") %>%
#         left_join(
#             dplyr::select(
#                 annotations, ensembl_transcript_id, ensembl_gene_id
#             ),
#             by = c("Feature" = "ensembl_transcript_id")
#         ) %>%
#         left_join(
#             biomaRt::getBM(
#                 attributes = c(
#                     "ensembl_transcript_id",
#                     "transcript_start", "transcript_end",
#                     "percentage_gene_gc_content"
#                 ),
#                 filters = c("ensembl_transcript_id"),
#                 values = .$Feature,
#                 mart = biomaRt::useMart(
#                     biomart = "ensembl",
#                     dataset = "hsapiens_gene_ensembl",
#                     host = "GRCh37.ensembl.org"
#                 )
#             ),
#             by = c("Feature" = "ensembl_transcript_id")
#         ) %>%
#         rename(
#             Feature_start = transcript_start,
#             Feature_end = transcript_end
#         )
# ) %>%
#     mutate(Feature_length = Feature_end - Feature_start)
# 
# module_info <- module_assigns %>%
#     distinct(Data, module_label, module_colour) %>%
#     bind_cols(
#         apply(
#             ., 1, function(x) {
#                 dat <- as.character(x[["Data"]])
#                 ml <- as.numeric(x[["module_label"]])
#                 mc <- as.character(x[["module_colour"]])
#                 mm <- filter(
#                     module_assigns,
#                     Data == dat & module_label == ml & module_colour == mc
#                 )
#                 data.frame(
#                     module_length = combine_lengths(
#                         mm$Feature_start, mm$Feature_end
#                     )
#                 )
#             }
#         ) %>%
#             bind_rows()
#     )
# module_assigns <- left_join(
#     module_assigns,
#     module_info,
#     by = c("Data", "module_label", "module_colour")
# )
# write_tsv(module_assigns, "data/ModuleInfo.tsv")
module_assigns <- read_tsv("data/ModuleInfo.tsv")

transcript_coordinates <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_exon_id",
                   "chromosome_name", "exon_chrom_start", "exon_chrom_end"),
    filters = c("ensembl_transcript_id"),
    values = module_assigns %>% filter(Data == "Isoform") %>% pull(Feature),
    mart = useMart(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        host = "GRCh37.ensembl.org"
    )
) %>%
    as_tibble()
splice_variants <- readxl::read_xlsx(
    "data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = 2
) %>%
    filter(str_detect(Consequence, "splice")) %>%
    filter(!str_detect(Consequence, "region"))

write_tsv(splice_variants, "data/splice_variant_impact.tsv")
splice_variants <- read_tsv("data/splice_variant_impact.tsv")

################################################################################
# Analysis                                                                     #
################################################################################

splice_affected_transcripts <- splice_variants %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Transcript ID`) %>%
    unique()
splice_unaffected_transcripts <- splice_variants %>%
    filter(`Affected status` == 1) %>%
    pull(`Ensembl Transcript ID`)
splice_affected_genes <- splice_variants %>%
    filter(`Affected status` == 2) %>%
    pull(`Ensembl Gene ID`) %>%
    unique()
splice_unaffected_genes <- splice_variants %>%
    filter(`Affected status` == 1) %>%
    pull(`Ensembl Gene ID`)

vpm <- module_assigns %>%
    dplyr::select(module_label, module_colour, Data, module_length) %>%
    distinct() %>%
    mutate(
        Affected_status = rep(list(c(1, 2)), nrow(.))
    ) %>%
    unnest(Affected_status) %>%
    mutate(Variants = mapply(
        function(dat, ml, mc, mlen, af) {
            if (af == 1) {
                gene_vec = splice_unaffected_genes
                iso_vec = splice_unaffected_transcripts
            } else {
                gene_vec = splice_affected_genes
                iso_vec = splice_affected_transcripts
            }
            if (dat == "Gene") {
                features = gene_vec
            } else {
                features = iso_vec
            }
            module_assigns %>%
                filter(module_label == ml) %>%
                filter(Data == dat) %>%
                filter(Feature %in% features) %>%
                nrow()
        },
        Data, module_label, module_colour, module_length, Affected_status
    )) %>%
    mutate(VPK = Variants / (module_length / 1000)) %>%
    mutate(
        ScaleFactor = sapply(Affected_status, function(af) {
            filter(splice_variants, `Affected status` == af) %>%
                distinct(Chromosome, `Variant start`, Reference, Alternate) %>%
                nrow()
        }) / 1000000
    ) %>%
    mutate(VPM = VPK / ScaleFactor) %>%
    filter(module_label != 0)

vpm_diff <- apply(distinct(dplyr::select(vpm, module_label, module_colour, Data)), 1, function(x) {
    vpm_case <- pull(filter(vpm, module_label == x[["module_label"]], module_colour == x[["module_colour"]], Data == x[["Data"]], Affected_status == 2), VPM)
    vpm_ctrl <- pull(filter(vpm, module_label == x[["module_label"]], module_colour == x[["module_colour"]], Data == x[["Data"]], Affected_status == 1), VPM)
    tibble(
        module_label = x[["module_label"]],
        module_colour = x[["module_colour"]],
        Data = x[["Data"]],
        VPM_diff = vpm_case - vpm_ctrl
    )
}) %>%
    bind_rows()

splice_vpm_plt <- ggplot(
    data = vpm,
    mapping = aes(
        x = module_label, y = VPM, fill = module_colour, 
        shape = as.character(Affected_status)
    )
) +
    facet_grid(. ~ Data, scales = "free_x", space = "free_x") +
    geom_segment(
        mapping = aes(
            x = module_label, xend = module_label, y = VPM, yend = 0
        ),
        size = 0.2
    ) +
    geom_point(size = 1.5, stroke = 0.1) +
    scale_shape_manual(
        values = c("1" = 21, "2" = 22),
        labels = c("1" = "ASD", "2" = "Control")
    ) +
    scale_fill_manual(
        values = setNames(
            nm = vpm$module_colour
        )
    ) +
    scale_x_continuous(
        breaks = seq(2, max(vpm$module_label), by = 2),
        expand = expand_scale(mult = 0.025)
    ) +
    guides(fill = FALSE) +
    labs(
        x = "Module", y = "Variants Per Million"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 8),
        axis.text.x = element_text(size = 5),
        legend.title = element_blank(),
        line = element_line(size = 0.1),
        panel.spacing = unit(x = 0, units = "cm"),
        legend.position = "top",
        panel.border = element_rect(size = 0.2, fill = NA),
        strip.background = element_rect(size = 0.2),
        legend.margin = margin(),
        legend.box.spacing = unit(0, units = "cm")
    )
ggsave("data/figures/SpliceVPM.pdf", splice_vpm_plt, device = "pdf", width = 8.5, height = 6, units = "cm", useDingbats = FALSE)

################################################################################
# Significance Testing                                                         #
################################################################################
vpm_permute <- list()
for (random_neighbor_set in list.files("data/", pattern = "module_assigns_Permutations*", full.names = TRUE)[1:1]) {
    print(random_neighbor_set)
    random_module_assigns_total <- readRDS(random_neighbor_set)
    number_permutations <- length(unlist(pull(random_module_assigns_total[1, ], random)))
    for (i in 1:number_permutations) {
        print(i)
        random_module_assigns <- random_module_assigns_total %>%
            mutate(Feature = sapply(random, function(r) r[[i]])) %>%
            dplyr::select(-random)
        this_vpm_permute <- random_module_assigns %>%
            dplyr::select(module_label, module_colour, Data, module_length) %>%
            distinct() %>%
            mutate(
                Affected_status = rep(list(c(1, 2)), nrow(.))
            ) %>%
            unnest(Affected_status) %>%
            mutate(Variants = mapply(
                function(dat, ml, mc, mlen, af) {
                    if (af == 1) {
                        gene_vec = splice_unaffected_genes
                        iso_vec = splice_unaffected_transcripts
                    } else {
                        gene_vec = splice_affected_genes
                        iso_vec = splice_affected_transcripts
                    }
                    if (dat == "Gene") {
                        features = gene_vec
                    } else {
                        features = iso_vec
                    }
                    random_module_assigns %>%
                        filter(module_label == ml) %>%
                        filter(Data == dat) %>%
                        filter(Feature %in% features) %>%
                        nrow()
                },
                Data, module_label, module_colour, module_length, Affected_status
            )) %>%
            mutate(VPK = Variants / (module_length / 1000)) %>%
            mutate(
                ScaleFactor = sapply(Affected_status, function(af) {
                    filter(splice_variants, `Affected status` == af) %>%
                        distinct(Chromosome, `Variant start`, Reference, Alternate) %>%
                        nrow()
                }) / 1000000
            ) %>%
            mutate(VPM = VPK / ScaleFactor) %>%
            filter(module_label != 0)
        vpm_diff_permute <- apply(distinct(dplyr::select(this_vpm_permute, module_label, module_colour, Data)), 1, function(x) {
            vpm_case <- pull(filter(this_vpm_permute, module_label == x[["module_label"]], module_colour == x[["module_colour"]], Data == x[["Data"]], Affected_status == 2), VPM)
            vpm_ctrl <- pull(filter(this_vpm_permute, module_label == x[["module_label"]], module_colour == x[["module_colour"]], Data == x[["Data"]], Affected_status == 1), VPM)
            tibble(
                module_label = x[["module_label"]],
                module_colour = x[["module_colour"]],
                Data = x[["Data"]],
                VPM_diff = vpm_case - vpm_ctrl
            )
        }) %>%
            bind_rows()
        vpm_permute[[length(vpm_permute)+1]] <- vpm_diff_permute
    }
}

saveRDS(vpm_permute, "data/SpliceVPMPermutations.rds")

pval <- apply(vpm_diff, 1, function(v) {
    my_module_label <- as.numeric(v[["module_label"]])
    my_module_colour <- as.character(v[["module_colour"]])
    my_module_data <- as.character(v[["Data"]])
    my_vpm <- as.numeric(v[["VPM_diff"]])
    all_vpm <- sapply(vpm_permute, function(vp) {
        vp %>% 
            filter(module_label == my_module_label) %>%
            filter(module_colour == my_module_colour) %>%
            filter(Data == my_module_data) %>%
            pull(VPM_diff)
    })
    1 - (sum(my_vpm >= all_vpm) / (length(all_vpm) + 1))
})
vpm_w_pval <- vpm_diff %>%
    mutate(pvalue = pval) %>%
    mutate(adj_pval = p.adjust(pvalue, method = "BH"))
write_tsv(vpm_w_pval, "data/SpliceVPM.tsv")
