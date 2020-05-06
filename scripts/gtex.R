# Check brainspan against gtex
library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores() - 1)

# gtex_samples <- readLines(
#     "data/source/GTEx_V8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct",
#     n = 3
# )[[3]] %>%
#     str_split("\t", simplify = TRUE)

phenotypes <- read_tsv(
    "data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
)
sampletypes <- read_tsv(
    "data/source/GTEx_V8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    col_types = cols(.default = "c")
)

sample_select <- tibble(
    SAMPID = sampletypes %>%
        filter(str_detect(SMTS, "Brain")) %>%
        pull(SAMPID)
) %>%
    mutate(
        SUBJID = sapply(
            str_split(SAMPID, "-"),
            function(x) {
                paste(x[[1]], x[[2]], sep = "-")
            }
        )
    ) %>%
    filter(
        SUBJID %in% (
            phenotypes %>%
                filter(AGE %in% c("20-29", "30-39")) %>%
                pull(SUBJID)
        )
    )

# gtex_brain_samples <- paste(c(1, 2, which(gtex_samples %in% sample_select$SAMPID)), collapse = ",")
# sample_cut_command <- paste0(
#     "tail -n +2 GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | cut -f",
#     gtex_brain_samples
# )

gtex_expr <- read_tsv("data/source/GTEx_V8/gtex_brainspan_20-40.txt", skip = 1) %>%
    mutate(ensembl_transcript_id = gsub("\\.[0-9]+", "", transcript_id))
itpm <- readRDS("data/iso_tpm_filter.rds")

metadata <- read_tsv("data/metadata.tsv")
p13samples <- metadata$Sample[metadata$Period == 13]

intersect_isoforms <- intersect(rownames(itpm), gtex_expr$ensembl_transcript_id)

sort(unique(sampletypes$SMTSD))

gtex_brain_samples <- c(
    "AMY" = "Brain - Amygdala",
    "CBC" = "Brain - Cerebellar Hemisphere",
    "HIP" = "Brain - Hippocampus"
)

amy_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "AMY") %>%
    pull(Sample)
cbc_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "CBC") %>%
    pull(Sample)
hip_bspan_samples <- metadata %>%
    filter(Period == 13) %>%
    filter(Regioncode == "HIP") %>%
    pull(Sample)
amy_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Amygdala") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))
cbc_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Cerebellar Hemisphere") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))
hip_gtex_samples <- sampletypes %>%
    filter(SMTSD == "Brain - Hippocampus") %>%
    pull(SAMPID) %>%
    intersect(sample_select$SAMPID) %>%
    intersect(colnames(gtex_expr))

data_cor <- sapply(
    setNames(nm = intersect_isoforms),
    function(tx) {
        amy_bspan <- mean(as.numeric(itpm[tx, amy_bspan_samples]))
        cbc_bspan <- mean(as.numeric(itpm[tx, cbc_bspan_samples]))
        hip_bspan <- mean(as.numeric(itpm[tx, hip_bspan_samples]))
        amy_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, amy_gtex_samples]))
        cbc_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, cbc_gtex_samples]))
        hip_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, hip_gtex_samples]))
        if (all(c(amy_bspan, cbc_bspan, hip_bspan, amy_gtex, cbc_gtex, hip_gtex) == 0)) {
            return(NULL)
        }
        cor(
            c(amy_bspan, cbc_bspan, hip_bspan),
            c(amy_gtex, cbc_gtex, hip_gtex),
        )
    }
)

data_cor_full <- sapply(
    setNames(nm = intersect_isoforms),
    function(tx) {
        amy_bspan <- as.numeric(itpm[tx, amy_bspan_samples])
        cbc_bspan <- as.numeric(itpm[tx, cbc_bspan_samples])
        hip_bspan <- as.numeric(itpm[tx, hip_bspan_samples])
        amy_gtex <- as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, amy_gtex_samples])
        cbc_gtex <- as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, cbc_gtex_samples])
        hip_gtex <- as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx, hip_gtex_samples])
        if (all(c(amy_bspan, cbc_bspan, hip_bspan, amy_gtex, cbc_gtex, hip_gtex) == 0)) {
            return(NULL)
        }
        cor(
            c(amy_bspan, cbc_bspan, hip_bspan),
            c(amy_gtex, cbc_gtex, hip_gtex),
        )
    }
)


# data_core_permutations <- mclapply(
#     1:100,
#     function(x) {
#         sapply(
#             intersect_isoforms,
#             function(tx_bspan) {
#                 tx_gtex <- sample(intersect_isoforms[intersect_isoforms != tx_bspan], 1)
#                 amy_bspan <- mean(as.numeric(itpm[tx_bspan, amy_bspan_samples]))
#                 cbc_bspan <- mean(as.numeric(itpm[tx_bspan, cbc_bspan_samples]))
#                 hip_bspan <- mean(as.numeric(itpm[tx_bspan, hip_bspan_samples]))
#                 amy_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx_gtex, amy_gtex_samples]))
#                 cbc_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx_gtex, cbc_gtex_samples]))
#                 hip_gtex <- mean(as.numeric(gtex_expr[gtex_expr$ensembl_transcript_id == tx_gtex, hip_gtex_samples]))
#                 if (all(c(amy_bspan, cbc_bspan, hip_bspan, amy_gtex, cbc_gtex, hip_gtex) == 0)) {
#                     return(NULL)
#                 }
#                 cor(
#                     c(amy_bspan, cbc_bspan, hip_bspan),
#                     c(amy_gtex, cbc_gtex, hip_gtex),
#                 )
#             }
#         )
#     }
# )
# saveRDS(data_core_permutations, "data/gtex_permutations.rds")
data_core_permutations <- readRDS("data/gtex_permutations.rds")


ggplot() +
    geom_density(
        data = tibble(x = unlist(data_core_permutations)), mapping = aes(x = x), alpha = 0.1, fill = "grey"
    ) +
    geom_density(
        data = tibble(x = unlist(data_cor)), mapping = aes(x = x), binwidth = 0.01, fill = "red", alpha = 0.5
    )

mean(unlist(data_cor), na.rm = TRUE)

#perm_mean <- lapply(data_core_permutations, function(x) mean(unlist(x)), na.rm = TRUE)

data_cor_obs <- unlist(data_cor)
data_cor_perm <- lapply(data_core_permutations, unlist)


1 - (sum(mean(unlist(data_cor), na.rm = TRUE) > lapply(data_core_permutations, function(x) mean(unlist(x), na.rm = TRUE))) / 101)

cor_bins <- tibble(
    bin_lo = seq(-1, 0.9, by = 0.1),
    bin_hi = seq(-0.9, 1, by = 0.1)
) %>%
    mutate(label = paste(bin_lo, bin_hi, sep = " ~ ")) %>%
    mutate(
        n = mapply(
            function(lo, hi) {
                cors <- unlist(data_cor)
                cors <- cors[!is.na(cors)]
                sum(lo <= cors & cors < hi)
            },
            bin_lo, bin_hi
        )
    ) %>%
    mutate(
        n_prop = n / sum(.$n)
    )
label_order <- cor_bins$label
cor_bins$label <- factor(cor_bins$label, levels = label_order)
cor_bins_permutations <- lapply(
    1:100,
    function(p) {
        tibble(
            bin_lo = seq(-1, 0.9, by = 0.1),
            bin_hi = seq(-0.9, 1, by = 0.1)
        ) %>%
            mutate(label = paste(bin_lo, bin_hi, sep = " ~ ")) %>%
            mutate(
                label = factor(
                    label, levels = factor(cor_bins$label, levels = label_order)
                )
            ) %>%
            mutate(
                n = mapply(
                    function(lo, hi) {
                        cors <- unlist(data_core_permutations[[p]])
                        cors <- cors[!is.na(cors)]
                        sum(lo <= cors & cors < hi)
                    },
                    bin_lo, bin_hi
                )
            ) %>%
            mutate(
                n_prop = n / sum(.$n)
            ) %>%
            mutate(Permutation = p)
    }
) %>%
    bind_rows()

null_colour = "grey"
graph_stroke = 0.5
gtex_plot <- ggplot() +
    geom_point(
        data = cor_bins_permutations,
        mapping = aes(x = label, y = n_prop),
        shape = 21, size = 1, alpha = 0.1, 
        stroke = 0.25, fill = null_colour, colour = "black"
    ) +
    geom_line(
        data = cor_bins_permutations,
        mapping = aes(x = label, y = n_prop, group = Permutation),
        stat = "smooth", method = "loess", 
        size = graph_stroke, colour = null_colour
    ) +
    geom_point(
        data = cor_bins,
        mapping = aes(x = label, y = n_prop),
        shape = 21, size = 1,
        stroke = 0.25,  fill = "red", colour = "black"
    ) +
    geom_line(
        data = cor_bins,
        mapping = aes(x = label, y = n_prop, group = 1),
        stat = "smooth", method = "loess", colour = "red", size = graph_stroke
    ) +
    annotate(
        geom = "text",
        label = "Mean Corr. Identical > Mean Corr. Random; P = 0.01, n = 100",
        x = 8, y = 0.4, size = 2, vjust = 0
    ) +
    annotate(
        geom = "text",
        label = "Correlation of random pairs of isoforms betw. BrainSpan and GTEx",
        x = 8, y = 0.125, size = 2, colour = "grey"
    ) +
    annotate(
        geom = "text",
        label = "Correlation of identical pairs of isoforms betw. BrainSpan and GTEx",
        x = 8, y = -0.01, size = 2, colour = "red"
    ) +
    labs(
        x = "Pearson Correlation Bin",
        y = "Proportion of Isoforms in\nEach Correlation Bin"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_line(size = 0.1)
    )

gtex_plot

ggsave(
    filename = "data/figures/GTExComparison.pdf",
    plot = gtex_plot,
    width = 12, height = 5.5, units = "cm",
    useDingbats = FALSE
)

