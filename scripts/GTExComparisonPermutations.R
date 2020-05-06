library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores())

permutation_prefix <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(permutation_prefix)) {
    permutation_prefix <- 0
}
permutations <- as.numeric(paste0(permutation_prefix, seq(0, 99)))

brainspan_group <- readRDS("data/gtex_analysis/BrainSpanTPM_matchGTEx.rds")
gtex_group <- readRDS("data/gtex_analysis/GTExTPM_matchBrainSpan.rds")
samples <- colnames(brainspan_group)[colnames(brainspan_group) != 'transcript_id']

null_correlations <- foreach(i=permutations) %dopar% {
    brainspan_input <- brainspan_group
    gtex_input <- gtex_group[sample(nrow(gtex_group)), ]
    correlation_data <- mclapply(1:nrow(brainspan_input), function(tx) {
        brainspan_data <- as.numeric(brainspan_input[tx, samples])
        gtex_data <- as.numeric(gtex_input[tx, samples])
        correlation <- cor.test(brainspan_data, gtex_data, method = "pearson") %>%
            broom::tidy()
        return(correlation)
    }) %>%
        bind_rows() %>%
        filter(!is.na(estimate))
    data.frame(table(cut(correlation_data$estimate, breaks = seq(-1, 1, by = 0.1)))) %>%
        mutate(permutation = i) %>%
        return()
} %>%
    bind_rows()
dir.create("data/gtex_analysis/null_correlations/")
saveRDS(
    null_correlations,
    paste0("data/gtex_analysis/null_correlations/null_correlations_permutation", permutation_prefix, ".rds")
)