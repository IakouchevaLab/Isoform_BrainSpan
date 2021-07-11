library(tidyverse)
library(doParallel)

registerDoParallel(cores = detectCores())

permutation_prefix <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(permutation_prefix)) {
    permutation_prefix <- 0
}
permutations <- as.numeric(paste0(permutation_prefix, seq(0, 99)))



#brainspan_group <- readRDS("data/gtex_analysis/BrainSpanTPM_matchGTEx.rds")
#gtex_group <- readRDS("data/gtex_analysis/GTExTPM_matchBrainSpan.rds")
#samples <- colnames(brainspan_group)[colnames(brainspan_group) != 'transcript_id']