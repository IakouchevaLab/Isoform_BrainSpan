library(tidyverse)

empirical_correlations <- tibble(
    Group = "Empirical",
    x = cut(as.numeric(readRDS("data/psychencode/BrainSpan-CMC-HBCC_empirical_correlations.rds")), breaks = seq(-1, 1, by = 0.1))
) %>%
    group_by(Group) %>%
    count(x)
null_correlations <- lapply(
    list.files("data/psychencode/", "*null_correlations_*", full.names = TRUE),
    function(rds) {
        readRDS(rds)
    }
)
null_correlations <- bind_rows(null_correlations)

ggplot(
    data = bind_rows(
        empirical_correlations %>% 
            mutate(Iteration = -1),
        null_correlations %>% 
            mutate(Iteration = as.numeric(IterGroup) * Iter) %>%
            select(x, n, Iteration) %>%
            mutate(Group = 'Null')
    ) %>%
        filter(!is.na(x)),
    aes(
        x = x, y = n, group = Iteration, colour = Group, 
        alpha = ifelse(Group == 'Empirical', 1, 0.5)
    )
) +
    geom_point(size = 5) +
    geom_line() +
    labs(x = 'Interval', y = 'Frequency') +
    scale_colour_manual(
        values = c('Empirical' = 'red', 'Null' = 'grey')
    ) +
    guides(alpha = FALSE, colour = guide_legend(title = '')) +
    theme_bw() +
    theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, hjust = 1)
    )

empirical_data <- as.numeric(readRDS("data/psychencode/BrainSpan-CMC-HBCC_empirical_correlations.rds"))
empirical_mean <- mean(empirical_data[!is.na(empirical_data)])
null_means <- as.numeric(sapply(seq(1, 10), function(x) {
    l <- readRDS(paste0("data/psychencode/BrainSpan-CMC-HBCC_null_correlations_", x, ".rds"))
    sapply(seq(1, 100), function(i) {
        mean(as.numeric(l[[i]])[!is.na(as.numeric(l[[i]]))])
    })
}))

sum(empirical_mean >= null_means)

