library(tidyverse)

empirical_data <- read_tsv("data/gtex_analysis/empirical_correlations.tsv") %>%
    filter(!is.na(estimate)) %>% {
        est <- pull(., estimate)
        return(data.frame(table(cut(est, breaks = seq(-1, 1, by = 0.1)))))
    } %>%
    mutate(permutation = -1)
null_data <- bind_rows(lapply(list.files("data/gtex_analysis/null_correlations/", full.names = T), readRDS))

total_data <- bind_rows(
    empirical_data %>%
        mutate(Group = "Empirical"),
    null_data %>%
        mutate(Group = "Null")
)
correlation_plot <- ggplot(
    data = total_data,
    mapping = aes(
        x = Var1, y = Freq, group = permutation, colour = Group, 
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
# ggsave(
#     filename = 'data/gtex_analysis/correlation_comparison.pdf',
#     plot = correlation_plot,
#     width = 12, height = 12, device = 'pdf'
# )

# Pvalues

read_tsv("data/gtex_analysis/empirical_correlations.tsv") %>%
    filter(!is.na(estimate)) %>%
    pull(estimate) %>%
    mean()
