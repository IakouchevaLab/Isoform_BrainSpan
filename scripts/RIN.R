library(tidyverse)

prefilter_metadata <- read_tsv("data/source/brainSpan.phenotype.meta.final.tsv")
metadata <- read_csv("data/brainspan_metadata.csv")

plt <- ggplot(
    data = metadata,
    mapping = aes(
        x = "", y = RIN
    )
) +
    geom_boxplot() +
    theme_bw() +
    theme(
        text = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    )
ggsave(filename = "data/figures/BrainSpanRIN.pdf",
       plot = plt,
       device = "pdf",
       width = 4, height = 3)
mean(metadata$RIN)
median(metadata$RIN)
ggplot(
    data = bind_rows(
        prefilter_metadata %>%
            mutate(Filter = "Pre-Filter"),
        metadata %>%
            mutate(Filter = "Post-Filter")
    ),
    mapping = aes(
        x = Filter, y = RIN
    )
) +
    geom_boxplot()
