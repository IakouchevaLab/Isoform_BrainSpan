library(tidyverse)
library(readxl)

# Master File
# Going to have to run VEP on this data
masterfile <- file.path(
    "data", "source",
    "MasterFile-AllGenomicNoIntergenic_mutations_cases_20181104_KC.xlsx"
)
master_case <- read_xlsx(
    masterfile, sheet = "DNM_Asd-Scz-ID-EE_29S", skip = 1, 
    col_types = "text"
) %>%
    filter(
        `Position(hg19)start_left_normalized` == `Position(hg19)end_left_normalized`
    ) %>%
    dplyr::select(
        Disorder, Chromosome, ends_with("normalized"),
        Reference_allele, Mutant_allele
    ) %>%
    filter(
        str_detect(Disorder, "ASD")
    ) %>%
    rename(
        Chromosome = Chromosome,
        start = `Position(hg19)start_left_normalized`,
        end = `Position(hg19)end_left_normalized`,
        ref = Reference_allele,
        mut = Mutant_allele
    )
master_control <- read_xlsx(
    masterfile, sheet = "controls", col_types = "text"
) %>%
    filter(
        `Position(hg19)start_left_normalized` == `Position(hg19)end_left_normalized`
    ) %>%
    rename(
        Chromosome = Chromosome,
        start = `Position(hg19)start_left_normalized`,
        end = `Position(hg19)end_left_normalized`,
        ref = Reference_allele,
        mut = Mutant_allele
    ) %>%
    dplyr::select(
        Chromosome, start, end, ref, mut
    ) %>%
    mutate(Disorder = NA)
master <- bind_rows(
    master_case %>%
        mutate(Phenotype = "Case"),
    master_control %>%
        mutate(Phenotype = "Control")
) %>%
    mutate(
        identifier = paste(
            Chromosome, start, end, ref, mut, Phenotype, sep = ":"
        )
    ) %>%
    mutate(Placeholder = ".") %>%
    mutate(Pass = "PASS")

# Write vcf

master_vcf <- master %>%
    dplyr::select(
        Chromosome, start, identifier, ref, mut, 
        Placeholder, Pass
    ) %>%
    mutate(Placeholder2 = ".")
colnames(master_vcf) <- c(
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
)

write_tsv(master_vcf, "data/masterfile_asd_vcf.txt")

# Set paths to VEP
# Call from shell
# source ~/.zshrc
# 
# vep -i data/masterfile_asd_vcf.txt -o data/masterfile_asd_vep.txt \
#     --force_overwrite \
#     --species homo_sapiens \
#     --cache \
#     --offline \
#     --fasta data/source/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
#     --plugin GeneSplicer,\
#           /Users/kkhaichau/bin/GeneSplicer/sources/genesplicer,\
#           /Users/kkhaichau/bin/GeneSplicer/human

masterfile_vep <- read_tsv(
    "data/masterfile_asd_vep.txt", skip = 41
)

# REACH/SSC cohorts

ssc <- read_tsv(
    "data/SNVs/annotated_snvs.tsv",
    col_types = cols(.default = "c")
)

# Full mutations list

mutations <- bind_rows(
    masterfile_vep %>%
        mutate(
            Phenotype_1ctrl_2case = ifelse(
                str_detect(`#Uploaded_variation`, "Case"),
                2, 1
            )
        ) %>%
        dplyr::select(
            Location, Allele, Gene, Feature, Feature_type, Consequence, 
            Phenotype_1ctrl_2case
        ) %>%
        mutate(Data = "Mastertable"),
    ssc %>%
        mutate(Phenotype_1ctrl_2case = as.numeric(Phenotype_1ctrl_2case)) %>%
        dplyr::select(
            Location, Allele, Gene, Feature, Feature_type, Consequence, 
            Phenotype_1ctrl_2case
        ) %>%
        mutate(Data = "REACH_SSC")
)

# mutations %>%
#     separate_rows(Consequence, sep = ",") %>%
#     dplyr::select(
#         Location, Consequence, Phenotype_1ctrl_2case, Data
#     ) %>%
#     filter(str_detect(Consequence, "acceptor|donor|stop_gain|start_lost")) %>%
#     count(Data, Phenotype_1ctrl_2case, Consequence) %>%
#     ggplot(
#         data = .,
#         mapping = aes(
#             x = Consequence, y = n, group = Data, fill = Data
#         )
#     ) +
#         facet_grid(. ~ Phenotype_1ctrl_2case) +
#         geom_bar(stat = "identity") +
#         theme(
#             text = element_text(size = 30),
#             axis.text.x = element_text(angle = 30, hjust = 1),
#             plot.margin = margin(l = 8, unit = "lines")
#         )

# Annotate features of mutations

anno <- read_csv("data/source/annotation.transcript.ensg75.txt")[, -1]

# Add extra information to annotation table
ensembl <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "GRCh37.ensembl.org"
)
anno_expand <- anno %>%
    left_join(
        getBM(
            attributes = c("ensembl_gene_id", "start_position", "end_position"),
            filters = "ensembl_gene_id",
            values = unique(.$ensembl_gene_id),
            mart = ensembl
        ),
        by = "ensembl_gene_id"
    ) %>%
    mutate(gene_length = abs(start_position - end_position))

mutations_annotate <- left_join(
    mutations,
    anno_expand,
    by = c("Feature" = "ensembl_transcript_id", "Gene" = "ensembl_gene_id")
) %>%
    distinct()

write_tsv(
    mutations_annotate,
    "data/SNVs/allSNVs_Masterfile_REACH_SSC.tsv"
)

################################################################################
# EDA                                                                          #
################################################################################

lof <- "acceptor|donor|stop_gain|start_lost|frameshift"
mutations <- read_tsv("data/SNVs/allSNVs_Masterfile_REACH_SSC.tsv")

lof_mutations <- mutations %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(str_detect(Consequence, lof)) %>%
    mutate(
        total_n = sapply(
            Phenotype_1ctrl_2case, function(phen) {
                nrow(filter(., Phenotype_1ctrl_2case == phen))
            }
        )
    )

lof_mutation_freq <- ggplot(
    data = lof_mutations %>%
        distinct(Location, Consequence, Phenotype_1ctrl_2case, total_n) %>%
        count(Consequence, Phenotype_1ctrl_2case, total_n) %>%
        mutate(rel_n = n / total_n),
    mapping = aes(
        x = Consequence, y = rel_n,
        group = Phenotype_1ctrl_2case, 
        fill = as.character(Phenotype_1ctrl_2case)
    )
) +
    coord_flip() +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        y = "Relative Mutation Frequency"
    ) +
    scale_fill_discrete(
        labels = c("1" = "Control", "2" = "Case"),
        breaks = c("2", "1")
    ) +
    guides(
        fill = guide_legend(title = "")
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
ggsave(
    filename = "data/figures/LoF_Frequency.pdf",
    plot = lof_mutation_freq,
    device = "pdf", width = 16, height = 12
)

################################################################################
# Enrichment of LoF mutations between case and control                         #
################################################################################

mutations_distinct <- mutations %>%
    filter(!str_detect(Consequence, "NMD|non_coding")) %>%
    separate_rows(Consequence, sep = ",") %>%
    dplyr::select(Location, Consequence, Phenotype_1ctrl_2case) %>%
    mutate(lof = ifelse(str_detect(Consequence, lof), TRUE, FALSE))

fisher.test(
    matrix(c(
        mutations_distinct %>%
            filter(Phenotype_1ctrl_2case == 2 & lof) %>%
            nrow(),
        mutations_distinct %>%
            filter(Phenotype_1ctrl_2case == 2 & ! lof) %>%
            nrow(),
        mutations_distinct %>%
            filter(Phenotype_1ctrl_2case == 1 & lof) %>%
            nrow(),
        mutations_distinct %>%
            filter(Phenotype_1ctrl_2case == 1 & ! lof) %>%
            nrow()
    ), byrow = TRUE, ncol = 2, nrow = 2),
    alternative = "greater"
)
