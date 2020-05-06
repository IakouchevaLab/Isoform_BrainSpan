#!/usr/bin/env Rscript
#
# Process Satterstrom ASD BioRXiv de novo variants

library(tidyverse)

variants <- readxl::read_xlsx(
    "data/source/Satterstrom2019_DeNovoVariants.xlsx",
    sheet = "De novo variants"
)

# Process to VCF format
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

variants_vcf <- variants %>%
    dplyr::select(
        Variant, Affected_Status, Variant_type
    ) %>%
    mutate(
        VariantID = paste(Variant, Affected_Status, Variant_type, sep = "_")
    ) %>%
    separate(
        Variant, into = c(
            "#CHROM", "POS", "REF", "ALT"
        ), 
        sep = ":"
    ) %>%
    rename(
        ID = VariantID
    ) %>%
    mutate(QUAL = ".", FILTER = "PASS", INFO = ".") %>%
    dplyr::select(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    )
write_tsv(variants_vcf, "data/SNVs/SatterstromDeNovoVariants.vcf")

# Call from shell
# 
# # vep -i data/masterfile_asd_vcf.txt -o data/masterfile_asd_vep.txt \
#     --force_overwrite \
#     --species homo_sapiens \
#     --cache \
#     --offline \
#     --fasta data/source/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
#     --plugin GeneSplicer,\
#           /Users/kkhaichau/bin/GeneSplicer/sources/genesplicer,\
#           /Users/kkhaichau/bin/GeneSplicer/human

# Process the VEP result

lof <- c(
    "splice_donor_variant" = "Splice Donor",
    "splice_acceptor_variant" = "Splice Acceptor",
    "frameshift_variant" = "Frameshift",
    "stop_gained" = "Stop Gain",
    "start_lost" = "Start Loss"
)

# Process extra column

extra <- read_tsv(
    "data/SNVs/SatterstromDeNovoVariants.vep",
    skip = sum(
        grepl("^##", readLines("data/SNVs/SatterstromDeNovoVariants.vep"))
    )
)$Extra %>%
    str_match_all("(\\w+)=") %>%
    lapply(as.data.frame) %>%
    bind_rows() %>%
    as.data.frame() %>%
    pull(2) %>%
    unique()

variant_vep <- read_tsv(
    "data/SNVs/SatterstromDeNovoVariants.vep",
    skip = sum(
        grepl("^##", readLines("data/SNVs/SatterstromDeNovoVariants.vep"))
    )
) %>%
    separate(
        `#Uploaded_variation`, 
        into = c(
            "Chromosome", "Variant_start", "Ref", "Alt", 
            "Affected_status", "Variant_type"
        ),
        sep = ":|_", remove = FALSE
    ) %>%
    mutate(Variant_start = as.numeric(Variant_start)) %>%
    mutate(
        Variant_end = Variant_start + (
            mapply(
                function(ref, alt) {
                    max(nchar(ref), nchar(alt)) - 1
                }, Ref, Alt
            )
        )
    ) %>%
    filter(Feature_type == "Transcript") %>%
    mutate(Extra = str_split(Extra, ";")) %>% 
    unnest(Extra) %>% 
    separate(Extra, into = c("Key", "Value"), sep = "=") %>% 
    rowid_to_column() %>% 
    spread(Key, Value) %>%
    dplyr::select(
        `#Uploaded_variation`, Chromosome, Variant_start, Variant_end, Ref, Alt,
        Variant_type, Affected_status, Gene, Feature, Consequence,
        cDNA_position, CDS_position, Protein_position, Amino_acids, Codons,
        Existing_variation, DISTANCE, FLAGS, GeneSplicer, starts_with("gnom"),
        IMPACT, STRAND
    ) %>%
    rename(ensembl_gene_id = Gene, ensembl_transcript_id = Feature) %>%
    mutate(
        LoF = sapply(
            Consequence, function(cons) any(str_detect(cons, names(lof)))
        )
    )
write_tsv(variant_vep, "data/SNVs/SatterstromProcessedVEP.txt")

# breakdown

lof_consequence_breakdown <- variant_vep %>%
    distinct(`#Uploaded_variation`, Affected_status, LoF, Consequence) %>%
    separate_rows(Consequence, sep = ",") %>%
    filter(LoF) %>%
    distinct(`#Uploaded_variation`, Affected_status, Consequence) %>%
    filter(Consequence %in% names(lof)) %>%
    ggplot(
        data = .,
        mapping = aes(
            x = Affected_status, fill = Consequence
        )
    ) +
    geom_bar(position = "fill", colour = "black") +
    labs(
        x = "Phenotype", y = "Proportion of Variant Consequences",
        title = "Breakdown of LoF Variants"
    ) +
    scale_x_discrete(
        breaks = c(1, 2),
        labels = c(
            "1" = "Control", "2" = "Case"
        )
    ) +
    scale_fill_discrete(
        labels = lof
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 30)
    )
ggsave(
    filename = "data/figures/LoFSNVBreakdown.pdf",
    plot = lof_consequence_breakdown,
    device = "pdf", width = 16, height = 16
)



variants <- readxl::read_xlsx("data/SupplementaryTables/Supplementary Table 7.xlsx", sheet = readxl::excel_sheets("data/SupplementaryTables/Supplementary Table 7.xlsx")[2])
dim(variants)
length(unique(pull(variants, `Ensembl Transcript ID`)))
write_csv(distinct(variants), "data/variants.csv")
