#!/usr/bin/env bash

# Check SNVs against gnomad exome database

GNOMAD=../../data/source/gnomad.exomes.r2.1.1.sites.vcf.bgz
SNVS=../../data/SNVs/allSNVs_Masterfile_REACH_SSC.tsv
GNOMAD_RESULTS=../../data/SNVs/allSNVs_gnomad.txt

DISTINCT_SNVS=../../data/SNVs/allSNVs_Masterfile_REACH_SSC_VariantsOnly.tsv

# cut -f1,2 $SNVS | sort | uniq > $DISTINCT_SNVS

while read snvs; do
    variant=$(echo $snvs | cut -f1 -d ' ' | cut -f1 | sed -E "s/-[0-9]+//g");
    location=$(echo $variant | grep -Eo "[0-9]+$");
    variant="$variant-$location";
    tabix $GNOMAD $variant;
done < $DISTINCT_SNVS > $GNOMAD_RESULTS
