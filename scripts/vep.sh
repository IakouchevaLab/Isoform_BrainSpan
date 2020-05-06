#!/usr/bin/env bash

 vep -i ../data/SNVs/SatterstromDeNovoVariants.vcf -o ../data/SNVs/SatterstromDeNovoVariants.vep \
     --force_overwrite --species homo_sapiens --cache --offline --fasta ../data/source/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
     --plugin GeneSplicer,/Users/kkhaichau/bin/GeneSplicer/sources/genesplicer,/Users/kkhaichau/bin/GeneSplicer/human \
     --custom ../data/source/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH
