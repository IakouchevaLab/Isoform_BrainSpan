#!/usr/bin/env bash

#SBATCH -t 00:30:00
#SBATCH -p shared
#SBATCH --cpus-per-task=12
#SBATCH --mem=64g
#SBATCH --job-name=SplicingFactorsPermutations
#SBATCH --output=log/SplicingFactorsPermutations_%a.out
#SBATCH --array=1-938

/home/kkchau/usr/bin/Rscript ~/BrainSpan/scripts/SplicingFactorsSelectedPermutations.R
