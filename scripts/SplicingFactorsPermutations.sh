#!/usr/bin/env bash

#SBATCH -t 24:00:00
#SBATCH -p shared
#SBATCH --mem=80g
#SBATCH --cpus-per-task=12
#SBATCH --job-name=SplicingFactorsPermutations
#SBATCH --output=log/SplicingFactorsPermutations_%a.log
#SBATCH --array=1-70

Rscript ~/BrainSpan/scripts/SplicingFactorsPermutations.R
