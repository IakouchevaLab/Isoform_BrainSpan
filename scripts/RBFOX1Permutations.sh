#!/usr/bin/env bash

#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --mem=32g
#SBATCH --cpus-per-task=12
#SBATCH --job-name=RBFOX1Permutations
#SBATCH --output=log/RBFOX1Permutations_%a.log
#SBATCH --array=1-70

Rscript ~/BrainSpan/scripts/RBFOX1Permutations.R
