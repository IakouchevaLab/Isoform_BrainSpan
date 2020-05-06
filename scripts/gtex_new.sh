#!/usr/bin/env bash

#SBATCH -t 05:00:00
#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -J gtex_new
#SBATCH -o logs/gtex_new.log

Rscript scripts/gtex_new.R