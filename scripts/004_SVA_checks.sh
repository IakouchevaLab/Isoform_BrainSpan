#!/usr/bin/env bash

#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --mem=50g
#SBATCH --mail-type=all
#SBATCH --mail-user=kkchau@ucsd.edu
#SBATCH --job-name=SVA_checks
#SBATCH --output=log/SVA_checks.out
#SBATCH --ntasks-per-node=12
#SBATCH --array=1-20

export PATH=$HOME/usr/R-3.5.1/bin:$PATH
export PATH="/home/kkchau/cellranger-3.0.0:$PATH"
export PATH=$PATH:$HOME/usr/hdf5-1.10.1/bin
export PATH=$PATH:$HOME/usr/openssl-1.0.2p/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/pcre-8.42/lib:$HOME/usr/xz-5.2.4/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/hdf5-1.10.1/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/miniconda3/lib
export R_LIBS=$HOME/usr/R-3.5.1/lib64/R/library
export BLAS=$HOME/usr/OpenBLAS/lib/libopenblas.a
export ATLAS=$HOME/usr/OpenBLAS/lib/libopenblas.a

#export PATH="/home/kkchau/miniconda3/bin:$PATH"

export OPENBLAS_USE_THREAD=0
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

Rscript ~/BrainSpan/scripts/004_SVA_checks.R

unset LD_PRELOAD
unset LD_LIBRARY_PATH
unset OPENBLAS_USE_THREAD
unset OMP_NUM_THREADS




