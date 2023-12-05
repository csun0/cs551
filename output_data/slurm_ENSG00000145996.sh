#!/bin/bash
#SBATCH --mincpus 4
#SBATCH --mem 24G
#SBATCH --time 1-0:00:00
#SBATCH --job-name csl_ENSG00000145996

Rscript R_ENSG00000145996.R
