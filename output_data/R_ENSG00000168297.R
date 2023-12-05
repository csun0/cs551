version
library(glue)
library(data.table)
library(stringr)
library(lme4)
library(mvnfast)
setwd("/scratch/network/cs9095/cs551")
source("walmart_robin.R")
print("prelim complete~")

robin_df <- fread("merged_data/ENSG00000168297_merge.tsv")
names(robin_df) <- c("snp", "gwas_b", "gwas_se", "eqtl_b", "eqtl_se")
setorder(robin_df, cols="snp")
n_replicates = 1000
output <- walmart_robin(robin_df, n_replicates, ld)
fwrite(list(output), "output_data/ENSG00000168297_pval")
print("done!")
