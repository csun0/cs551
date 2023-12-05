library(data.table)
setwd("/scratch/network/cs9095/cs551/")
library(dplyr)
library(tidyverse)
library(glue)

joint <- fread("eQTL_data/QTD000085.all.tsv.gz")
joint <- joint[complete.cases(joint[, 'rsid'])]
joint <- joint %>% select(gene_id, rsid, beta, se, pvalue)
dimnames(joint) <- list(NULL, c("gene_id", "rsid", "beta_0", "SE_0", "pvalue_0"))
joint <- joint %>%
    distinct(gene_id, rsid, .keep_all = TRUE)

files <- c('QTD000090.all.tsv.gz',
    'QTD000116.all.tsv.gz',
    'QTD000121.all.tsv.gz',
    'QTD000226.all.tsv.gz',
    'QTD000231.all.tsv.gz',
    'QTD000266.all.tsv.gz',
    'QTD000276.all.tsv.gz',
    'QTD000281.all.tsv.gz',
    'QTD000296.all.tsv.gz',
    'QTD000321.all.tsv.gz',
    'QTD000331.all.tsv.gz',
    'QTD000341.all.tsv.gz',
    'QTD000534.all.tsv.gz',
    'QTD000554.all.tsv.gz')
files <- paste("eQTL_data/", files, sep="")

for (i in seq_along(files)) {
    dt <- fread(files[i])
    dt <- dt[complete.cases(dt[, 'rsid'])]
    dt <- dt %>% select(gene_id, rsid, beta, se, pvalue)

    dimnames(dt) <- list(NULL, c("gene_id", "rsid", glue("beta_{i}"), glue("SE_{i}"), glue("pvalue_{i}")))
    dt <- dt %>%
        distinct(gene_id, rsid, .keep_all = TRUE)

    joint <- merge(joint,dt,by=c("gene_id","rsid"))
}

fwrite(joint, "eQTL_data/eQTL_joint_rsid.csv.gz", compress="gzip", row.names=FALSE)
