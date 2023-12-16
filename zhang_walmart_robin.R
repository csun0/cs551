
library(stringr)
library(lme4)
library(mvnfast)
library(glue)
# MR-MtRobin

walmart_robin = function(robin_df, n_replicates, ld) {
  #' Perform a (weighted) mixed-effects linear regression relating GWAS summary statistics 
  #' to eQTL summary statistics for any particular gene, assuming the necessary filtering 
  #' steps for valid cross-tissue IVs have already been carried out. Also calculates a p-value
  #' through repeated sampling of null GWAS effect vectors. 
  #' 
  #' @param robin_df a dataframe where each row is a unique SNP-tissue combination, with
  #' colNames(robin_df) = c('snp', 'gwas_b', 'gwas_se', 'eqtl_b', eqtl_se')
  #' @param n_replicates number of null replicates to sample
  #' @param ld I x I matrix of linkage disequilibrium scores between each unique SNP, assuming
  #' there are I SNPs 
  #' 
  #' @return p_val p-value for gene-trait association
  
  eqtl_se = robin_df$eqtl_se
  eqtl_b = robin_df$eqtl_b
  gwas_b = robin_df$gwas_b
  snp = robin_df$snp
  
  ld_len <- length(unique(snp))
  # LD <- diag(ncol=ld_len, nrow=ld_len) test with diagonal covariance matrix

  # Model does not fit an intercept, either random or fixed
  me_fit = lmer(eqtl_b~-1+gwas_b+(gwas_b-1|snp), 
                weights=1/eqtl_se^2)
  t_stat = summary(me_fit)$coefficients[3]
  
  # Simulate null vectors of GWAS summary stats
  temp = subset(robin_df, select = -c(eqtl_b, eqtl_se))
  temp <- temp[!duplicated(temp)]
  gwas_se = temp$gwas_se

  sigma_gwas = (gwas_se %*% t(gwas_se)) * ld
  null_gwas_b = rmvn(n=n_replicates, mu=rep(0, length(gwas_se)), 
                     sigma=sigma_gwas) # n_replicates x p (unique SNPs) matrix
  colnames(null_gwas_b) = unique(robin_df$snp)
  t_nulls = numeric(n_replicates)
  # print(null_gwas_b)
  
  for (r in 1:n_replicates) {
    print(glue("run {r}"))
    # Suppress convergence issue warnings 
    suppressWarnings({null_fit = lmer(eqtl_b~-1+null_gwas_b[r,snp]+
                                        (null_gwas_b[r,snp]-1|snp),
                                      weights=1/eqtl_se^2)})
    t_nulls[r] = summary(null_fit)$coefficients[3]
  }
  
  p_val = sum(abs(t_nulls) > abs(t_stat))/n_replicates
  return(p_val)
}
