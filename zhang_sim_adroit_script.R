#!/usr/bin/env Rscript
args = commandArgs(TRUE)
args = as.numeric(args)

### Run population structure simulations for MR-MtRobin

# Imports 
.libPaths('/scratch/network/az2747/software/R/lib')
library(MASS)
library(Matrix)
library(bnpsd)
library(mvnfast)
library(reshape2)
library(glmnet)
library(stats)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(lme4)

# Functions
sim_indv_data = function(gamma, N_g, N_R, K, m_loci=100, Q=10, M=2, 
                         k_subpops=6, inbr_subpops=c(0.15, 0.3, 0.5, 0.25, 0.2, 0.21)) {
  #' Simulate individual-level data for a GWAS and eQTL study. This includes
  #' genotype matrices, as well as phenotypes (trait vector and gene expression matrix
  #' for GWAS and eQTL respectively).
  
  N = N_g + N_R
  gwas_id = sample(1:N, N_g)
  eqtl_id = setdiff(1:N, gwas_id)
  
  # Simulate admixture model with 1D geography, k=6
  p_anc = draw_p_anc(m_loci, p_min=0.1)
  p_subpops = draw_p_subpops(p_anc, inbr_subpops)
  admix_proportions =  admix_prop_1d_linear(N, k_subpops, sigma = 0.5)
  p_ind = make_p_ind_admix(p_subpops, admix_proportions) 
  L = t(draw_genotypes_admix(p_ind)) # Full genotype matrix
  
  # Redraw genotypes if any locus is entirely fixed 
  while ((length(which(apply(L[gwas_id,], 2, function(x) var(x)) > 0)) < m_loci) | 
         (length(which(apply(L[eqtl_id,], 2, function(x) var(x)) > 0)) < m_loci)) {
    p_anc = draw_p_anc(m_loci, p_min=0.1)
    p_subpops = draw_p_subpops(p_anc, inbr_subpops)
    admix_proportions =  admix_prop_1d_linear(N, k_subpops, sigma = 0.5)
    p_ind = make_p_ind_admix(p_subpops, admix_proportions) 
    L = t(draw_genotypes_admix(p_ind))
  }
  
  ## GWAS data generation
  signal = sample(1:m_loci, Q, replace=FALSE) 
  G = L[gwas_id, signal] # Matrix of true eSNPs
  invalid_id = sample(1:Q, M)
  H = G[,invalid_id] # Matrix of invalid IVs 
  
  ### Effect sizes of eSNPs
  S_mu_x = diag(x=rep(0.5,Q), ncol=Q)
  mu_x = rmvn(n=1, mu=rep(0,Q), sigma=S_mu_x)
  
  ### Latent confounder
  Z = rnorm(N_g)
  
  ### Gene expression levels
  eta_x = runif(N_g, min=0, max=0.1)
  e_x = as.matrix(rnorm(N_g), ncol=1)
  X = G%*%t(mu_x) + as.matrix(eta_x*Z, ncol=1) + e_x
  
  ### Trait vector
  S_mu_y = diag(x=rep(0.5,M), ncol=M)
  mu_y = rmvn(n=1, mu=rep(0,M), sigma=S_mu_y)
  eta_y = runif(N_g, min=0, max=0.1)
  e_y = as.matrix(rnorm(N_g), ncol=1)
  
  Y = gamma*X + H%*%t(mu_y) + as.matrix(eta_y*Z) + e_y
  
  ## eQTL data generation 
  G_R = L[eqtl_id, signal]
  X_R = matrix(nrow=N_R, ncol=K)
  
  for (k in 1:K) {
    mu_x_k = t(mu_x) + as.matrix(rnorm(n=K, sd=0.3), ncol=1)
    theta_k = rbinom(n=K, size=1, prob=0.5)
    X_R_k = G_R%*%(mu_x_k*theta_k) + as.matrix(rnorm(N_R, sd=0.3), ncol=1)
    X_R[,k] = X_R_k
  }
  
  return(list(Y=Y, X_R=X_R, eSNP_gwas=L[gwas_id,], eSNP_eqtl=L[eqtl_id,]))
}


fix_na = function(v){
  
  sum = summary(v)
  eqtl_b_k = as.numeric(sum$coefficients[-1,1])
  eqtl_se_k = as.numeric(sum$coefficients[-1,2])
  eqtl_p_k = as.numeric(sum$coefficients[-1,4])
  
  if(!any(is.na(v))){
    idx = as.numeric(which(is.na(v$coefficients)))
    for (id in idx) {
      id = id - 1
      eqtl_b_k = c(eqtl_b_k[1:(id-1)], 0, tail(eqtl_b_k, -(id-1)))
      eqtl_se_k = c(eqtl_se_k[1:(id-1)], 0, tail(eqtl_se_k, -(id-1)))
      eqtl_p_k = c(eqtl_p_k[1:(id-1)], 1, tail(eqtl_p_k, -(id-1)))
    }
  } 
  
  return(list(eqtl_b_k=eqtl_b_k, eqtl_se_k=eqtl_se_k, eqtl_p_k=eqtl_p_k))
}


est_eqtl_effects = function(X_R, eSNP_eqtl) {
  #' Produce eQTL summary statistics and keep significant IVs. 
  
  K = ncol(X_R)
  m_loci = ncol(eSNP_eqtl)
  eqtl_df = data.frame(matrix(ncol=5, nrow=0))
  names(eqtl_df) = c('snp', 'beta', 'se', 'pvalue', 'tissue')
  
  for (k in 1:K) {
    ols_eqtl = lm(X_R[,k] ~ eSNP_eqtl)
    coef = fix_na(ols_eqtl)
    eqtl_b_k = coef$eqtl_b_k
    eqtl_se_k = coef$eqtl_se_k
    eqtl_p_k = coef$eqtl_p_k
    df_tissue = data.frame(snp=seq(1:m_loci),
                           eqtl_b=eqtl_b_k, eqtl_se=eqtl_se_k, 
                           pvalue=eqtl_p_k,
                           tissue=rep(k, m_loci))
    eqtl_df = rbind(eqtl_df, df_tissue)
  }
  
  valid_iv_neg = eqtl_df[(eqtl_df$pvalue < 0.001 & eqtl_df$eqtl_b < 0),]
  valid_iv_neg = valid_iv_neg %>% group_by(snp) %>% filter(n() >= 3)
  
  valid_iv_pos = eqtl_df[(eqtl_df$pvalue < 0.001 & eqtl_df$eqtl_b > 0),]
  valid_iv_pos = valid_iv_pos %>% group_by(snp) %>% filter(n() >= 3)
  
  iv_negs = unique(valid_iv_neg$snp)
  iv_pos = unique(valid_iv_pos$snp)
  
  overlap_iv = intersect(iv_negs, iv_pos)
  if (length(overlap_iv) > 0) {
    for (i in overlap_iv) {
      n_pos = length(which(valid_iv_pos$snp == overlap_iv))
      n_neg = length(which(valid_iv_neg$snp == overlap_iv))
      
      if (n_pos > n_neg) {
        valid_iv_neg = valid_iv_neg[valid_iv_neg$snp != i,]
      } else {
        valid_iv_pos = valid_iv_pos[valid_iv_pos$snp != i,]
      }
    }
  }
  
  full_iv = rbind(valid_iv_neg, valid_iv_pos)
  full_iv = full_iv %>% arrange(snp)
  
  return(full_iv)
}


est_gwas_effects = function(Y, eSNP_gwas) {
  #' Produce GWAS summary statistics.
  
  ols_gwas = lm(Y ~ eSNP_gwas)
  ols_gwas_out = summary(ols_gwas)
  gwas_b = ols_gwas_out$coefficients[-1,1]
  gwas_se = ols_gwas_out$coefficients[-1,2]
  gwas_df = data.frame(snp=seq(1:ncol(eSNP_gwas)), gwas_b=gwas_b, gwas_se=gwas_se)
  rownames(gwas_df) = NULL
  
  return(gwas_df)
}


source('./walmart_robin.R')
run_sim_for_k = function(n_subpop, n_sim, M, n_replicates) {
  
  N_g = 3000 # GWAS cohort
  N_R = 200 # eQTL cohort
  K = 10 # Tissue types
  k_subpops = seq(8)
  inbr_subpops = c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  p_signal = 0.2 # Proportion of causal genes 
  
  fdr = 0.1
  
  sig_idx = sample(seq(n_sim), as.integer(p_signal*n_sim))
  null_idx = setdiff(seq(n_sim), sig_idx)
  p_vec = numeric(n_sim)
  
  for (tr in 1:n_sim) {
    if (tr %in% sig_idx) {
      gamma = 0.5
    } else {
      gamma = 0
    }
    indv_data = sim_indv_data(gamma=gamma, N_g=N_g, N_R=N_R, K=K, M=M,
                              k_subpops=n_subpop, inbr_subpops=inbr_subpops[1:n_subpop])
    eqtl = est_eqtl_effects(indv_data$X_R, indv_data$eSNP_eqtl)
    gwas = est_gwas_effects(indv_data$Y, indv_data$eSNP_gwas)
    comb = eqtl %>% inner_join(gwas, by=c('snp'))
    robin_df = comb[,c('snp', 'gwas_b', 'gwas_se', 'eqtl_b', 'eqtl_se')]
    
    # Calculate linkage disequilibrium 
    ld = cor(indv_data$eSNP_gwas)
    snp_ids = sort(unique(robin_df$snp), decreasing=FALSE)
    ld = ld[snp_ids, snp_ids]
    
    snp_ids_str = as.character(snp_ids)
    robin_df$snp = as.character(robin_df$snp)
    rownames(ld) = snp_ids_str
    colnames(ld) = snp_ids_str
    
    p = walmart_robin(robin_df=robin_df, n_replicates, ld=ld)
    p_vec[tr] = p
    print(tr)
  }
  
  p_bh_vec = p.adjust(p_vec, method='BH')
  gene_disc = which(p_bh_vec < fdr)
  
  if (length(gene_disc) > 0) {
    FDP_bh = length(intersect(gene_disc, null_idx))/(length(gene_disc) + 1)
    power_bh = length(intersect(gene_disc, sig_idx))/max(1, length(sig_idx))
  } else {
    FDP_bh = 0
    power_bh = 0
  }

  return(list(FDP=FDP_bh, power=power_bh))
}


# Run script
run_sim_for_k(n_subpop=args[1], n_sim=args[2], M=args[3], n_replicates=args[4])
