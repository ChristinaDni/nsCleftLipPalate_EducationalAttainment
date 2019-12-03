#Perform two-sample MR analysis: nsCLP-EA

setwd("~/Cleft_EA/Datasets")

library(TwoSampleMR)

#Format the data to be processed from MRBase R package

nsCLP_i<- format_data(nsCLP_i, type = "exposure", snps = NULL, header = TRUE, 
                      snp_col = "SNP", beta_col = "Effect", se_col = "StdErr", effect_allele_col = "Allele1", 
                      other_allele_col = "Allele2", pval_col = "P.value")

#Estimate mean F statistic of the instruments

BetaXG   = nsCLP_i$beta.exposure
seBetaXG = nsCLP_i$se.exposure

F   = BetaXG^2/seBetaXG^2
mF  = mean(F)


#Extract effect sizes and SEs from the outcome

EA_out<- read_outcome_data(snps = nsCLP_i$SNP, filename = "EA_Lee", 
                           sep = ",", phenotype_col= "", snp_col = "MarkerName", beta_col = "Beta", 
                           se_col = "SE", eaf_col = "", effect_allele_col = "A1", 
                           other_allele_col = "A2", pval_col = "P")

#Data harmonisation

nsCLP_EA<- harmonise_data(nsCLP_i, EA_out)

#Analysis

res_nsCLP_EA<- mr(nsCLP_EA)


#Perform two-sample MR analysis: nsCLP-IQ

#Extract effect sizes and SEs from the outcome

IQ_out<- read_outcome_data(snps = nsCLP_i$SNP, filename = "IQ_Lee", 
                           sep = ",", phenotype_col= "", snp_col = "MarkerName", beta_col = "Beta", 
                           se_col = "SE", eaf_col = "", effect_allele_col = "A1", 
                           other_allele_col = "A2", pval_col = "P")

#Data harmonisation

nsCLP_IQ<- harmonise_data(nsCLP_i, IQ_out)

#Analysis

res_nsCLP_IQ<- mr(nsCLP_IQ)

