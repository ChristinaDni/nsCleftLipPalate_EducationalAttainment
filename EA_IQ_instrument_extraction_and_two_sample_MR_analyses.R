#Instrument extraction process and two sample MR analyses: i. EA on nsCLP, ii. IQ on nsCLP


#Import the GWAS data files#

EA_Lee <- read.csv("~/Cleft_EA/Datasets/EA_Lee.gz")
OFC_metal_het_20191 <- read.delim("~/Cleft_EA/Datasets/OFC_metal_het_20191.txt")

#Identify SNPs in EA overlaping with nsCLP~

EA_common<- EA_Lee[(EA_Lee$MarkerName%in%OFC_metal_het_20191$MarkerName),]

#Format the data to be processed from MR Base R package#

library(TwoSampleMR)

EA_common<- format_data(EA_common, type = "exposure", snps = NULL, 
                        header = TRUE, snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", 
                        eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", 
                        pval_col = "Pval")

#Identify Instruments#

EA_i<- EA_common [which (EA_common$eaf.exposure >= 0.01), ]
EA_i<- EA_i [which (EA_i$pval.exposure < 5e-08), ]
EA_i<- clump_data(EA_i, clump_r2 = 0.01)

#Estimate the mean F statistic of the instruments#

BetaXG   = EA_i$beta.exposure
seBetaXG = EA_i$se.exposure

F   = BetaXG^2/seBetaXG^2
mF  = mean(F)

#Extract effect sizes and standard errors from the outcome#

nsCLP_out<- read_outcome_data(snps = EA_i$SNP, filename = "OFC_metal.csv", 
                              sep = ",", snp_col = "MarkerName", beta_col = "Effect", se_col = "StdErr", 
                              effect_allele_col = "Allele1", other_allele_col = "Allele2",
                              pval_col = "P.value")

#Data harmonization#

EA_nsCLP<- harmonise_data(EA_i, nsCLP_out)

#Analysis#

res_EA_nsCLP<- mr(EA_nsCLP)

#Import the GWAS data files#

GWAS_CP_all <- read.delim("~/Cleft_EA/Datasets/GWAS_CP_all.txt")
OFC_metal_het_20191 <- read.delim("~/Cleft_EA/Datasets/OFC_metal_het_20191.txt")

#Identify SNPs in EA overlaping with nsCLP~

IQ_common<- GWAS_CP_all[(GWAS_CP_all$MarkerName%in%OFC_metal_het_20191$MarkerName),]

#Format the data to be processed from MR Base R package#

library(TwoSampleMR)

IQ_common<- format_data(IQ_common, type = "exposure", snps = NULL, 
                        header = TRUE, snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", 
                        eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", 
                        pval_col = "Pval")

#Identify Instruments#

IQ_i<- IQ_common [which (IQ_common$eaf.exposure >= 0.01), ]
IQ_i<- IQ_i [which (IQ_i$pval.exposure < 5e-08), ]
IQ_i<- clump_data(IQ_i, clump_r2 = 0.01)

#Estimate the mean F statistic of the instruments#

BetaXG   = IQ_i$beta.exposure
seBetaXG = IQ_i$se.exposure

F   = BetaXG^2/seBetaXG^2
mF  = mean(F)

#Extract effect sizes and standard errors from the outcome#

nsCLP_out<- read_outcome_data(snps = IQ_i$SNP, filename = "OFC_metal.csv", 
                              sep = ",", snp_col = "MarkerName", beta_col = "Effect", se_col = "StdErr", 
                              effect_allele_col = "Allele1", other_allele_col = "Allele2",
                              pval_col = "P.value")

#Data harmonization#

IQ_nsCLP<- harmonise_data(IQ_i, nsCLP_out)

#Analysis#

res_IQ_nsCLP<- mr(IQ_nsCLP)

