#LDSC Analyses#

#Format the nsCLP GWAS file#
./munge_sumstats.py \
--out nsCLP \
--merge-alleles w_hm3.snplist \
--N 5951.0 \
--sumstats OFC_metal_het_20191.txt 

#Format the Educational attainment GWAS file#
./munge_sumstats.py \
--out EA \
--merge-alleles w_hm3.snplist \
--N 766345.0 \
--sumstats GWAS_EA_excl23andMe.txt 

#Format the Intelligencce GWAS file#
./munge_sumstats.py \
--out IQ \
--merge-alleles w_hm3.snplist \
--N 257828.0 \
--sumstats GWAS_CP_all.txt 

#Format the Philtrum width GWAS file#
./munge_sumstats.py \
--out philtrum \
--merge-alleles w_hm3.snplist \
--N 6136.0 \
--sumstats PhiltrumSummaryStats.txt 


#Estimate the genetic correlation between nsCLP and EA#
./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out EA_cleft \
--rg cleft.sumstats.gz,EA.sumstats.gz \
--w-ld-chr eur_w_ld_chr/ 

  
#Estimate the genetic correlation between nsCLP and IQ#
  ./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out IQ_cleft \
--rg cleft.sumstats.gz,IQ.sumstats.gz \
--w-ld-chr eur_w_ld_chr/ 


#Estimate the genetic correlation between nsCLP and Philtrum width (positive control analysis)#
  ./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out cleft_philtrum \
--rg cleft.sumstats.gz,philtrum.sumstats.gz \
--w-ld-chr eur_w_ld_chr/ 
  
  
