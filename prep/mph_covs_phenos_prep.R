#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("readr")
#install.packages("data.table")
#install.packages("tidyverse", lib=)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(rrcov)
library(RNOmni)
library(tidyverse)


### 1. PREPARE COVARIATES
qc_df <- read_tsv("samples_with_qc_issues_for_exclusion")

samples <- read.table("genomic_metrics.tsv", sep="\t", header=T)
colnames(samples)[1] <- 'FID'

ancestries <- read.table("ancestry_preds.tsv", sep="\t", header=T)
rm_ancestr = filter(ancestries, ancestries$ancestry_pred != ancestries$ancestry_pred_other)

samples <- samples %>% 
    filter(!(FID %in% qc_df$id)) %>% 
    filter(!(FID %in% rm_ancestr[,1])) %>% 
    mutate(sample_type = as.integer(factor(sample_source)) - 1) 

sexage <- merge(sexage, samples[,c(1,13)], by='FID')


# AFR COVARIATES
afr <- read.table("keep_umap1_afr_fidiid.txt", sep='', col.names=c('FID', 'IID'))
afr_covs <- merge(afr, sexage, by=c('FID', 'IID'))

# add PCs to covs df
pcs <- read.table("cleaned_afr_umap1_pcs.csv", sep='\t', header=T)
pc_cols <- list()
for (i in 1:120) {pc_cols <- append(pc_cols, paste('pc', i, sep=''))}
colnames(pcs) <- pc_cols
pcs$FID <- afr[,1] # label IIDs from list of unfiltered post-umap PCs
afr_covs <- merge(afr_covs, pcs, by='FID') 

write.table(afr_covs, 'afr_covs.txt', sep='\t', quote=F, row.names=F)


# AMR COVARIATES
amr <- read.table("amrQC.postimputation.fam", sep='', header=F, col.names=c('FID', 'IID'))
amr_covs <- merge(amr, sexage, by=c('FID', 'IID'))

# add PCs to covs df
pcs <- read.table("cleaned_amr_umap0_pcs.csv", sep='\t', header=T)
colnames(pcs) <- pc_cols
pcs$FID <- amr[,1] # label IIDs from list of unfiltered post-umap PCs
amr_covs <- merge(amr_covs, pcs, by='FID') 

write.table(amr_covs, 'afr_covs.txt', sep='\t', quote=F, row.names=F)


# EUR COVARIATES
eur <- read.table("keep_umap2_eur_fidiid.txt", sep='', header=F, col.names=c('FID', 'IID'))
eur_covs <- merge(eur, sexage, by=c('FID', 'IID'))

# add PCs to covs df
pcs <- read.table("cleaned_eur_umap2_pcs.csv", sep='\t', header=T)colnames(pcs) <- pc_cols
pcs$FID <- eur[,1] # label IIDs from list of unfiltered post-umap PCs
eur_covs <- merge(eur_covs, pcs, by='FID') 

write.table(eur_covs, 'eur_covs.txt', sep='\t', quote=F, row.names=F)



# 2. REGRESS PHENOTYPES (ONLY FOR QUANTITATIVE PHENOTYPES)
pheno <- read.csv("{anc}_{pheno}_raw.csv", header=T)  # {anc}_{pheno}_raw.csv is the per-ancestry raw phenotype file. Columns are IID, FID, named phenotype
covs <- read.table("{anc}_covs.txt", sep='\t', header=T) # {anc}_covs.txt can be afr_covs.txt, amr_covs.txt, or eur_covs.txt
pheno_covs <- merge(covs, pheno, by=c('FID','IID'))

# rank-inverse normalise phenotype separately by sex
pheno_f <- dplyr::filter(pheno_covs, pheno_covs$sex==2)
pheno_m <- dplyr::filter(pheno_covs, pheno_covs$sex==1)
pheno_f$rint_pheno <- RankNorm(pheno_f[,127])
pheno_m$rint_pheno <- RankNorm(pheno_m[,127])

# then combine
pheno_all <- rbind(pheno_f, pheno_m)
pheno_all <- pheno_all %>% arrange(IID)

# residualise on covariates
# note: "pheno_rintres" should be "{name_of_phenotype}_rintres" to match the final phenotype file (see below)
pheno_all$pheno_rintres <- lm(rint_pheno ~ as.factor(sex) + as.factor(sample_type) 
                         + age + age_sq + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 
                         + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + pc26 + pc27 
                         + pc28 + pc29 + pc30 + pc31 + pc32 + pc33 + pc34 + pc35 + pc36 + pc37 + pc38 + pc39 
                         + pc40 + pc41 + pc42 + pc43 + pc44 + pc45 + pc46 + pc47 + pc48 + pc49 + pc50 + pc51 
                         + pc52 + pc53 + pc54 + pc55 + pc56 + pc57 + pc58 + pc59 + pc60 + pc61 + pc62 + pc63 
                         + pc64 + pc65 + pc66 + pc67 + pc68 + pc69 + pc70 + pc71 + pc72 + pc73 + pc74 + pc75 
                         + pc76 + pc77 + pc78 + pc79 + pc80 + pc81 + pc82 + pc83 + pc84 + pc85 + pc86 + pc87 
                         + pc88 + pc89 + pc90 + pc91 + pc92 + pc93 + pc94 + pc95 + pc96 + pc97 + pc98 + pc99 
                         + pc100 + pc101 + pc102 + pc103 + pc104 + pc105 + pc106 + pc107 + pc108 + pc109 + pc110 
                         + pc111 + pc112 + pc113 + pc114 + pc115 + pc116 + pc117 + pc118 + pc119 + pc120, 
                         data = pheno_all)$residuals
# histogram plot of {pheno} should be almost perfectly normally-distributed

# make rintres phenotype file
# phenotype MPH format require only IID and phenotype columns. Phenotype column must match file name for run_mph_aou.wdl to work
write.table(pheno_all[,c(1,129)], "{pheno}_rintres_{anc}_mph.csv", row.names=F, quote=F, sep=',')


