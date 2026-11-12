source("wdls-jsons/mph_functs.R")
library(tidyr)
library(dplyr)
library(data.table)

# using code from https://github.com/jiang18/mph/blob/main/docs/docs/examples.md#by-functional-annotations

anc <- commandArgs(trailingOnly=T)[1]
pheno <- commandArgs(trailingOnly=T)[2]
vcfile_suff <- commandArgs(trailingOnly=T)[3] # suffix for the MPH vc file eg "newvep2_alpha_sepcod"

print("Recalculating fold-enrichments for ", anc, pheno, "...")

# Read the variant info file to get the variant weighting matrix <sw>
sw = fread(paste0("../cat_annots/merged_", anc, "_sepRNA_afternorm.sw"))
head(sw)

sw <- sw[,2:41]
colnames(sw)

# Get the variant incidence matrix (whether variant is present or absent)
# by converting sw into binary format
si = (sw != 0) * 1

# Calculate the crossproduct of 'si' and 'sw'.
# The columns of 'sw' should match the rows of 'vcfile'; otherwise, they should be reordered.
#index = 1:28; sw = sw[, index]
cp = crossprod(as.matrix(si), as.matrix(sw))


# Number of variants in each annotation category of interest
annot_size = colSums(si)
length(annot_size)
# Total number of variants
nvars = nrow(si)

print("SNP/variant weighing matrix dimensions: ", dim(sw))
print("Number of variants: ", nvars)

vcfile = "mph_outputs/", anc, "_", pheno, "_mph_", vcfile_suff, ".mq.vc.csv"
result = recalculate_enrichments(vcfile, cp, nsnps=nvars, annot.size=annot_size) 
result$vc_name = rownames(result)

# reformat result to fit structure of vcfile and save
result_reformat = merge(vcfile[,1:6], result[,c(2:length(colnames(result)))], how="left")

write.table(result_reformat, "mph_outputs/", anc, "_", pheno, "_mph_", vcfile_suff, "_recalculated.mq.vc.csv", 
            quote=F, row.names=F, sep='\t')