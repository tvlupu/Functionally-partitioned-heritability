library(bigsnpr)

ancestry <- commandArgs(trailingOnly=T)[1]
umap_no <- commandArgs(trailingOnly=T)[2]

# options so bigsnpR doesn't crash
options(bigstatsr.ncores.max = 32)
options(bigstatsr.check.parallel.blas = FALSE)
try(bigparallelr::set_blas_ncores(32))

rds <- snp_readBed2(paste0(ancestry,"_umap", umap_no, ".bed"), paste0(ancestry,"_umap", umap_no), ncores=31) # convert to bigsnpR

full <- snp_attach(paste0(ancestry,"_umap", umap_no, ".rds")) # load bigsnpr file
G <- full$genotypes # define genotype matrix
CHR <- full$map$chromosome 
POS <- full$map$physical.pos

# impute array data
infos <- snp_fastImputeSimple(G, method="mode", ncores = 31) 

# modify .rds to contain imputed genotypes
full$genotypes$code256 <- CODE_IMPUTE_PRED
full <- snp_save(full)

full$code256 <- bigsnpr:::CODE_IMPUTE_PRED    
full$genotypes$code256 <- bigsnpr:::CODE_IMPUTE_PRED

saveRDS(full, paste0(ancestry,"_umap", umap_no, ".rds"))
print('imputation complete')

full <- snp_attach(paste0(ancestry,"_umap", umap_no, ".rds")) # load bigsnpr file
G <- full$genotypes # define genotype matrix
CHR <- full$map$chromosome 
POS <- full$map$physical.pos

# run PCA
svd <- snp_autoSVD(G, CHR, POS, k= 120, min.maf = 0.01, ncores = 31)

# export
write.table(svd$u, paste0("svd_", ancestry, "_umap", umap_no, ".csv"), sep='/t', quote=FALSE, row.names=FALSE)
