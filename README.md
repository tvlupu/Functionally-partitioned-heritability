# Functionally-partitioned-heritability
This repository contains analysis and plotting code for the manuscript *Quantifying heritability in functional regions of the coding and non-coding genome across diverse individuals*, [available as a preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2026.10.30.763088v1).


## Data preparation: `prep/`
### Cohort filtering and covariates
- `cohort_filtering.py` separates individuals into three genetic similarity groups based on AoU labels and performs kinship filtering using maximum independent set methods 
- `calc_pcs.R` / `run_scat_PCs_aou.wdl` / `pcs_umap.R` calculate principal components and generates UMAP plots per genetic similarity group
- `mph_covs_phenos_prep.R` generates covariates files (genotype PCs, sex, age, age squared, and sample type) and rank-inverse normalises + residualises quantitative phenotypes onto the covariates 

### Generating LD scores
- `run_scat_LD_aou.wdl` generates variant LD scores per genetic similarity group 
- `view_LD.R` graphs and views LD scores per genetic similarity group and/or functional annotation group 

### Functional annotations
- `make_annots_vep.sh` slade script to make per-chromosome functional annotation files (coding, ncRNA, proximal UTR, proximal non-UTR, intronic, distal regulatory)
- `run_make_snplists_aou.wdl` using updated annotation file, makes 32 .varlist files (lists of variants per bin) based on annotation, MAF, LD per ancestry. "other" variant list created by exclusion 

### Generating GRMs
- `run_make_grm_charAnnots.wdl` makes 32 GRMs using the 28 autosomal + 4 X-chromosome variant lists. Code uses GCTA-alpha (see **Additional**)

## Data analysis: `aou_h2_workflows/`
### Heritability and fold-enrichment estimation
- `run_mph_aou.wdl` runs [MPH](https://github.com/jiang18/mph/tree/main). Uses `mph_convert_gcta.py` to convert GCTA GRMs to MPH format
- `process_mph_outputs.Rmd` aggregates heritability and fold-enrichment per functional annotation category and for rare versus common categories

### Simulations
- `simulate_mph_polars.py` simulates phenotypes by assigning a set heritability per annotation bin. Assumes a normal distribution of effect sizes within each bin, and uses PLINK2 scoring to scale per-SNP betas. Code also rank-inverse normalises + residualises phenotypes
- phenotypes generated in this way can be used with `run_mph_aou.wdl`
  
## Additional: `aux/`
- `all_plotting_code.R` contains code for plotting main and supplementary figures
- `data/h2_results/` contains MPH variance component estimate output files which can be used with `process_mph_outputs.Rmd` and `all_plotting_code.R`
- `gcta64_alphaGRM.zip` binary file for GCTA-alpha, a modified version of [GCTA](https://github.com/jianyangqt/gcta) which includes an optional `--alpha` flag for GRM calculation. Alpha is used to down-weight the contribution of rare variants to GRMs, see [Speed et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28530675/)
