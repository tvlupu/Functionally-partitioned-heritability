# Functionally-partitioned-heritability
This repository contains analysis and plotting code for the manuscript *Quantifying heritability in functional regions of the coding and non-coding genome across diverse individuals*, [available as a preprint in bioRxiv](https://www.biorxiv.org/content/10.64898/2026.10.30.763088v1).


## Data preparation: `prep/`
### Cohort filtering and covariates
- `cohort_filtering.py` separates individuals into three genetic similarity groups based on AoU labels and performs kinship filtering using maximum independent set methods 
- `calc_pcs.R` / `pcs_umap.R` / `run_scat_PCs_aou.wdl` calculate principal components and generates UMAP plots per genetic similarity group
- `mph_covs_phenos_prep.R` generates covariates files (genotype PCs, sex, age, age squared, and sample type) and rank-inverse normalises quantitative phenotypes onto the covariates 

### Generating LD scores
- `run_scat_LD_aou.wdl` generates variant LD scores per genetic similarity group 
- `view_LD.R` graphs and views LD scores per genetic similarity group and/or functional annotation group 

### Functional annotations
- `make_annots_vep.sh` slade script to make per-chromosome functional annotation files (coding, ncRNA, proximal UTR, proximal non-UTR, intronic, distal regulatory)
- `run_make_snplists_aou.wdl` using updated annotation file, makes 32 .varlist files (lists of variants per bin) based on annotation, MAF, LD per ancestry. "other" variant list created by exclusion 

### Generating GRMs
- `run_make_grm_charAnnots.wdl` makes 32 GRMs using the 28 autosomal + 4 X-chromosome variant lists. Code uses GCTA-alpha (see **Additional**)

## Data analysis: `aou_h2_workflows/`
### h2SNP and enrichment calculations
- `run_greml_aou.wdl` runs GCTA GREML-LDMS
- `run_mph_aou.wdl` runs [MPH](https://github.com/jiang18/mph/tree/main). Uses `mph_convert_gcta.py` to convert GCTA GRMs to MPH format
- `process_mph_outputs.Rmd` aggregates heritability and enrichment per functional annotation category and for rare versus common categories 

## Additional: `aux/`
- `all_plotting_code.R` contains code for plotting main and supplementary figures
- `data/h2_results/` contains MPH variance component estimate output files which can be used with `process_mph_outputs.Rmd` and `all_plotting_code.R`
- `gcta64_alphaGRM.zip` binary file for GCTA-alpha, a modified version of [GCTA](https://github.com/jianyangqt/gcta) which includes an optional `--alpha` flag for GRM calculation. Alpha is used to down-weight the contribution of rare variants to GRMs, see [Speed et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28530675/)
