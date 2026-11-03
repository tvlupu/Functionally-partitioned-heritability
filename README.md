# Functionally-partitioned-heritability
This repository contains analysis and plotting code for the paper: Quantifying heritability in functional regions of the coding and non-coding genome across diverse individuals, available as a preprint here: https://www.biorxiv.org/content/10.64898/2026.10.30.763088v1.

Code is organised into the following sub-directories.

## Cohort filtering and PC calculation: `prep/cohort/`
`cohort_filtering.py` creates three genetic similarity groups based on AoU labels and performs kinship filtering using maximum independent set methods \
`calc_pcs.R`, `pcs_umap.R`, and `run_scat_PCs_aou.wdl` calculate principal components and generates UMAP plots per genetic similarity group \

## Generating LD scores: `prep/ldscores/`
`run_scat_LD_aou.wdl` generates variant LD scores per genetic similarity group \
`view_LD.R` graphs and views LD scores per genetic similarity group and/or functional annotation group \

## Other covariates and phenotype processing: `prep/covs_phenos/`

## Functional annotations:  `prep/annots/`
`make_annots.wdl` simplifies annotations to cod/prox/reg/other with `make_annot_files.py`. Uploads updated annotation file to gcloud bucket. (no longer used) \
`make_annots_vep.sh` slade script to make per-chromosome annotation files (coding, RNA, proximal UTR, proximal not UTR, coding intronic, regulatory) \
`run_make_snplists_aou.wdl` using updated annotation file, makes 28 .snplist files (lists of SNPs per bin) based on annotation, MAF, LD per ancestry. "other" snplist created by exclusion \

## Generating GRMs: `aou_h2_workflows/`
`run_make_grm_charAnnots.wdl` makes 28 GCTA GRMs using the 28 snplists \
`run_make_grm_mph.wdl` makes 28 MPH GRMs using the 28 snplists \

## h2SNP and enrichment calculations: `aou_h2_workflows/`
`run_greml_aou.wdl` runs GREML-LDMS \
`run_mph_aou.wdl` runs MPH. Uses `mph_convert_gcta.py` to convert GCTA GRMs to MPH format \
`run_he_aou.wdl` runs HE-reg with GCTA \

## Figure generation and aggregating heritability per bin: `results_processing/`
