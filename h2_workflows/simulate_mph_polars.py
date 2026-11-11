#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path
from statistics import NormalDist

import numpy as np
import polars as pl


NORMAL = NormalDist()


# -------------------------------------------------------------------------
# Phenotype Normalisation & Residualization Helpers
# -------------------------------------------------------------------------

def rank_inverse_normal(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    mask = np.isfinite(values)
    if not mask.any():
        return out

    x = values[mask]
    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    ranks_sorted = np.arange(1, len(sorted_x) + 1, dtype=float)

    unique_vals, starts, counts = np.unique(sorted_x, return_index=True, return_counts=True)
    for start, count in zip(starts, counts):
        if count > 1:
            avg_rank = ranks_sorted[start : start + count].mean()
            ranks_sorted[start : start + count] = avg_rank

    ranks = np.empty_like(ranks_sorted)
    ranks[order] = ranks_sorted

    n = len(x)
    probs = (ranks - 0.375) / (n + 0.25)
    out[mask] = np.array([NORMAL.inv_cdf(p) for p in probs], dtype=float)
    return out


def build_design_matrix(df: pl.DataFrame, pc_max: int = 120) -> np.ndarray:
    n = df.height
    parts = [np.ones((n, 1), dtype=float)]

    for cat_col in ["sex", "sample_type"]:
        if cat_col in df.columns:
            dummy_df = df.select(cat_col).to_dummies(columns=[cat_col])
            if len(dummy_df.columns) > 1:
                for col in dummy_df.columns[1:]:
                    parts.append(dummy_df[col].cast(pl.Float64).to_numpy().reshape(-1, 1))

    numeric_cols = [c for c in ["age", "age_sq"] if c in df.columns]
    numeric_cols.extend([f"pc{i}" for i in range(1, pc_max + 1) if f"pc{i}" in df.columns])

    for col in numeric_cols:
        parts.append(df[col].cast(pl.Float64).to_numpy().reshape(-1, 1))

    return np.hstack(parts)


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    mask = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    out = np.full(y.shape, np.nan, dtype=float)
    if not mask.any():
        return out
    beta, *_ = np.linalg.lstsq(x[mask], y[mask], rcond=None)
    out[mask] = y[mask] - x[mask] @ beta
    return out


# -------------------------------------------------------------------------
# Dynamic Binning and Calibration Core
# -------------------------------------------------------------------------

def pick_score_column(df: pl.DataFrame) -> str:
    preferred = ["gKpost", "gKpre", "SCORE1_SUM", "SCORE1_AVG"]
    for col in preferred:
        if col in df.columns:
            return col
    return df.columns[-1]


def get_or_generate_afreq(plink2_path: Path, pfile: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "snps_freqs"
    afreq_path = out_prefix.with_suffix(".afreq")
    if afreq_path.exists():
        print(f"Using existing afreq file: {afreq_path}")
        return afreq_path

    print("Generating allele frequencies (plink2 --freq)...")
    cmd = [
        str(plink2_path), "--pfile", pfile,
        "--freq", "cols=+reffreq,+alt1freq,-maybeprovref",
        "--nonfounders", "--out", str(out_prefix)
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return afreq_path


def read_ld_scores(ld_path: Path) -> pl.DataFrame:
    print(f"Reading LD scores from: {ld_path}")
    return pl.read_csv(
        ld_path, has_header=False, separator=" ",
        new_columns=["SNP", "CHR", "MAF_LD", "ldscore"],
        ignore_errors=True)


def read_annotations(annot_path: Path) -> pl.DataFrame:
    print(f"Reading annotations from: {annot_path}")
    return pl.read_csv(
        annot_path, has_header=False, separator=" ",
        #new_columns=["SNP", "ANNOT"], schema_overrides={"SNP": pl.Utf8, "ANNOT": pl.Utf8},
        new_columns=["SNP"],
        ignore_errors=True
    ).unique(subset=["SNP"], keep="first")


def calibrate_bin(
    args: argparse.Namespace, rng: np.random.Generator, afreq_path: Path
) -> Path:
    ld_df = read_ld_scores(args.ld_file)
    annot_df = read_annotations(args.annot_file)


    # 1. Join LD & Annotations
    df = ld_df.join(annot_df, on="SNP", how="right")
    #df = df.with_columns(pl.col("ANNOT").fill_null("other"))

    # 2. Subset to relevant annotation and LD MAF (as per original logic)
    #df = df.filter(
        #(pl.col("ANNOT") == args.annot) &
        #(pl.col("MAF_LD") > args.maf_low) &
        #(pl.col("MAF_LD") <= args.maf_high)
    #)

    # 3. LD Bin split
    #if args.ld_bin in {"lowbin", "highbin"}:
        #q2 = df.select(pl.col("ldscore").median()).item()
        #if args.ld_bin == "lowbin":
            #df = df.filter(pl.col("ldscore") < q2)
        #else:
            #df = df.filter(pl.col("ldscore") >= q2)
    
    print(f"Num SNPs after filtering MAF & Annotation & LD bin: {df.height}")
    if df.height == 0:
        raise ValueError(f"No SNPs remaining for bin {args.bin_name}.")

    # 4. Join accurate PLINK Freq (for calibration math)
    afreq_df = pl.read_csv(afreq_path, separator="\t", has_header=True)
    afreq_df = afreq_df.with_columns([
        pl.col("ALT1_FREQ").cast(pl.Float64),
        pl.when(pl.col("ALT1_FREQ") > 0.5).then(pl.col("REF")).otherwise(pl.col("ALT")).alias("ALLELE"),
        pl.when(pl.col("ALT1_FREQ") > 0.5).then(1.0 - pl.col("ALT1_FREQ")).otherwise(pl.col("ALT1_FREQ")).alias("EMPIRICAL_MAF")
    ]).rename({"ID": "SNP"}).unique(subset=["SNP"], maintain_order=True)

    df = df.join(afreq_df, on="SNP", how="inner").unique(subset=["SNP"], maintain_order=True)

    # 5. Generate Raw Beta with Var(beta) scaling -> z_j * [2p(1-p)]^(alpha/2)
    #z = rng.normal(0, 1, df.height)
    z = rng.gamma(shape=0.186, scale=0.05516129032, size=df.height)
    # randomly flip sign for exactly half the entries (floor if odd) for symmetrical effect size distribution
    idx = rng.permutation(df.height)[: df.height // 2]
    z[idx] *= -1

    p_col = pl.col('EMPIRICAL_MAF').clip(1e-6, 0.5)
    
    df = df.with_columns([
    (pl.lit(z)* (2.0 * p_col * (1.0 - p_col)) ** (args.alpha / 2.0)).alias("BETA_RAW")
])
    
    args.out_dir.mkdir(parents=True, exist_ok=True)
    
    # Save unscaled raw beta for PLINK scoring
    raw_weights = args.out_dir / f"{args.bin_name}_raw_weights.tsv"
    df.select(["SNP", "ALLELE", "BETA_RAW"]).write_csv(raw_weights, separator="\t", include_header=False)
    
    # 6. Score with Plink2
    print(f"Scoring with PLINK (raw effects) for {args.bin_name}...")
    out_prefix = args.out_dir / f"{args.bin_name}_raw"
    cmd = [
        str(args.plink2), "--pfile", args.pfile, "--read-freq", "sim_outputs/snps_freqs.afreq",
        "--score", str(raw_weights), "center", 
        "--score-col-nums", "3", "--out", str(out_prefix), "--extract", args.annot_file
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    
    # 7. Heritability Calibration Post-Score
    scores = pl.read_csv(f"{out_prefix}.sscore", separator="\t", has_header=True)
    score_col = pick_score_column(scores)

    # Use ddof=0 for strictly large population limit tracking
    var_raw = scores.select(pl.col(score_col).var(ddof=0)).item() 
    c_k = math.sqrt(args.hk / var_raw)
    
    print(f"Bin {args.bin_name} raw Var(G_k*) = {var_raw}. Calibration Factor (C_k) = {c_k}")

    # Scale Betas & Scores
    df.with_columns(
        (pl.col("BETA_RAW") * c_k).alias("BETA_CALIBRATED")
    ).select(["SNP", "ALLELE", "BETA_CALIBRATED"]).write_csv(
        args.out_dir / f"{args.bin_name}_calibrated_weights.tsv", separator="\t", include_header=False
    )
    
    keep_cols = [c for c in ["#FID", "FID", "IID"] if c in scores.columns]
    cal_scores = scores.with_columns((pl.col(score_col) * c_k).alias("gKpost")).select(keep_cols + ["gKpost"])
    cal_path = args.out_dir / f"{args.bin_name}_calibrated_scores.tsv"
    cal_scores.write_csv(cal_path, separator="\t")
    
    return cal_path


# -------------------------------------------------------------------------
# Step 3: Phenotype Combination
# -------------------------------------------------------------------------

def combine_and_residualize(
    score_paths: list[Path], 
    covs_path: Path, 
    sim_pheno_name: str, 
    anc: str, 
    h_tot: float, 
    out_dir: Path, 
    rng: np.random.Generator
) -> None:
    # 1. Sum up all provided gKpost bins
    frames: list[pl.DataFrame] = []
    for path in score_paths:
        df = pl.read_csv(path, separator="\t", has_header=True)
        col = pick_score_column(df)
        keep = ["IID", col]
        if "FID" in df.columns: keep = ["FID", "IID", col]
        elif "#FID" in df.columns: keep = ["#FID", "IID", col]
        frames.append(df.select(keep).rename({col: path.stem}))

    merged = frames[0]
    for df in frames[1:]:
        join_keys = [c for c in ["#FID", "FID", "IID"] if c in merged.columns and c in df.columns]
        merged = merged.join(df, on=join_keys, how="inner")

    score_cols = [c for c in merged.columns if c not in {"#FID", "FID", "IID"}]
    merged = merged.with_columns(pl.sum_horizontal(score_cols).alias("gKpost_total"))

    # 2. Add residual noise scaled by sum total non-genetic variance
    noise_sd = math.sqrt(max(0.0, 1.0 - h_tot))
    raw = merged.select(["IID", "gKpost_total"]).with_columns(
        (pl.col("gKpost_total") + pl.Series("noise", rng.normal(0.0, noise_sd, merged.height))).alias("pheno")
    ).select(["IID", "pheno"])

    out_dir.mkdir(parents=True, exist_ok=True)
    raw.write_csv(out_dir / f"{sim_pheno_name}_{anc}_raw.csv")

    # 3. Residualise & Rank Inverse Normalize
    print("Residualising against covariates & Rank-Inverse Normalising...")
    covs = pl.read_csv(covs_path, separator="\t", has_header=True, try_parse_dates=False)
    pheno_covs = covs.join(raw, on="IID", how="inner").sort("IID")

    sex = pheno_covs["sex"].to_numpy()
    pheno = pheno_covs["pheno"].to_numpy()
    rint = np.full(pheno.shape, np.nan, dtype=float)

    for sex_value in pheno_covs.select("sex").drop_nulls().unique().to_series().to_list():
        mask = sex == sex_value
        rint[mask] = rank_inverse_normal(pheno[mask])

    pheno_covs = pheno_covs.with_columns(pl.Series("rint_pheno", rint))
    x = build_design_matrix(pheno_covs, pc_max=120)
    resid = residualize(pheno_covs["rint_pheno"].to_numpy(), x)

    mph_out = pl.DataFrame({"IID": pheno_covs["IID"], sim_pheno_name: resid})
    mph_out.write_csv(out_dir / f"{sim_pheno_name}_{anc}_mph.csv")
    print(f"Final output written to {out_dir}/{sim_pheno_name}_{anc}_mph.csv")


def parse_args():
    p = argparse.ArgumentParser(description="Iterative heritability bin simulator")
    
    # Modes
    p.add_argument("--skip-pheno", action="store_true", help="Just generate bin scores, skip combining phenotype.")
    p.add_argument("--only-pheno", action="store_true", help="Skip bin generation, only aggregate listed score files.")
    
    # Setup
    p.add_argument("--bin-name", required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-dir", type=Path, default=Path("sim_outputs"))
    
    # Bin Gen Limits
    p.add_argument("--pfile", help="PLINK2 pfile prefix")
    p.add_argument("--plink2", type=Path, default=Path("./plink2"))
    p.add_argument("--ld-file", type=Path)
    p.add_argument("--annot", help="Annotation snplist prefix")
    p.add_argument("--annot-file", type=Path)
    p.add_argument("--maf-low", type=float)
    p.add_argument("--maf-high", type=float)
    p.add_argument("--ld-bin", choices=["lowbin", "highbin", "all"], default="all")
    p.add_argument("--alpha", type=float, default=-0.25)
    p.add_argument("--hk", type=float, help="Heritability slice matching this specific bin.")
    
    # Phenotype Limits
    p.add_argument("--covs", type=Path)
    p.add_argument("--anc")
    p.add_argument("--sim-pheno-name")
    p.add_argument("--htot", type=float, help="Total Heritability target (for noise scaling).")
    p.add_argument("--extra-score-files", type=Path, nargs="*", default=[])
    
    return p.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    calibrated_path: Path | None = None

    if not args.only_pheno:
        afreq = get_or_generate_afreq(args.plink2, args.pfile, args.out_dir)
        calibrated_path = calibrate_bin(args, rng, afreq)
        
    if not args.skip_pheno:
        files = args.extra_score_files.copy()
        if calibrated_path:
            files.insert(0, calibrated_path)
            
        combine_and_residualize(
            files, args.covs, args.sim_pheno_name, args.anc, 
            args.htot, args.out_dir, rng
        )

if __name__ == "__main__":
    main()
