#!/usr/bin/env python3
"""
extract_high_confidence.py — multi-model consensus SNP integration.

For each trait, find SNPs that are significant in >=N models (default >=2).
A "high-confidence locus" is the union of nearby consensus SNPs (LD-style
windowing on physical distance, default 100kb).

Inputs:  final_results/trait_XX.{GLM,MLM,FarmCPU,BLINK}.csv
Outputs: post_gwas/summary/
  - high_confidence_snps.tsv         per-SNP rows with model_count, models_hit, min_p
  - high_confidence_loci.tsv         clustered loci (chr, start, end, lead_snp, n_snps)
  - model_agreement_matrix.tsv       trait x (1-model, 2-model, 3-model, 4-model) counts
"""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import json

MODELS_DEFAULT = ["GLM", "MLM", "FarmCPU", "BLINK"]


def load_pcounts(map_file: Path) -> int:
    with map_file.open() as f:
        return sum(1 for _ in f) - 1


def read_one(path: Path):
    df = pd.read_csv(path)
    pcol = df.columns[-1]
    df = df.rename(columns={pcol: "p", "CHROM": "chr", "POS": "pos"})
    df["chr"] = df["chr"].astype(str)
    df = df[(df["p"] > 0) & (df["p"] <= 1)]
    return df[["SNP", "chr", "pos", "MAF", "Effect", "SE", "p"]]


def cluster_loci(df: pd.DataFrame, window: int) -> pd.DataFrame:
    """Greedy left-to-right clustering by chr+pos within `window` bp."""
    if df.empty:
        return df.assign(locus_id=pd.Series(dtype=int))
    df = df.sort_values(["chr", "pos"]).reset_index(drop=True)
    locus_id = np.zeros(len(df), dtype=int)
    cur = 0
    last_chr, last_pos = None, -10**12
    for i, (c, p) in enumerate(zip(df["chr"], df["pos"])):
        if c != last_chr or p - last_pos > window:
            cur += 1
        locus_id[i] = cur
        last_chr, last_pos = c, p
    df["locus_id"] = locus_id
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="path to <proj>/05.gwas")
    ap.add_argument("--res-dir", default="rMVP/final_results",
                    help="results dir relative to --base (default: rMVP/final_results)")
    ap.add_argument("--models", nargs="+", default=MODELS_DEFAULT)
    ap.add_argument("--min-models", type=int, default=2,
                    help="SNP must be significant in >= this many models (default: 2)")
    ap.add_argument("--threshold", choices=["bonferroni", "suggestive"], default="bonferroni")
    ap.add_argument("--window", type=int, default=100_000,
                    help="bp window for clustering nearby consensus SNPs (default: 100000)")
    args = ap.parse_args()

    base = Path(args.base)
    res_dir = base / args.res_dir
    out = base / "post_gwas" / "summary"
    out.mkdir(parents=True, exist_ok=True)

    m_total = load_pcounts(base / "rMVP" / "mvp.geno.map")
    bonf = 0.05 / m_total
    sugg = 1.0 / m_total
    cutoff = bonf if args.threshold == "bonferroni" else sugg
    print(f"M_total={m_total:,} | cutoff={args.threshold}={cutoff:.3e} | "
          f"min_models={args.min_models} | window={args.window:,}bp")

    # trait_zh map (optional)
    cmap = base / "pheno_category_map.json"
    trait_zh = {}
    if cmap.exists():
        cm = json.loads(cmap.read_text())
        trait_zh = {v: k for k, v in cm.get("original_to_safe", cm).items()}

    # discover traits from CSV listing
    traits = sorted({p.stem.rsplit(".", 1)[0]
                     for p in res_dir.glob("trait_*.csv")
                     if "_signals" not in p.stem})
    print(f"Found {len(traits)} traits")

    consensus_rows = []
    agreement_rows = []

    for trait in traits:
        # collect significant SNP set per model
        per_model = {}
        for model in args.models:
            path = res_dir / f"{trait}.{model}.csv"
            if not path.exists():
                continue
            try:
                df = read_one(path)
            except Exception as e:
                print(f"[ERR] {path.name}: {e}")
                continue
            sig = df[df["p"] <= cutoff].copy()
            per_model[model] = sig

        if not per_model:
            continue

        # pivot: SNP -> set of models that hit it
        all_sig = pd.concat(
            [d.assign(model=m) for m, d in per_model.items()],
            ignore_index=True,
        )
        if all_sig.empty:
            agreement_rows.append({
                "trait": trait, "trait_zh": trait_zh.get(trait, trait),
                **{f"n_{k}_model_hits": 0 for k in range(1, len(args.models) + 1)},
            })
            continue

        # group by SNP
        grp = (all_sig
               .groupby(["SNP", "chr", "pos"], as_index=False)
               .agg(model_count=("model", "nunique"),
                    models_hit=("model", lambda s: ",".join(sorted(set(s)))),
                    min_p=("p", "min"),
                    mean_maf=("MAF", "mean")))

        # agreement counts (how many SNPs were hit by exactly k models)
        ag = {f"n_{k}_model_hits": int((grp["model_count"] == k).sum())
              for k in range(1, len(args.models) + 1)}
        agreement_rows.append({"trait": trait, "trait_zh": trait_zh.get(trait, trait), **ag})

        # high-confidence subset
        hc = grp[grp["model_count"] >= args.min_models].copy()
        if hc.empty:
            continue
        hc.insert(0, "trait", trait)
        hc.insert(1, "trait_zh", trait_zh.get(trait, trait))
        consensus_rows.append(hc)

    if consensus_rows:
        snps = pd.concat(consensus_rows, ignore_index=True).sort_values(
            ["trait", "chr", "pos"])
        snps.to_csv(out / "high_confidence_snps.tsv", sep="\t", index=False)
        print(f"  wrote {len(snps)} high-confidence SNP rows")

        # cluster into loci per trait
        locus_rows = []
        for trait, sub in snps.groupby("trait"):
            sub = cluster_loci(sub.copy(), window=args.window)
            for lid, lsub in sub.groupby("locus_id"):
                lead = lsub.loc[lsub["min_p"].idxmin()]
                locus_rows.append({
                    "trait": trait,
                    "trait_zh": lsub["trait_zh"].iloc[0],
                    "chr": lead["chr"],
                    "start": int(lsub["pos"].min()),
                    "end": int(lsub["pos"].max()),
                    "n_snps": len(lsub),
                    "lead_snp": lead["SNP"],
                    "lead_pos": int(lead["pos"]),
                    "lead_min_p": float(lead["min_p"]),
                    "max_model_count": int(lsub["model_count"].max()),
                    "models_at_lead": lead["models_hit"],
                })
        loci = pd.DataFrame(locus_rows).sort_values(
            ["trait", "chr", "start"])
        loci.to_csv(out / "high_confidence_loci.tsv", sep="\t", index=False)
        print(f"  wrote {len(loci)} high-confidence loci")
    else:
        print("  no high-confidence SNPs at given threshold/min-models")

    pd.DataFrame(agreement_rows).to_csv(
        out / "model_agreement_matrix.tsv", sep="\t", index=False)
    print("DONE")


if __name__ == "__main__":
    main()
