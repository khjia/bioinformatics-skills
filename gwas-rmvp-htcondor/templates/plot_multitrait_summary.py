#!/usr/bin/env python3
"""
plot_multitrait_summary.py — multi-trait summary visualizations.

Generates:
  1. multitrait_manhattan.png   stacked Manhattan rows (one per trait), single model
  2. trait_chrom_hotspot.png    heatmap of -log10(min_p) per (trait, chrom) — finds hotspots
  3. trait_chrom_hits.png       heatmap of #significant SNPs per (trait, chrom)

Usage:
    python plot_multitrait_summary.py --base /abs/path/to/05.gwas \
        --model FarmCPU --threshold suggestive

For Manhattan rows we use ONE model (default FarmCPU — most balanced on this project).
The hotspot heatmap aggregates across all models that exist.
"""
from pathlib import Path
import argparse
import math
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODELS_ALL = ["GLM", "MLM", "FarmCPU", "BLINK"]


def load_pcounts(map_file: Path) -> int:
    with map_file.open() as f:
        return sum(1 for _ in f) - 1


def read_one(path: Path):
    df = pd.read_csv(path)
    pcol = df.columns[-1]
    df = df.rename(columns={pcol: "p", "CHROM": "chr", "POS": "pos"})
    df["chr"] = df["chr"].astype(str)
    df = df[(df["p"] > 0) & (df["p"] <= 1)]
    return df[["SNP", "chr", "pos", "p"]]


def stacked_manhattan(traits, res_dir, model, chrom_order, bonf, sugg, out_path,
                      trait_zh):
    n = len(traits)
    fig, axes = plt.subplots(n, 1, figsize=(14, 1.6 * n + 1), dpi=130, sharex=True)
    if n == 1:
        axes = [axes]
    rng = np.random.default_rng(1)
    colors = ["#2f5597", "#9e480e"]

    # global x-axis offsets (use first available trait to get max positions)
    maxpos = {}
    for trait in traits:
        path = res_dir / f"{trait}.{model}.csv"
        if not path.exists():
            continue
        df = read_one(path)
        for c in chrom_order:
            if (df["chr"] == c).any():
                maxpos[c] = max(maxpos.get(c, 0),
                                int(df.loc[df["chr"] == c, "pos"].max()))

    offsets = {}; ticks = []; labels = []; off = 0
    for c in chrom_order:
        if c in maxpos:
            offsets[c] = off
            ticks.append(off + maxpos[c] / 2)
            labels.append(c)
            off += maxpos[c]

    for ax, trait in zip(axes, traits):
        path = res_dir / f"{trait}.{model}.csv"
        if not path.exists():
            ax.set_ylabel(trait, rotation=0, ha="right", va="center", fontsize=8)
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                    ha="center", va="center", color="grey", fontsize=8)
            continue
        df = read_one(path)
        df = df[df["chr"].isin(chrom_order)]
        # downsample non-significant
        sig_mask = df["p"].values <= sugg
        nonsig = np.flatnonzero(~sig_mask)
        if len(nonsig) > 80_000:
            keep = rng.choice(nonsig, size=80_000, replace=False)
            df = df.iloc[np.r_[np.flatnonzero(sig_mask), keep]]
        df = df.copy()
        df["x"] = df["pos"] + df["chr"].map(offsets)
        df["mlogp"] = -np.log10(df["p"])
        for i, c in enumerate(chrom_order):
            sub = df[df["chr"] == c]
            if len(sub):
                ax.scatter(sub["x"], sub["mlogp"],
                           s=1.2, c=colors[i % 2], alpha=.6, linewidths=0)
        ax.axhline(-math.log10(bonf), color="red", lw=.6, ls="--")
        ax.axhline(-math.log10(sugg), color="orange", lw=.5, ls=":")
        label = trait
        ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(axis="y", alpha=.15)

    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels)
    axes[-1].set_xlabel("Chromosome")
    fig.suptitle(f"Multi-trait Manhattan ({model})", y=1.0)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def hotspot_heatmaps(traits, res_dir, models, chrom_order, cutoff,
                     out_minp, out_count):
    """Two heatmaps: min_p (-log10) and #hits per (trait, chrom)."""
    minp_mat = np.full((len(traits), len(chrom_order)), np.nan)
    count_mat = np.zeros((len(traits), len(chrom_order)), dtype=int)

    for ti, trait in enumerate(traits):
        for model in models:
            path = res_dir / f"{trait}.{model}.csv"
            if not path.exists():
                continue
            df = read_one(path)
            for ci, c in enumerate(chrom_order):
                sub = df[df["chr"] == c]
                if sub.empty:
                    continue
                mp = float(sub["p"].min())
                cur = minp_mat[ti, ci]
                if np.isnan(cur) or mp < cur:
                    minp_mat[ti, ci] = mp
                count_mat[ti, ci] += int((sub["p"] <= cutoff).sum())

    # heatmap 1: -log10(min_p)
    mlog = -np.log10(minp_mat)
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(chrom_order) + 3),
                                    max(4, 0.3 * len(traits))), dpi=140)
    im = ax.imshow(mlog, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(chrom_order)))
    ax.set_xticklabels(chrom_order)
    ax.set_yticks(range(len(traits)))
    ax.set_yticklabels(traits, fontsize=7)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Trait")
    ax.set_title("Hotspot heatmap: -log10(min P) across models")
    fig.colorbar(im, ax=ax, label="-log10(min P)")
    fig.tight_layout()
    fig.savefig(out_minp); plt.close(fig)

    # heatmap 2: hit counts
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(chrom_order) + 3),
                                    max(4, 0.3 * len(traits))), dpi=140)
    # log1p for display
    im = ax.imshow(np.log1p(count_mat), aspect="auto", cmap="Blues")
    ax.set_xticks(range(len(chrom_order)))
    ax.set_xticklabels(chrom_order)
    ax.set_yticks(range(len(traits)))
    ax.set_yticklabels(traits, fontsize=7)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Trait")
    ax.set_title("Significant SNP counts per (trait, chrom) [log1p scaled]")
    fig.colorbar(im, ax=ax, label="log1p(count)")
    # annotate cells with raw count if small
    for ti in range(len(traits)):
        for ci in range(len(chrom_order)):
            v = count_mat[ti, ci]
            if v > 0:
                ax.text(ci, ti, str(v), ha="center", va="center",
                        fontsize=6, color="black" if v < count_mat.max() / 2 else "white")
    fig.tight_layout()
    fig.savefig(out_count); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--res-dir", default="rMVP/final_results")
    ap.add_argument("--model", default="FarmCPU",
                    help="model used for stacked Manhattan (default: FarmCPU)")
    ap.add_argument("--models-for-hotspot", nargs="+", default=MODELS_ALL)
    ap.add_argument("--threshold", choices=["bonferroni", "suggestive"],
                    default="suggestive",
                    help="cutoff for hit-count heatmap (default: suggestive)")
    ap.add_argument("--chroms", nargs="+",
                    default=[str(i) for i in range(1, 8)],
                    help="chromosome list in plotting order")
    args = ap.parse_args()

    base = Path(args.base)
    res_dir = base / args.res_dir
    out = base / "post_gwas" / "plots"
    out.mkdir(parents=True, exist_ok=True)

    m_total = load_pcounts(base / "rMVP" / "mvp.geno.map")
    bonf = 0.05 / m_total
    sugg = 1.0 / m_total
    cutoff = bonf if args.threshold == "bonferroni" else sugg

    cmap = base / "pheno_category_map.json"
    trait_zh = {}
    if cmap.exists():
        cm = json.loads(cmap.read_text())
        trait_zh = {v: k for k, v in cm.get("original_to_safe", cm).items()}

    traits = sorted({p.stem.rsplit(".", 1)[0]
                     for p in res_dir.glob("trait_*.csv")
                     if "_signals" not in p.stem})
    print(f"traits={len(traits)} | model={args.model} | cutoff={cutoff:.3e}")

    stacked_manhattan(traits, res_dir, args.model, args.chroms,
                      bonf, sugg, out / f"multitrait_manhattan.{args.model}.png",
                      trait_zh)
    print(f"  wrote multitrait_manhattan.{args.model}.png")

    hotspot_heatmaps(traits, res_dir, args.models_for_hotspot, args.chroms,
                     cutoff,
                     out / "trait_chrom_hotspot.png",
                     out / "trait_chrom_hits.png")
    print("  wrote trait_chrom_hotspot.png + trait_chrom_hits.png")


if __name__ == "__main__":
    main()
