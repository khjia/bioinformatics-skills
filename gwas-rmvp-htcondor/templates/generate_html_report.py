#!/usr/bin/env python3
"""
generate_html_report.py — single-file HTML GWAS report.

Reads:
  post_gwas/summary/per_trait_model_summary.tsv
  post_gwas/summary/high_confidence_loci.tsv     (optional)
  post_gwas/summary/model_agreement_matrix.tsv   (optional)
  post_gwas/plots/*.png                          (embedded as base64)

Writes:
  post_gwas/gwas_report.html

No template engine — pure stdlib + base64 embedding so the HTML is portable
(can be emailed / opened anywhere with a browser, no extra files needed).

Usage:
    python generate_html_report.py --base /abs/path/to/05.gwas
"""
from pathlib import Path
import argparse
import base64
import datetime as dt
import html
import pandas as pd


CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       max-width: 1200px; margin: 24px auto; padding: 0 20px; color: #222; }
h1 { border-bottom: 2px solid #2f5597; padding-bottom: 6px; }
h2 { color: #2f5597; margin-top: 32px; border-left: 4px solid #2f5597; padding-left: 8px; }
h3 { color: #555; margin-top: 24px; }
table { border-collapse: collapse; margin: 12px 0; font-size: 13px; }
th, td { border: 1px solid #ddd; padding: 4px 8px; text-align: left; }
th { background: #f4f6fa; }
tr:nth-child(even) { background: #fafbfc; }
.kv { margin: 4px 0; }
.kv b { display: inline-block; min-width: 180px; color: #555; }
.warn { color: #b95000; }
.ok { color: #1b7a1b; }
.bad { color: #b00020; }
img { max-width: 100%; border: 1px solid #eee; margin: 6px 0; }
.toc a { display: block; padding: 2px 0; color: #2f5597; text-decoration: none; }
.toc a:hover { text-decoration: underline; }
details { margin: 8px 0; }
summary { cursor: pointer; font-weight: 500; color: #2f5597; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
"""


def b64img(path: Path) -> str:
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" alt="{path.name}">'


def df_to_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df is None or len(df) == 0:
        return "<p><i>(empty)</i></p>"
    truncated = ""
    if len(df) > max_rows:
        truncated = f'<p class="warn">Showing top {max_rows} of {len(df)} rows.</p>'
        df = df.head(max_rows)
    return truncated + df.to_html(index=False, escape=True, classes="data")


def lambda_class(v: float) -> str:
    if pd.isna(v):
        return ""
    if 0.85 <= v <= 1.15:
        return "ok"
    if 0.7 <= v <= 1.3:
        return "warn"
    return "bad"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", default=None,
                    help="output HTML path (default: <base>/post_gwas/gwas_report.html)")
    ap.add_argument("--embed-plots", action="store_true",
                    help="embed PNG plots inline as base64 (larger file, fully portable)")
    ap.add_argument("--max-traits-with-plots", type=int, default=10,
                    help="limit how many traits get all 4-model plots embedded (default: 10)")
    args = ap.parse_args()

    base = Path(args.base)
    summ_dir = base / "post_gwas" / "summary"
    plot_dir = base / "post_gwas" / "plots"
    out_path = Path(args.out) if args.out else base / "post_gwas" / "gwas_report.html"

    summary_tsv = summ_dir / "per_trait_model_summary.tsv"
    if not summary_tsv.exists():
        raise SystemExit(f"missing {summary_tsv} — run plot_rmvp_all.py first")
    summary = pd.read_csv(summary_tsv, sep="\t")

    loci_tsv = summ_dir / "high_confidence_loci.tsv"
    loci = pd.read_csv(loci_tsv, sep="\t") if loci_tsv.exists() else pd.DataFrame()

    agree_tsv = summ_dir / "model_agreement_matrix.tsv"
    agree = pd.read_csv(agree_tsv, sep="\t") if agree_tsv.exists() else pd.DataFrame()

    bonf_tsv = summ_dir / "bonferroni_significant_hits.tsv"
    bonf_hits = pd.read_csv(bonf_tsv, sep="\t") if bonf_tsv.exists() else pd.DataFrame()

    # ---- header / overview --------------------------------------------------
    n_traits = summary["trait"].nunique()
    models = sorted(summary["model"].unique())
    lam_med = summary.groupby("model")["lambda_gc"].median().to_dict()
    n_bonf_per_model = summary.groupby("model")["n_bonferroni"].sum().to_dict()

    parts = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
             "<title>GWAS Report</title>", f"<style>{CSS}</style></head><body>"]
    parts.append("<h1>GWAS Report</h1>")
    parts.append(f"<p class='kv'><b>Generated:</b> "
                 f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M')}</p>")
    parts.append(f"<p class='kv'><b>Project:</b> {html.escape(str(base))}</p>")
    parts.append(f"<p class='kv'><b>Traits analyzed:</b> {n_traits}</p>")
    parts.append(f"<p class='kv'><b>Models:</b> {', '.join(models)}</p>")

    # TOC
    parts.append("<h2>Contents</h2><div class='toc'>")
    parts.append("<a href='#overview'>1. Overview</a>")
    parts.append("<a href='#calibration'>2. Model calibration (lambda_GC)</a>")
    parts.append("<a href='#high-confidence'>3. High-confidence loci</a>")
    parts.append("<a href='#per-trait'>4. Per-trait summary</a>")
    parts.append("<a href='#multitrait'>5. Multi-trait visualizations</a>")
    parts.append("<a href='#detail'>6. Per-trait plots (top traits)</a>")
    parts.append("</div>")

    # ---- overview -----------------------------------------------------------
    parts.append("<h2 id='overview'>1. Overview</h2>")
    overview_rows = []
    for m in models:
        sub = summary[summary["model"] == m]
        overview_rows.append({
            "model": m,
            "median_lambda_gc": round(float(sub["lambda_gc"].median()), 3),
            "mean_lambda_gc": round(float(sub["lambda_gc"].mean()), 3),
            "traits_in_band_0.85_1.15": int(((sub["lambda_gc"] >= 0.85) &
                                              (sub["lambda_gc"] <= 1.15)).sum()),
            "total_bonferroni_hits": int(sub["n_bonferroni"].sum()),
            "total_suggestive_hits": int(sub["n_suggestive"].sum()),
            "traits_with_any_hit": int((sub["n_bonferroni"] > 0).sum()),
        })
    parts.append(df_to_html(pd.DataFrame(overview_rows), max_rows=10))

    # ---- calibration --------------------------------------------------------
    parts.append("<h2 id='calibration'>2. Model calibration (lambda_GC)</h2>")
    parts.append("<p>Well-calibrated: 0.85 &le; &lambda;<sub>GC</sub> &le; 1.15. "
                 "&lambda; &gt; 1.2 indicates inflation (residual structure); "
                 "&lambda; &lt; 0.85 over-conservative.</p>")

    cal = summary.pivot_table(
        index="trait", columns="model", values="lambda_gc", aggfunc="first")
    rows = ["<table><thead><tr><th>trait</th>"]
    for m in cal.columns:
        rows.append(f"<th>{html.escape(m)}</th>")
    rows.append("</tr></thead><tbody>")
    for trait, row in cal.iterrows():
        rows.append(f"<tr><td>{html.escape(str(trait))}</td>")
        for m in cal.columns:
            v = row[m]
            cls = lambda_class(v)
            txt = "" if pd.isna(v) else f"{v:.3f}"
            rows.append(f"<td class='{cls}'>{txt}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    parts.append("".join(rows))

    # ---- high-confidence loci -----------------------------------------------
    parts.append("<h2 id='high-confidence'>3. High-confidence loci</h2>")
    if len(loci) == 0:
        parts.append("<p><i>No high-confidence loci file found "
                     "(run extract_high_confidence.py).</i></p>")
    else:
        parts.append(f"<p>{len(loci)} loci where >=2 models agree at the chosen "
                     "significance threshold (after windowed clustering).</p>")
        parts.append(df_to_html(loci.sort_values("lead_min_p"), max_rows=50))

    if len(agree) > 0:
        parts.append("<h3>Model agreement matrix</h3>")
        parts.append("<p>Counts of SNPs significant in exactly k models, per trait.</p>")
        parts.append(df_to_html(agree, max_rows=50))

    # ---- per-trait summary --------------------------------------------------
    parts.append("<h2 id='per-trait'>4. Per-trait summary</h2>")
    parts.append(df_to_html(summary, max_rows=200))

    # ---- multitrait plots ---------------------------------------------------
    parts.append("<h2 id='multitrait'>5. Multi-trait visualizations</h2>")
    candidates = ["multitrait_manhattan.FarmCPU.png",
                  "multitrait_manhattan.MLM.png",
                  "trait_chrom_hotspot.png",
                  "trait_chrom_hits.png"]
    any_found = False
    for name in candidates:
        p = plot_dir / name
        if not p.exists():
            continue
        any_found = True
        parts.append(f"<h3>{html.escape(name)}</h3>")
        if args.embed_plots:
            parts.append(b64img(p))
        else:
            rel = p.relative_to(out_path.parent) if out_path.parent in p.parents \
                else p
            parts.append(f"<img src='{html.escape(str(rel))}' alt='{name}'>")
    if not any_found:
        parts.append("<p><i>No multitrait plots found "
                     "(run plot_multitrait_summary.py).</i></p>")

    # ---- per-trait plots (limited) ------------------------------------------
    parts.append("<h2 id='detail'>6. Per-trait plots</h2>")
    # rank traits by smallest min_p across all models
    rank = (summary.groupby("trait")["min_p"].min()
            .sort_values().head(args.max_traits_with_plots).index.tolist())
    parts.append(f"<p>Showing top {len(rank)} traits by minimum p-value across models.</p>")
    for trait in rank:
        parts.append(f"<details><summary>{html.escape(trait)}</summary>")
        parts.append("<div class='grid'>")
        for m in models:
            for kind in ("manhattan", "qq"):
                p = plot_dir / f"{trait}.{m}.{kind}.png"
                if not p.exists():
                    continue
                if args.embed_plots:
                    parts.append(b64img(p))
                else:
                    rel = p.relative_to(out_path.parent) if out_path.parent in p.parents \
                        else p
                    parts.append(f"<img src='{html.escape(str(rel))}' alt='{p.name}'>")
        parts.append("</div></details>")

    parts.append(f"<hr><p><small>Generated by gwas-rmvp-htcondor skill, "
                 f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M')}</small></p>")
    parts.append("</body></html>")

    out_path.write_text("\n".join(parts), encoding="utf-8")
    size_mb = out_path.stat().st_size / 1e6
    print(f"wrote {out_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
