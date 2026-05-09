#!/usr/bin/env python3
"""
rMVP post-processing: Manhattan + QQ for all traits x 3 models + top-hits summary.
Reads: <proj>/05.gwas/rMVP/results/trait_XX.{GLM,MLM,FarmCPU}.csv
Writes: post_gwas/plots/*.png and post_gwas/summary/*.tsv

EDIT the BASE path and TRAITS list to match your project.
"""
from pathlib import Path
import math, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
from multiprocessing import Pool

BASE    = Path('/abs/path/to/project/05.gwas')   # EDIT
# If adaptive PC tuning was used, point to final_results/ (hardlinked winners).
# If only the single PC=3 run was done, point to rMVP/results/ instead.
RES_DIR = BASE / 'rMVP' / 'final_results'
OUT     = BASE / 'post_gwas'
PLOT    = OUT / 'plots'
SUMM    = OUT / 'summary'
PLOT.mkdir(parents=True, exist_ok=True)
SUMM.mkdir(parents=True, exist_ok=True)

MAP_FILE = BASE / 'rMVP' / 'mvp.geno.map'
with MAP_FILE.open() as f:
    M_TOTAL = sum(1 for _ in f) - 1
BONF = 0.05 / M_TOTAL
SUGG = 1.0 / M_TOTAL
print(f'M_total={M_TOTAL:,}  Bonferroni={BONF:.3e}  Suggestive={SUGG:.3e}')

# optional: Chinese / alternate name map
cmap_path = BASE / 'pheno_category_map.json'
if cmap_path.exists():
    CMAP = json.loads(cmap_path.read_text())
    TRAIT_ZH = {v: k for k, v in CMAP.get('original_to_safe', {}).items()}
else:
    TRAIT_ZH = {}

MODELS = ['GLM', 'MLM', 'FarmCPU', 'BLINK']
# EDIT: list of traits that were actually run (exclude any skipped with n_valid<50)
TRAITS = [f'trait_{i:02d}' for i in range(1, 30) if i != 3]

# chromosomes to plot (x-axis). For plants commonly 1..N.
CHROM_ORDER = [str(i) for i in range(1, 8)]


def lambda_gc(p):
    if len(p) == 0:
        return float('nan')
    return float(np.nanmedian(chi2.isf(p, 1)) / 0.454936423119572)


def read_rmvp(path):
    # last column is the p-value (named trait_XX.MODEL).
    df = pd.read_csv(path)
    pcol = df.columns[-1]
    df = df.rename(columns={pcol: 'p', 'CHROM': 'chr', 'POS': 'pos'})
    df['chr'] = df['chr'].astype(str)
    df = df[df['chr'].isin(CHROM_ORDER)]
    df = df[(df['p'] > 0) & (df['p'] <= 1)]
    return df[['SNP', 'chr', 'pos', 'MAF', 'Effect', 'SE', 'p']].reset_index(drop=True)


def plot_pair(df, title, outprefix):
    if len(df) == 0:
        return 0, float('nan'), float('nan')
    rng = np.random.default_rng(1)
    sig_mask = df['p'].values <= SUGG
    nonsig = np.flatnonzero(~sig_mask)
    if len(nonsig) > 350_000:
        sampled = rng.choice(nonsig, size=350_000, replace=False)
        idx = np.r_[np.flatnonzero(sig_mask), sampled]
        dplot = df.iloc[idx].copy()
    else:
        dplot = df.copy()
    maxpos = {c: int(df.loc[df['chr'] == c, 'pos'].max()) for c in CHROM_ORDER if (df['chr'] == c).any()}
    offsets = {}; ticks = []; labels = []; off = 0
    for c in CHROM_ORDER:
        if c in maxpos:
            offsets[c] = off; ticks.append(off + maxpos[c] / 2); labels.append(c); off += maxpos[c]
    dplot['x'] = dplot['pos'] + dplot['chr'].map(offsets)
    dplot['mlogp'] = -np.log10(dplot['p'])
    lam = lambda_gc(df['p'].values)
    minp = float(df['p'].min())

    # Manhattan
    fig, ax = plt.subplots(figsize=(14, 5), dpi=140)
    colors = ['#2f5597', '#9e480e']
    for i, c in enumerate(CHROM_ORDER):
        sub = dplot[dplot['chr'] == c]
        if len(sub):
            ax.scatter(sub['x'], sub['mlogp'], s=1.5, c=colors[i % 2], alpha=.65, linewidths=0)
    ax.axhline(-math.log10(BONF), color='red',    lw=1, ls='--', label=f'Bonferroni ({BONF:.2e})')
    ax.axhline(-math.log10(SUGG), color='orange', lw=1, ls=':',  label=f'Suggestive ({SUGG:.2e})')
    ax.set_xticks(ticks); ax.set_xticklabels(labels)
    ax.set_xlabel('Chromosome'); ax.set_ylabel('-log10(P)')
    ax.set_title(f'{title}\nminP={minp:.2e}  lambdaGC={lam:.3f}  n={len(df):,}')
    ax.legend(fontsize=8, loc='upper right'); ax.grid(axis='y', alpha=.2)
    fig.tight_layout()
    fig.savefig(str(outprefix) + '.manhattan.png'); plt.close(fig)

    # QQ
    p = np.sort(df['p'].values); n = len(p)
    obs = -np.log10(p); exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    if n > 250_000:
        idx = np.unique(np.r_[np.linspace(0, n - 1, 220_000, dtype=int), np.arange(max(0, n - 30_000), n)])
        exp2, obs2 = exp[idx], obs[idx]
    else:
        exp2, obs2 = exp, obs
    fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=140)
    ax.scatter(exp2, obs2, s=2, c='#2f5597', alpha=.55, linewidths=0)
    m = max(float(exp.max()), float(obs.max()))
    ax.plot([0, m], [0, m], color='red', lw=1)
    ax.set_xlabel('Expected -log10(P)'); ax.set_ylabel('Observed -log10(P)')
    ax.set_title(f'QQ: {title}\nlambdaGC={lam:.3f}')
    ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(str(outprefix) + '.qq.png'); plt.close(fig)
    return len(df), minp, lam


def process_one(args):
    trait, model = args
    path = RES_DIR / f'{trait}.{model}.csv'
    if not path.exists():
        return None
    try:
        df = read_rmvp(path)
    except Exception as e:
        print(f'[ERR] read {path}: {e}', flush=True); return None
    zh = TRAIT_ZH.get(trait, trait)
    # English-only title — avoids CJK font dependency in matplotlib.
    # Chinese name is kept in the output TSV (`trait_zh` column) for cross-reference.
    title = f'{trait} - {model}'
    outprefix = PLOT / f'{trait}.{model}'
    n_p, minp, lam = plot_pair(df, title, outprefix)
    tops = df.nsmallest(10, 'p').copy()
    tops.insert(0, 'trait', trait); tops.insert(1, 'trait_zh', zh); tops.insert(2, 'model', model)
    return {
        'trait': trait, 'trait_zh': zh, 'model': model,
        'n_p': n_p, 'min_p': minp, 'lambda_gc': lam,
        'n_bonferroni': int((df['p'] <= BONF).sum()),
        'n_suggestive':  int((df['p'] <= SUGG).sum()),
        'manhattan': f'{trait}.{model}.manhattan.png',
        'qq':        f'{trait}.{model}.qq.png',
        'tops': tops,
    }


def main():
    tasks = [(t, m) for t in TRAITS for m in MODELS]
    print(f'Total tasks: {len(tasks)} (traits={len(TRAITS)} x models={len(MODELS)})')
    results = []
    with Pool(processes=16) as pool:
        for i, r in enumerate(pool.imap_unordered(process_one, tasks), 1):
            if r is not None:
                results.append(r)
                print(f'  [{i}/{len(tasks)}] {r["trait"]}.{r["model"]}  '
                      f'minP={r["min_p"]:.2e}  lam={r["lambda_gc"]:.3f}  bonf={r["n_bonferroni"]}', flush=True)
            else:
                print(f'  [{i}/{len(tasks)}] MISSING', flush=True)

    summ_rows = [{k: v for k, v in r.items() if k != 'tops'} for r in results]
    summ = pd.DataFrame(summ_rows).sort_values(['trait', 'model'])
    summ.to_csv(SUMM / 'per_trait_model_summary.tsv', sep='\t', index=False)
    all_tops = pd.concat([r['tops'] for r in results], ignore_index=True)
    all_tops.to_csv(SUMM / 'top10_per_trait_model.tsv', sep='\t', index=False)
    sig  = all_tops[all_tops['p'] <= BONF].sort_values('p')
    sugg = all_tops[all_tops['p'] <= SUGG].sort_values('p')
    sig.to_csv (SUMM / 'bonferroni_significant_hits.tsv', sep='\t', index=False)
    sugg.to_csv(SUMM / 'suggestive_hits.tsv',             sep='\t', index=False)

    overlap = []
    for t in TRAITS:
        sub = sig[sig['trait'] == t]
        row = {'trait': t, 'trait_zh': TRAIT_ZH.get(t, t)}
        for m in MODELS:
            row[f'{m}_n'] = int((sub['model'] == m).sum())
        overlap.append(row)
    pd.DataFrame(overlap).to_csv(SUMM / 'model_overlap_by_trait.tsv', sep='\t', index=False)

    print(f'\nDONE: {len(results)} pairs, M_total={M_TOTAL:,}, BONF={BONF:.3e}')


if __name__ == '__main__':
    main()
