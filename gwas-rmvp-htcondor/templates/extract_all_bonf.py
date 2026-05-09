#!/usr/bin/env python3
"""
Scan every rMVP results CSV and export ALL SNPs passing the Bonferroni threshold
(not just the top-10 subset that plot_rmvp_all.py keeps).

Writes: post_gwas/summary/bonferroni_all_hits.tsv
"""
from pathlib import Path
from multiprocessing import Pool
import pandas as pd
import json

BASE    = Path('/abs/path/to/project/05.gwas')   # EDIT
RES_DIR = BASE / 'rMVP' / 'results'
SUMM    = BASE / 'post_gwas' / 'summary'
SUMM.mkdir(parents=True, exist_ok=True)

with (BASE / 'rMVP' / 'mvp.geno.map').open() as f:
    M_TOTAL = sum(1 for _ in f) - 1
BONF = 0.05 / M_TOTAL

cmap_path = BASE / 'pheno_category_map.json'
TRAIT_ZH = {}
if cmap_path.exists():
    m = json.loads(cmap_path.read_text())
    TRAIT_ZH = {v: k for k, v in m.get('original_to_safe', m).items()}


def process(path: Path):
    stem = path.stem                  # trait_XX.MODEL
    trait, model = stem.rsplit('.', 1)
    try:
        df = pd.read_csv(path, usecols=lambda c: True)
    except Exception as e:
        print(f'[ERR] {path.name}: {e}'); return None
    pcol = df.columns[-1]
    df = df.rename(columns={pcol: 'p'})
    df = df[(df['p'] > 0) & (df['p'] <= BONF)].copy()
    if df.empty:
        return None
    df.insert(0, 'trait', trait)
    df.insert(1, 'trait_zh', TRAIT_ZH.get(trait, trait))
    df.insert(2, 'model', model)
    return df


def main():
    paths = sorted(RES_DIR.glob('trait_*.csv'))
    paths = [p for p in paths if '_signals' not in p.stem]
    print(f'{len(paths)} CSV files; BONF={BONF:.3e}')
    with Pool(processes=8) as pool:
        frames = [f for f in pool.imap_unordered(process, paths) if f is not None]
    if not frames:
        print('No Bonferroni hits found.'); return
    out = pd.concat(frames, ignore_index=True).sort_values(['trait', 'model', 'p'])
    out_path = SUMM / 'bonferroni_all_hits.tsv'
    out.to_csv(out_path, sep='\t', index=False)
    print(f'{len(out)} hits -> {out_path}')


if __name__ == '__main__':
    main()
