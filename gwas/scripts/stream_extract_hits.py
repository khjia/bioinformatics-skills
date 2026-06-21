#!/usr/bin/env python3
"""Streaming post-GWAS hit extraction from rMVP CSVs.
Memory-efficient: uses csv.DictReader, only stores significant SNPs.
For datasets with 100+ CSVs × ~700MB each (6M+ SNPs).
"""
import csv, json, sys, os, math
from pathlib import Path
from collections import defaultdict

# === CONFIGURE ===
RESULTS_DIR = "rmvp/results"        # path to rMVP CSV directory
OUT_DIR = "post_gwas"               # output directory
BONF = 0.05 / 6446649              # Bonferroni threshold
SUGGESTIVE = 1.0 / 6446649         # Suggestive threshold
# =================

results = Path(RESULTS_DIR)
out = Path(OUT_DIR)
out.mkdir(parents=True, exist_ok=True)

csv_files = [f for f in results.glob("*.csv") if "_signals" not in f.name]
print(f"Found {len(csv_files)} CSV files")

all_hits = []
snp_models = defaultdict(set)
trait_stats = []

for f in csv_files:
    stem = f.stem
    parts = stem.rsplit(".", 1)
    if len(parts) != 2:
        continue
    trait, model = parts
    if model not in ("GLM", "MLM", "FarmCPU"):
        continue
    
    n_sig = n_sug = 0
    min_p = 1.0
    pvals_for_lambda = []
    
    try:
        with open(f) as fh:
            reader = csv.DictReader(fh)
            pcol = reader.fieldnames[-1]
            for row in reader:
                try:
                    pval = float(row[pcol])
                except (ValueError, KeyError):
                    continue
                if pval <= 0 or pval > 1:
                    continue
                if pval < min_p:
                    min_p = pval
                pvals_for_lambda.append(pval)
                
                chr_val = row.get('CHROM', '')
                pos_val = int(row.get('POS', 0))
                snp_val = row.get('SNP', '')
                
                if pval < BONF:
                    n_sig += 1
                    snp_models[(trait, chr_val, pos_val, snp_val)].add(model)
                    all_hits.append({
                        'trait': trait, 'model': model,
                        'CHROM': chr_val, 'POS': pos_val, 'SNP': snp_val,
                        'pval': pval, 'significance': 'Bonferroni'
                    })
                elif pval < SUGGESTIVE:
                    n_sug += 1
                    all_hits.append({
                        'trait': trait, 'model': model,
                        'CHROM': chr_val, 'POS': pos_val, 'SNP': snp_val,
                        'pval': pval, 'significance': 'Suggestive'
                    })
    except Exception as e:
        print(f"  Error: {f.name}: {e}")
        continue
    
    if pvals_for_lambda:
        pvals_for_lambda.sort()
        n = len(pvals_for_lambda)
        lambda_gc = (-math.log10(pvals_for_lambda[n//2])) / (-math.log10(0.5))
    else:
        lambda_gc = float('nan')
    
    trait_stats.append({
        'trait': trait, 'model': model,
        'lambda_GC': round(lambda_gc, 4),
        'n_sig_bonf': n_sig, 'n_sig_suggestive': n_sug, 'min_p': min_p
    })
    
    if len(trait_stats) % 50 == 0:
        print(f"  ... {len(trait_stats)} files")

# Write outputs
with open(out / "all_hits.csv", "w", newline='') as f:
    w = csv.DictWriter(f, fieldnames=['trait','model','CHROM','POS','SNP','pval','significance'])
    w.writeheader()
    w.writerows(all_hits)
bonf_hits = sum(1 for h in all_hits if h['significance'] == 'Bonferroni')
print(f"Total hits: {len(all_hits)} ({bonf_hits} Bonferroni)")

with open(out / "summary.csv", "w", newline='') as f:
    w = csv.DictWriter(f, fieldnames=['trait','model','lambda_GC','n_sig_bonf','n_sig_suggestive','min_p'])
    w.writeheader()
    w.writerows(trait_stats)

consensus = {k: v for k, v in snp_models.items() if len(v) >= 2}
print(f"Consensus SNPs: {len(consensus)}")

WINDOW = 100000
loci = []
for (trait, chr_val, pos, snp), models in sorted(consensus.items()):
    loci.append({'trait': trait, 'CHROM': chr_val, 'POS': pos, 'SNP': snp,
                 'models': ','.join(sorted(models)), 'n_models': len(models)})
if loci:
    loci.sort(key=lambda x: (str(x['CHROM']), x['POS']))
    cur_chr, cur_loc, last_pos = None, 0, -WINDOW-1
    for l in loci:
        if l['CHROM'] != cur_chr or l['POS'] - last_pos > WINDOW:
            cur_loc += 1; cur_chr = l['CHROM']
        l['locus_id'] = f"{l['CHROM']}_{cur_loc}"
        last_pos = l['POS']
    with open(out / "consensus_loci.csv", "w", newline='') as f:
        w = csv.DictWriter(f, fieldnames=['trait','CHROM','POS','SNP','models','n_models','locus_id'])
        w.writeheader(); w.writerows(loci)
    print(f"Consensus loci: {cur_loc}")

WINDOW2 = 500000
trait_sig = defaultdict(set)
for h in all_hits:
    if h['significance'] == 'Bonferroni':
        w = (h['POS'] // WINDOW2) * WINDOW2
        trait_sig[(h['CHROM'], w)].add(h['trait'])
hotspots = [(c, w, len(t)) for (c, w), t in trait_sig.items() if len(t) >= 3]
with open(out / "hotspots.tsv", "w") as f:
    f.write("CHROM\tWINDOW_START\tN_TRAITS\n")
    for c, w, n in sorted(hotspots, key=lambda x: -x[2]):
        f.write(f"{c}\t{w}\t{n}\n")
print(f"Hotspots: {len(hotspots)}")
print("Done.")
