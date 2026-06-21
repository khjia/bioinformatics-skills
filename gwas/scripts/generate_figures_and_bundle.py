#!/usr/bin/env python3
"""Generate post-GWAS figures and package report_bundle.tar.gz

Generates: QTL density heatmap, hotspot scatter, stacked Manhattan, bundles.
Usage: adjust BASE, POST, FIG paths below, then run.
"""
import csv, json, math, os, re, tarfile, shutil
from pathlib import Path
from collections import defaultdict, Counter
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ===== CONFIGURE THESE PATHS =====
BASE = Path("04.gwas")
POST = BASE / "post_gwas"
FIG = BASE / "figures"
FAST3_RESULTS = BASE / "fast3vmrmlm" / "results"
MAP_FILE = BASE / "rmvp" / "mvp.geno.map"
N_SNP = 6446649
DPI = 600
# =================================

BONF_LOG = -math.log10(0.05 / N_SNP)
FIG.mkdir(parents=True, exist_ok=True)

def save(fig, name):
    for fmt in ['png','pdf']:
        fig.savefig(FIG/f"{name}.{fmt}", dpi=DPI, bbox_inches='tight')
    plt.close(fig)

# Load chromosome order
chr_sizes, chr_order = {}, []
with open(MAP_FILE) as f:
    f.readline()
    for line in f:
        c, pos = line.strip().split('\t')[0], int(line.strip().split('\t')[2])
        if c not in chr_sizes: chr_order.append(c); chr_sizes[c] = 0
        chr_sizes[c] = max(chr_sizes[c], pos)
cum = {}; off = 0
for c in chr_order: cum[c] = off; off += chr_sizes[c]
colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2']

# 1. QTL density heatmap
loci = [row for row in csv.DictReader(open(POST/"consensus_loci.csv"))]
trait_chr = Counter((l['trait'][:60], l['CHROM']) for l in loci)
chrs = sorted(set(l['CHROM'] for l in loci), key=lambda x: int(x) if x.isdigit() else 999)
top_traits = [t for t,_ in Counter(l['trait'][:60] for l in loci).most_common(20)]
matrix = np.zeros((len(top_traits), len(chrs)))
for i,t in enumerate(top_traits):
    for j,c in enumerate(chrs): matrix[i,j] = trait_chr.get((t,c),0)
fig,ax = plt.subplots(figsize=(14,8))
im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(chrs))); ax.set_xticklabels(chrs, fontsize=9)
ax.set_yticks(range(len(top_traits))); ax.set_yticklabels(top_traits, fontsize=7)
ax.set_title(f"QTL Density (N={len(loci)} consensus SNPs)", fontsize=12)
plt.colorbar(im, ax=ax, label="SNPs"); plt.tight_layout()
save(fig, "qtl_density_heatmap")

# 2. Hotspot scatter
hotspots = [{'chr':r['CHROM'],'pos':int(r['WINDOW_START'])+250000,'n':int(r['N_TRAITS'])}
            for r in csv.DictReader(open(POST/"hotspots.tsv"), delimiter='\t')]
fig,ax = plt.subplots(figsize=(18,5))
hcolors = ['#e74c3c','#e67e22','#f1c40f','#2ecc71','#3498db']
for h in hotspots:
    if h['chr'] in cum:
        ax.scatter(cum[h['chr']]+h['pos'], h['n'], s=h['n']*2, alpha=0.6,
                   c=hcolors[min(h['n']-3,4)], edgecolors='none', rasterized=True)
ticks,labels = [],[]
for c in chr_order:
    if c in cum and c in chr_sizes: ticks.append(cum[c]+chr_sizes[c]/2); labels.append(c)
ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("N traits"); ax.set_title(f"Hotspots (N={len(hotspots)})", fontsize=12)
save(fig, "hotspot_scatter")

# 3. Stacked Manhattan (top 3, sampled every 10th SNP)
trait_bonf = defaultdict(int)
for r in csv.DictReader(open(POST/"summary.csv")): trait_bonf[r['trait']] += int(r['n_sig_bonf'])
top3 = [t for t,_ in sorted(trait_bonf.items(), key=lambda x:-x[1])[:3]]
fig,axes = plt.subplots(len(top3),1,figsize=(18,3*len(top3)),sharex=True)
if len(top3)==1: axes=[axes]
for ax,trait in zip(axes,top3):
    mid_f = FAST3_RESULTS/f"{trait}_midresult.csv"
    if not mid_f.exists(): continue
    pv_by_chr = defaultdict(list)
    for i,row in enumerate(csv.DictReader(open(mid_f))):
        if i%10!=0: continue
        try:
            pv=float(row['pval'])
            if 0<pv<=1: pv_by_chr[row['CHR']].append((int(row['POS']),pv))
        except: pass
    for i,c in enumerate(chr_order):
        if c not in pv_by_chr or c not in cum: continue
        pos=np.array([p[0] for p in pv_by_chr[c]])+cum[c]
        pv=np.clip(np.array([p[1] for p in pv_by_chr[c]]),1e-300,1)
        ax.scatter(pos,-np.log10(pv),s=0.3,c=colors[i%7],alpha=0.5,rasterized=True)
    ax.axhline(BONF_LOG,color='red',ls='--',lw=0.5)
    ax.set_ylabel("-log10(p)",fontsize=8); ax.set_title(trait[:80],fontsize=9)
axes[-1].set_xlabel("Chromosome",fontsize=10); plt.tight_layout()
save(fig, "stacked_manhattan_top3")

# 4. Package report_bundle.tar.gz
BUNDLE=POST/"report_bundle.tar.gz"; TMP=Path("/tmp/report_bundle")
for d in [TMP/"附图",TMP/"附表"]: d.mkdir(parents=True,exist_ok=True)
for src,dst,sd in [
    ("stacked_manhattan_top3","图1_堆叠曼哈顿图","附图"),
    ("qtl_density_heatmap","图2_QTL密度热图","附图"),
    ("hotspot_scatter","图3_共定位热点分布","附图"),
]:
    for fmt in ['png','pdf']:
        s=FIG/f"{src}.{fmt}"
        if s.exists(): shutil.copy(s,TMP/sd/f"{dst}.{fmt}")
for hp in sorted(FIG.glob("haplotype_*.png")):
    for fmt in ['png','pdf']:
        s=FIG/f"{hp.stem}.{fmt}"
        if s.exists(): shutil.copy(s,TMP/"附图"/f"图{4+int(hp.stem.split('_')[0])}_{hp.stem}.{fmt}")
for src,dst in [("summary.csv","表1_rMVP模型统计"),("consensus_loci.csv","表2_共识位点"),("hotspots.tsv","表3_共定位热点")]:
    s=POST/src
    if s.exists(): shutil.copy(s,TMP/"附表"/f"{dst}.tsv")
if (POST/"report.html").exists(): shutil.copy(POST/"report.html",TMP/"report.html")
with tarfile.open(BUNDLE,"w:gz") as tar: tar.add(TMP,arcname="report_bundle")
shutil.rmtree(TMP)
print(f"Bundle: {BUNDLE} ({BUNDLE.stat().st_size:,} bytes)")
