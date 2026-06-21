---
name: gwas
description: "Complete GWAS pipeline: Fast3VmrMLM v2.0 (multi-locus), post-GWAS colocalization/haplotype/enrichment, and publication-quality report generation. One-stop for plant GWAS on CNS cluster (ssh+nohup)."
---

# GWAS — Complete Pipeline (rMVP + Fast3VmrMLM + Post-GWAS)

Covers the full GWAS lifecycle: input prep → rMVP multi-model → Fast3VmrMLM multi-locus → colocalization hotspots → haplotype boxplots → GO/KEGG enrichment → figure generation → HTML report.

## When to use

Triggers:
- "跑 GWAS" / "Run GWAS" / "全基因组关联分析"
- User has PLINK BED + phenotype TSV + PCA covariates
- Multi-trait GWAS on plant populations (inbreds or selfing species)

## Tool selection quick guide

| Tool | Use case | N markers | Speed |
|---|---|---|---|
| **rMVP** (primary) | Multi-model (GLM/MLM/FarmCPU), single-locus p-values | up to 20M | ~2-5 min/trait on 16 cores |
| **Fast3VmrMLM** (complement) | Multi-locus QTL detection, controls false positives | up to 10M | ~3.5 min/trait on 40 cores |
| mrMLM R package | **DO NOT USE** — fails at >50k markers, internal bugs | <10k only | hours/trait |

**Recommendation**: Run both rMVP and Fast3VmrMLM. rMVP gives per-SNP significance (Manhattan/QQ), Fast3VmrMLM gives clean QTL lists. Cross-reference the results.

**When rMVP underperforms**: For structured plant populations with mixed binary/quantitative traits, rMVP's GLM may have inflated λ (>1.4), MLM may be overly conservative (0 Bonf hits), and FarmCPU may find only a handful of SNPs. In these cases, **Fast3VmrMLM alone can dramatically outperform** — a real example: 12M SNPs × 146 samples × 4 traits, rMVP FarmCPU found only 3-13 sig SNPs/trait, while Fast3VmrMLM found 164-16,649 sig SNPs/trait with reasonable λ (0.83-1.93).

## Execution preference

User has **disabled HTCondor/SLURM** (2026-06-10). All long tasks run via:

```bash
ssh CNS2 "cd /path && nohup /media/nfs1/hermes/miniforge3/bin/Rscript script.R > log 2>&1 &"
```

## ⚠️ PRE-FLIGHT: Verify the correct VCF/BED before starting

**This is the #1 mistake.** Projects often have multiple VCF files at different filtering stages. Before any GWAS, always:

1. **List all VCF files** in the project:
```bash
ls -lh <project>/vcf/   # or */vcf/
```

2. **Check record counts** for each VCF (use bcftools index --nrecords if CSI exists):
```bash
/home/khjia/.local/share/mamba/envs/sv/bin/bcftools index --nrecords <file>.vcf.gz
```

3. **Identify the final filtered VCF** — look for the one with ALL filters applied:
   - `pv25.raw.vcf.gz` → raw (DO NOT USE for GWAS)
   - `pv25.norm.vcf.gz` → normalized (DO NOT USE — still has indels, multi-allelic)
   - `pv25.snp.vcf.gz` → SNP-only (DO NOT USE — no MAF/missing filter)
   - `pv25.gt.vcf.gz` → GT-masked (DO NOT USE — no MAF/missing filter)
   - `pv25.final.vcf.gz` → **✅ THIS ONE** — GT mask + MAF + missingness

4. **Verify filtering matches user's standards** (from memory): `--minGQ 20 --minDP 4 --maxDP 1000 --max-missing 0.8 --maf 0.05 --max-alleles 2`

5. **Confirm SNP count with user** before converting to BED. Example expected counts:
   - Pea: 6.4M SNPs after filtering
   - Common bean (Pv25): 4.1M SNPs after filtering
   - If SNP count seems wrong (e.g. 12M vs expected 4M), STOP and verify which VCF was used.

6. **If converting from VCF, use plink2** with `--allow-extra-chr --threads 40`

Key rules:
- `Rscript` MUST use full path (`which Rscript` shows system default, not miniforge3)
- `$!` in double quotes → escape as `\\$!`
- `nohup` ensures survival after SSH disconnect
- `&` goes inside SSH quotes

### SSH batch submission (multiple jobs)

When submitting multiple `nohup` jobs via separate `ssh CNS2` calls, later SSH connections may timeout. Use a **single SSH connection** with a loop inside:

```bash
ssh CNS2 "
cd /abs/path/rmvp
for i in 1 2 3 4 5 6 7 8 9 10; do
    nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_rmvp.R \$i 10 \
        > ../logs/batch\${i}.log 2>&1 &
    echo \"Batch \$i PID=\$!\"
    sleep 2
done
"
```

When submitting multiple `nohup` jobs via separate `ssh CNS2` calls, CNS2 may reject or timeout later connections. Use a **single SSH connection** with a loop inside:

```bash
ssh CNS2 "
cd /abs/path/rmvp
for i in 1 2 3 4 5 6 7 8 9 10; do
    nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_rmvp.R \$i 10 \
        > ../logs/batch\${i}.log 2>&1 &
    echo \"Batch \$i PID=\$!\"
    sleep 2
done
"
```

Monitor with periodic checks:
```bash
# Check if still running
ssh CNS2 "ps -p PID -o etime,pid,cmd --no-headers"
# Tail log
ssh CNS2 "tail -30 /path/to/log"
# Count output files
ssh CNS2 "ls /path/to/results/ | wc -l"
```

## Project layout

```
<project>/04.gwas/          (or 05.gwas/ or 06.gwas/)
├── inputs/                 # source data
│   ├── pheno.tsv           # FID IID trait_01 trait_02 ...
│   ├── pca.tsv             # FID IID PC1 PC2 PC3
│   └── pheno_map.json      # {original: safe_name, ...}
├── scripts/
│   ├── prepare_gwas_inputs.R
│   ├── run_rmvp.R
│   ├── run_fast3vmrmlm.R
│   └── post_gwas/          # Python post-processing
│       ├── extract_hits.py
│       ├── plot_manhattan_qq.py
│       ├── colocalization.py
│       ├── haplotype_boxplots.py
│       └── build_html_report.py
├── logs/
├── rmvp/
│   ├── mvp.geno.{bin,desc,ind,map}
│   ├── kinship.rds
│   ├── kinship.fp          # genotype fingerprint for cache invalidation
│   └── results/            # trait_XX.{GLM,MLM,FarmCPU}.csv
├── fast3vmrmlm/
│   ├── results/            # trait_XX_midresult.csv, trait_XX_result.xlsx
│   └── preKinship.csv
├── post_gwas/
│   ├── hotspots.tsv
│   ├── high_confidence_loci.tsv
│   ├── plots/              # Manhattan, QQ, hotspots, haplotype boxplots
│   ├── enrichment/         # GO/KEGG results
│   └── report.html         # self-contained HTML
└── figures/                # final delivery figures (600 DPI PNG + PDF)
```

---

# Part 1: rMVP — Multi-Model GWAS

## Environment

```bash
# rMVP is in miniforge3 base R
/media/nfs1/hermes/miniforge3/bin/Rscript
```

Version: rMVP 1.4.6 (supports GLM/MLM/FarmCPU only; NO BLINK).

## Step 1: Convert BED → rMVP big.matrix

```r
library(rMVP)
MVP.Data(fileBed="genotype_prefix", fileOut="mvp", verbose=TRUE)
```

Produces `mvp.geno.{bin,desc,ind,map}`. 5-20 min for ~16M SNPs × 140 samples.

## Step 2: Run GWAS

### CRITICAL: phe must be data.frame with ALL genotype individuals

rMVP requires the phenotype data.frame to contain **every individual** in the genotype (with NA for missing values). Passing only valid samples causes:
> "The number of individuals in phenotype and genotype doesn't match!"

```r
# WRONG — only includes valid samples → "individuals mismatch" error
valid <- !is.na(pheno[[trait_name]])
phe_df <- data.frame(taxon = pheno$IID[valid], ...)

# WRONG — phe is a vector → "argument is of length zero"
y <- pheno[[trait_name]]
names(y) <- pheno$IID
MVP(phe = y, ...)

# CORRECT — ALL genotype individuals, NA for missing phenotypes
ind <- readLines("mvp.geno.ind")
idx <- match(ind, pheno$IID)          # reorder pheno to match genotype
pheno_all <- pheno[idx, ]
y <- as.numeric(pheno_all[[trait_name]]) # includes NAs
phe_df <- data.frame(taxon = ind, trait = y, stringsAsFactors = FALSE)
colnames(phe_df)[2] <- trait_name
MVP(phe = phe_df, ...)
```

### Full driver script

```r
library(rMVP)

# Load genotype
genotype <- attach.big.matrix("mvp.geno.desc")
map <- read.table("mvp.geno.map", header=TRUE)
ind <- readLines("mvp.geno.ind")     # individual IDs in genotype order
cat(sprintf("Genotype: %d SNPs x %d samples\n", nrow(map), length(ind)))

# Load phenotype — MUST use sep="\t" explicitly (NFS-safe)
pheno <- read.table("inputs/pheno.tsv", header=TRUE, sep="\t",
                    na.strings=c("", "NA", "."), stringsAsFactors=FALSE)
# Reorder pheno to match genotype ind order
idx <- match(ind, pheno$IID)
if (any(is.na(idx))) stop("Some genotype individuals missing from pheno!")
pheno <- pheno[idx, ]

# Load PCA — MUST use sep="\t" explicitly
pc <- read.table("inputs/pca.tsv", header=TRUE, sep="\t")
pc_idx <- match(ind, pc$IID)
pc <- pc[pc_idx, ]
CV <- as.matrix(pc[, c("PC1","PC2","PC3")])

# Kinship with fingerprint cache
fp_file <- "kinship.fp"
fp_current <- paste(file.info("mvp.geno.bin")$size,
                    file.info("mvp.geno.bin")$mtime,
                    tools::md5sum("mvp.geno.bin"), sep="_")

if (file.exists("kinship.rds") && file.exists(fp_file) &&
    readLines(fp_file) == fp_current) {
    K <- readRDS("kinship.rds")
    cat("Loaded cached K\n")
} else {
    K <- MVP.K.VanRaden(genotype, cpu=16, verbose=TRUE)  # NOTE: cpu not ncpus!
    saveRDS(K, "kinship.rds")
    writeLines(fp_current, fp_file)
    cat("Computed new K\n")
}

results_dir <- "results"
dir.create(results_dir, showWarnings=FALSE)

# Run per trait
traits <- colnames(pheno)[-(1:2)]
for (trait_name in traits) {
    cat(sprintf("\n=== %s ===\n", trait_name))

    y <- as.numeric(pheno[[trait_name]])
    n_valid <- sum(!is.na(y))
    if (n_valid < 50) { cat(sprintf("SKIP: %d valid\n", n_valid)); next }

    # ALL genotype individuals — NAs for missing phenotypes
    phe_df <- data.frame(taxon = ind, trait = y, stringsAsFactors = FALSE)
    colnames(phe_df)[2] <- trait_name

    for (attempt in 1:3) {
        res <- tryCatch({
            MVP(
                phe = phe_df, geno = genotype, map = map, K = K,
                CV.GLM = CV, CV.MLM = CV, CV.FarmCPU = CV,
                nPC.GLM = 3, nPC.MLM = 3, nPC.FarmCPU = 3,
                ncpus = 16,
                method = c("GLM", "MLM", "FarmCPU"),
                maxLoop = 10, method.bin = "static",
                threshold = 0.05 / nrow(map),
                file.output = TRUE, file.type = "csv",
                outpath = results_dir      # write directly to results/
            )
            TRUE
        }, error = function(e) {
            cat(sprintf("[RETRY %d] %s: %s\n", attempt, trait_name, e$message))
            gc(); Sys.sleep(30); return(FALSE)
        })
        if (isTRUE(res)) break
        if (attempt == 3) cat(sprintf("[FAIL] %s\n", trait_name))
    }
}
```

## rMVP Critical Bugs (v1.4.6)

| Bug | Symptom | Fix |
|---|---|---|
## rMVP Critical Bugs (v1.4.6)

| Bug | Symptom | Fix |
|---|---|---|
| phe as vector | "argument is of length zero" | Wrap in data.frame with taxon column |
| phe as vector | "argument is of length zero" | Wrap in data.frame with taxon column |
| phe only valid samples | "number of individuals... doesn't match!" | Include ALL genotype individuals, NA for missing phenotype |
| MVP.K.VanRaden ncpus | "unused argument (ncpus = ncpus)" | Use `cpu=N`, not `ncpus` |
| BLINK not available | "Unknow method: BLINK" | Only `c("GLM","MLM","FarmCPU")` |
| file.output=FALSE | Empty 3-byte CSV output | Always `file.output=TRUE` |
| file.cols parameter | "unused argument (file.cols = ...)" | Remove it — not in v1.4.6 |
| plot parameter | "unused argument (plot = FALSE)" | Cannot disable built-in plotting in v1.4.6 |
| CSV column naming | p-value header is `trait.MODEL` | Use `df.columns[-1]` |
| CSV missing CHROM/POS | Output has only Effect, SE, p_value | Merge with map by row index in post-processing |
| CSV missing CHROM/POS (from fwrite) | rMVP results written with `fwrite()` don't include map columns | Python: `pd.concat([map[['SNP','CHROM','POS']], result], axis=1)` |
| NFS file.rename | Silently fails | Use `outpath=results_dir` to write directly |
| NFS file.copy+remove | Files can be deleted before copy completes | Use `outpath` instead, never copy+remove on NFS |
| NFS read.table sep auto-detect | "line 1 did not have N elements" | Always use `sep="\\t"` explicitly |
| FarmCPU on binary | λ=0, all p≈1 | Recode to 0/1, drop rare category |
| FarmCPU on 3+ cats | λ=0, no hits | Collapse to binary or use MLM only |
| Graphics device leak | R hangs/stalls after ~10 traits | `while (length(dev.list()) > 0) try(dev.off())` before & after MVP() |
| pdf("/dev/null") | Still computes plots, just redirects | Let rMVP plot normally, clean up devices after |
| FarmCPU on 3+ cats | λ=0, no hits | Collapse to binary or use MLM only |
| Graphics device leak | R hangs/stalls after ~10 traits | `while (length(dev.list()) > 0) try(dev.off())` before & after MVP() |
| pdf("/dev/null") | Still computes plots, just redirects | Use `while(dev.list()>0) dev.off()` instead — let rMVP plot but clean up |

## Parameters

| Param | Value | Rationale |
|---|---|---|
| Models | GLM, MLM, FarmCPU | rMVP 1.4.6 supported |
| PC covariates | 3 | Standard for structured populations |
| K method | VanRaden | Default, robust |
| Bonferroni | 0.05 / nSNP | Standard genome-wide |
| Suggestive | 1 / nSNP | Discovery threshold |
| ncpus | 16 | Match CNS2 per-job allocation |
| λ_GC pass | [0.85, 1.15] | Well-calibrated range |
| min valid samples | 50 | Too few = unreliable fitting |

## Post-processing: Manhattan + QQ

Generate from rMVP CSV output using Python/matplotlib:

```python
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

DPI = 600  # REQUIRED

def save_both(fig, path_no_ext):
    fig.savefig(f'{path_no_ext}.png', dpi=DPI, bbox_inches='tight')
    fig.savefig(f'{path_no_ext}.pdf', bbox_inches='tight')

# Manhattan
fig, ax = plt.subplots(figsize=(14, 4))
for chr_num in chr_order:
    chr_df = df[df['CHROM'] == chr_num].sort_values('POS')
    pv = -np.log10(chr_df.iloc[:, -1].values.clip(1e-300))
    ax.scatter(pos, pv, s=0.3, c=color, alpha=0.6, rasterized=True)

# QQ
pvals = df.iloc[:, -1].dropna().values
pvals = pvals[pvals > 0]
obs = -np.log10(np.sort(pvals))
exp = -np.log10((np.arange(1, len(pvals)+1) - 0.5) / len(pvals))
lambda_gc = np.median(obs) / (-np.log10(0.5))
ax.scatter(exp, obs, s=0.5, alpha=0.3, rasterized=True)
```

---

# Part 2: Fast3VmrMLM v2.0 — Multi-Locus GWAS

Standalone R package (different from mrMLM R package). GitHub: `YuanmingZhang65/Fast3VmrMLM`.

## Installation & TBB Fix

Pre-compiled binary needs old TBB 2020.2:

```bash
# Old TBB setup (one-time)
mkdir -p /media/nfs1/hermes/lib/tbb2020
cd /tmp && curl -sL --max-time 120 \
  -o tbb_old.tar.bz2 \
  "https://conda.anaconda.org/conda-forge/linux-64/tbb-2020.2-h4bd325d_1.tar.bz2"
tar xjf tbb_old.tar.bz2
cp lib/libtbb* /media/nfs1/hermes/lib/tbb2020/

# Install Fast3VmrMLM in R base env
LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2 \
  /media/nfs1/hermes/miniforge3/bin/Rscript -e \
  'install.packages("Fast3VmrMLM", repos=NULL)'
```

## Input Format & Preparation

### Phenotype CSV — MUST match FAM sample order

**Robust method** (handles any FAM ordering, any ID format):

```bash
{
    echo -n "<Phenotype>"
    head -1 pheno.tsv | cut -f3- | tr '\t' ',' | sed 's/^/,/'
    awk -F'\t' 'NR==FNR && FNR>1 {pheno[$2]=$0; next}
               $1 in pheno {
                   split(pheno[$1], p, "\t")
                   printf "%s", p[2]
                   for(i=3; i<=length(p); i++) printf ",%s", (p[i]==""?"NA":p[i])
                   printf "\n"
               }' \
        pheno.tsv \
        <(cut -f2 genotype.fam)
} > pheno_fast3.csv
```

**Why this over the one-liner sort+join approach**: The sort-based method fails when FAM and pheno have different ID sort orders (common with numeric IDs like `124,049,117` vs `001,002,003`). The hash-based lookup is order-independent.
```

Result format:
```
<Phenotype>,trait1,trait2,...
sample1_in_fam_order,val1,val2,...
```

- Header MUST be `<Phenotype>,trait1,...` (angle brackets)
- Sample IDs must match FAM **and** be in FAM order
- Missing values: use empty or `NA`

### PCA CSV — MUST match FAM sample order

```bash
# Reorder PCA to match FAM order
{
    echo "<pca>,pc_1,pc_2,pc_3"
    awk 'NR==FNR{ids[$2]=NR; order[NR]=$2; next}
         FNR>1{val[$2]=$3","$4","$5}
         END{for(i=1;i<=NR;i++){sid=order[i]; if(sid in val) print sid","val[sid]; else print sid",0,0,0"}}' \
        genotype.fam pca.tsv
} > pca_fast3.csv
```

Result format:
```
<pca>,pc_1,pc_2,pc_3
sample1_in_fam_order,pc1_val,pc2_val,pc3_val
```

- Header MUST be `<pca>,pc_1,pc_2,pc_3`
- Samples in FAM order, matching pheno file exactly

## Running

```r
library(Fast3VmrMLM)

Fast3VmrMLM(
    fileGen    = "genotype",       # PLINK prefix, no .bed
    filePhe    = "trait.csv",
    filePS     = "pca.csv",
    PopStrType = "PC",
    fileOut    = "results/",       # MUST end with / and dir must exist
    genoType   = "SNP",
    trait      = 1:17,             # column indices
    svrad      = 20,               # search radius (kb)
    svpal      = 0.01,             # critical P-value
    svmlod     = 3,                # critical LOD
    nThreads   = 40,
    DrawPlot   = TRUE,
    Plotformat = "*.tiff"
)
```

### Full launch

```bash
ssh CNS2 "cd /abs/path/04.gwas/fast3vmrmlm && \
  mkdir -p results && \
  export LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2 && \
  nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_fast3vmrmlm.R > ../logs/fast3vmrmlm.log 2>&1 & \
  echo PID=\$!"
```

## Output

- `{trait}_midresult.csv` — genome-wide scan (MarkerID, CHR, POS, Waldst, pval)
- `{trait}_result.xlsx` — significant QTL (LOD, add, dom, r², P-value, SIG/SUG)
- `{trait}_Manhattan_plot.tiff` — built-in (huge, use PIL MAX_IMAGE_PIXELS=None)

## Fast3VmrMLM Pitfalls

| Issue | Fix |
|---|---|
## Fast3VmrMLM Pitfalls

| Issue | Fix |
|---|---|
| TBB symbol error | `LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2` |
| Output dir doesn't exist | `mkdir -p` before running |
| XLSX cols are character | `as.numeric()` before use |
| `Position.(bp)` in pandas | Use `[c for c in df if 'Position' in c][0]` | **WRONG** — actual columns are `CHR,POS,MarkerID,Waldst,pval`. The position column is `POS` (not `Position` or `Position.(bp)`), chromosome is `CHR` (not `Chromosome`), SNP ID is `MarkerID` (not `RS#`). Post-GWAS scripts must use these exact names. |
| TIFF DecompressionBombError | `Image.MAX_IMAGE_PIXELS = None` |
| Full Rscript path needed | SSH doesn't inherit PATH — use `/media/nfs1/hermes/miniforge3/bin/Rscript` |
| First trait slow (18+ min) | Kinship computed on first trait; subsequent traits reuse it (~2-3 min each) |
| "Drawing manplot failed" | **Benign** — plot generation fails on 6M+ SNPs but doesn't crash; process continues to next trait |
| Process dies mid-batch (no error) | Restart from breakpoint using trait indices (see below) |
| Pheno/PCA not in FAM order | Fast3VmrMLM matches by row position, not sample name — MUST reorder to FAM order |

## Fast3VmrMLM restart from breakpoint

Processes can die silently mid-run (no OOM, no error in log). Since Fast3VmrMLM writes midresult CSV per trait and the `trait` parameter accepts column indices, just restart from where each batch left off:

```r
# Restart script — takes start and end trait column indices
args <- commandArgs(trailingOnly = TRUE)
trait_start <- as.integer(args[1])
trait_end <- as.integer(args[2])

Fast3VmrMLM(
    fileGen = "../inputs/genotype",
    filePhe = "pheno_fast3.csv",
    filePS  = "pca_fast3.csv",
    PopStrType = "PC",
    fileOut = "results/",
    genoType = "SNP",
    trait   = trait_start:trait_end,  # restart from breakpoint
    svrad   = 20, svpal = 0.01, svmlod = 3,
    nThreads = 40,
    DrawPlot = TRUE, Plotformat = "*.tiff"
)
```

```bash
# Example: batch originally did traits 1-19, completed 9, restart with 10-19
ssh CNS2 "
export LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2
nohup /media/nfs1/hermes/miniforge3/bin/Rscript restart.R 10 19 > log 2>&1 &
"
```

Monitor with: `strings log | grep -c 'completes'` to count finished traits.

## Fast3VmrMLM batch submission (114 traits example)

For 114 traits on 224 samples (6.4M SNPs), use **6 batches** (~19 traits each):
```bash
ssh CNS2 "
for i in 1 2 3 4 5 6; do
    mkdir -p /abs/path/fast3vmrmlm/results
    export LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2
    nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_fast3vmrmlm.R \$i 6 > ../logs/fast3vmrmlm_batch\${i}.log 2>&1 &
    sleep 3
done
"
```
6 batches × ~19 traits = 114 total. Each batch completes 4-5 traits in the first 30 min. Full run ~4-5 hours. Monitor: `ls results/*_midresult.csv | wc -l`

## Fast3VmrMLM Input Format & Preparation

### Phenotype CSV — MUST match FAM sample order

**Robust method** (handles any FAM ordering, any ID format):

```bash
{
    echo -n "<Phenotype>"
    head -1 pheno.tsv | cut -f3- | tr '\t' ',' | sed 's/^/,/'
    awk -F'\t' 'NR==FNR && FNR>1 {pheno[$2]=$0; next}
               $1 in pheno {
                   split(pheno[$1], p, "\t")
                   printf "%s", p[2]
                   for(i=3; i<=length(p); i++) printf ",%s", (p[i]==""?"NA":p[i])
                   printf "\n"
               }' \
        pheno.tsv \
        <(cut -f2 genotype.fam)
} > pheno_fast3.csv
```

**Why this over the one-liner sort+join approach**: The sort-based method fails when FAM and pheno have different ID sort orders (common with numeric IDs like `124,049,117` vs `001,002,003`). The hash-based lookup is order-independent.
```

### PCA CSV — MUST match FAM sample order

```bash
{
    echo "<pca>,pc_1,pc_2,pc_3"
    awk 'NR==FNR{ids[$2]=NR; order[NR]=$2; next}
         FNR>1{val[$2]=$3","$4","$5}
         END{for(i=1;i<=NR;i++){sid=order[i]; if(sid in val) print sid","val[sid]; else print sid",0,0,0"}}' \
        genotype.fam pca.tsv
} > pca_fast3.csv
```

Result format:
```
<pca>,pc_1,pc_2,pc_3
sample1_in_fam_order,pc1_val,pc2_val,pc3_val
```
| First trait slow (18+ min) | Kinship computed on first trait; subsequent traits reuse it (~2-3 min each) |
| "Drawing manplot failed" | **Benign** — plot generation fails on 6M+ SNPs but doesn't crash; process continues to next trait |
| Process dies mid-batch (no error) | Restart from breakpoint using trait indices (see below) |
| Pheno/PCA not in FAM order | Fast3VmrMLM matches by row position, not sample name — MUST reorder to FAM order |

## Fast3VmrMLM restart from breakpoint

Processes can die silently mid-run (no OOM, no error in log). Since Fast3VmrMLM writes midresult CSV per trait and the `trait` parameter accepts column indices, just restart from where each batch left off:

```r
# Restart script — takes start and end trait column indices
args <- commandArgs(trailingOnly = TRUE)
trait_start <- as.integer(args[1])
trait_end <- as.integer(args[2])

Fast3VmrMLM(
    fileGen = "../inputs/genotype",
    filePhe = "pheno_fast3.csv",
    filePS  = "pca_fast3.csv",
    PopStrType = "PC",
    fileOut = "results/",
    genoType = "SNP",
    trait   = trait_start:trait_end,  # restart from breakpoint
    svrad   = 20, svpal = 0.01, svmlod = 3,
    nThreads = 40,
    DrawPlot = TRUE, Plotformat = "*.tiff"
)
```

```bash
# Example: batch originally did traits 1-19, completed 9, restart with 10-19
ssh CNS2 "
export LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2
nohup /media/nfs1/hermes/miniforge3/bin/Rscript restart.R 10 19 > log 2>&1 &
"
```

Monitor with: `strings log | grep -c 'completes'` to count finished traits.

## Fast3VmrMLM batch submission (114 traits example)

For 114 traits on 224 samples (6.4M SNPs), use **6 batches** (~19 traits each):
```bash
ssh CNS2 "
for i in 1 2 3 4 5 6; do
    mkdir -p /abs/path/fast3vmrmlm/results
    export LD_PRELOAD=/media/nfs1/hermes/lib/tbb2020/libtbb.so.2
    nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_fast3vmrmlm.R \$i 6 > ../logs/fast3vmrmlm_batch\${i}.log 2>&1 &
    sleep 3
done
"
```
6 batches × ~19 traits = 114 total. Each batch completes 4-5 traits in the first 30 min. Full run ~4-5 hours. Monitor: `ls results/*_midresult.csv | wc -l`

## Parameter Tuning

| Parameter | Default | Guideline |
|---|---|---|
| svpal | 0.01 | <1000 samples: 0.01-0.05; >4000: 0.01; >20000: 1e-5 |
| svrad | 20 | Search radius kb around significant markers |
| svmlod | 3 | Higher = more stringent |

---

# Part 3: Post-GWAS Analysis

## 3a. Multi-Model Consensus (rMVP) — USE STREAMING, NOT PANDAS

**CRITICAL**: For >100 rMVP CSV files (each ~700MB for 6M SNPs), loading all into pandas
causes OOM or >10 min load time. Use Python's csv.DictReader to stream one file at a time.

```python
import csv
from collections import defaultdict

snp_models = defaultdict(set)  # (trait, chr, pos, snp) -> set of models
all_hits = []

for csv_file in results_dir.glob("*.csv"):
    trait, model = parse_filename(csv_file)
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        pcol = reader.fieldnames[-1]  # last column = p-value
        for row in reader:
            pval = float(row[pcol])
            if pval < bonf:
                key = (trait, row['CHROM'], int(row['POS']), row['SNP'])
                snp_models[key].add(model)
                all_hits.append({...})

# After streaming: consensus = SNPs significant in >=2 models
consensus = {k: v for k, v in snp_models.items() if len(v) >= 2}
```

Memory-efficient: only stores significant SNPs (~tens of thousands), not all 6.4M × 482 rows.

## 3b. Manhattan + QQ from Fast3VmrMLM Midresult (Standalone)

When running Fast3VmrMLM standalone (no rMVP), generate Manhattan/QQ directly from midresult CSVs.

**CRITICAL: User requires ALL SNPs on Manhattan plots — NO SAMPLING.** Sampling to 500K-800K from 6M-16M SNPs is not acceptable. Every SNP must be plotted. Specification: Three common mistakes:

| Mistake | Symptom | Why wrong |
|---------|---------|-----------|
| Plot only Bonf-sig SNPs | Empty figure (4 dots for stem_color) | Background density IS the Manhattan plot |
| Mix `ax.plot(',')` + `ax.scatter` for sig overlay | "Significant and non-significant look different" | Two different marker types = visual inconsistency |
| Skip Y-axis capping | Extremely stretched y-axis from few ultra-low p-values | Makes the figure unreadable |

**PREFERRED approach: `ax.scatter` with chunking, Y-axis capped at 15 (yuanbaofeng style).**
This is what users expect visually — uniform small circles, alternating chromosome colors, red threshold line.

```python
# Manhattan — ALL SNPs, uniform scatter, numpy + chunking, Y capped at 15
from matplotlib.ticker import MaxNLocator

fig, axes = plt.subplots(len(traits), 1, figsize=(18, 3.5*len(traits)), sharex=True)
fig.subplots_adjust(hspace=0.15)
chr_colors = ['#1f77b4', '#ff7f0e']  # alternating blue/orange

for i, trait in enumerate(traits):
    ax = axes[i]
    # Read ALL SNPs into numpy arrays
    chrs, poss, pvs = [], [], []
    with open(f"{trait}_midresult.csv") as fh:
        for row in csv.DictReader(fh):
            try:
                chrs.append(int(row['CHR']))  # Note: CHR not Chromosome
                poss.append(int(row['POS']))  # Note: POS not Position
                pvs.append(float(row['pval']))
            except: pass
    chrs = np.array(chrs); pvs = np.array(pvs)
    
    # Genomic positions, -log10(p) capped at 15
    x_all = np.array([cum_pos.get(c, 0) + p for c, p in zip(chrs, poss)])
    y_all = -np.log10(np.clip(pvs, 1e-15, 1))
    y_all = np.minimum(y_all, 15)
    
    # Scatter by chromosome — 500k point chunks (memory safety)
    for color_idx, color in enumerate(chr_colors):
        mask = (chrs % 2 == color_idx)
        idx = np.where(mask)[0]
        for start in range(0, len(idx), 500000):
            end = min(start + 500000, len(idx))
            chunk = idx[start:end]
            ax.scatter(x_all[chunk], y_all[chunk], s=0.3, c=color,
                       alpha=0.6, rasterized=True, edgecolors='none')
    
    ax.axhline(-np.log10(BONF), color='red', linestyle='--', linewidth=0.6, alpha=0.7)
    ax.set_ylabel(trait.replace('_', '\\_'), fontsize=9)
    ax.set_ylim(0, 16)
    ax.yaxis.set_major_locator(MaxNLocator(4))
```

Key principles:
- **ALL SNPs rendered IDENTICALLY** — no separate sig overlay. Significance = position above threshold line.
- **Y-axis capped at 15** (ylim 0-16) — prevents ultra-low p-values from stretching the plot.
- **`edgecolors='none'`** — critical for performance with millions of points.
- **500k point chunks** — prevents single scatter call from OOM on 4M+ points.
- **`rasterized=True`** — converts to bitmap, keeps vector axes sharp. File size ~20MB for 4M SNPs.
- Even/odd chromosomes get alternating colors.

**FALLBACK for >10M SNP datasets**: If scatter is too slow, use `ax.plot(..., ',')` pixel marker for speed:

**PREFERRED method: R/ggplot2 on CNS2.** Proven to handle 16M SNPs with `data.table::fread` + `geom_point(size=0.25)`. Run on CNS2 for local NFS speed — Python cross-node NFS stalls on 500MB files.

```r
# Manhattan — ALL SNPs, R/ggplot2, CNS2
library(ggplot2)
dt <- data.table::fread(f, select = c("CHR","POS","pval"))
dt <- dt[is.finite(dt$pval) & dt$pval > 0 & dt$pval <= 1, ]
dt$cum_pos <- dt$POS + cum_start[as.character(dt$CHR)]
dt$neglogp <- -log10(dt$pval)

ggplot(dt, aes(x = cum_pos, y = neglogp, color = CHR)) +
    geom_point(size = 0.25, alpha = 0.4, stroke = 0, shape = 16) +
    scale_color_manual(values = colors, guide = "none") +
    geom_hline(yintercept = bonf_log, color = "red", linewidth = 0.3, linetype = "dashed") +
    labs(title = trait, x = "Chromosome", y = expression(-log[10](p))) +
    theme_bw(base_size = 9) + theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylim(0, max(dt$neglogp, bonf_log + 1, na.rm = TRUE))

ggsave(sprintf("manhattan_full_%s.png", trait), p, width = 14, height = 3.5, dpi = 600)
```

Key specs: `size=0.25, alpha=0.4, stroke=0, shape=16`, PNG only (no PDF with 6M+ points), 600 DPI, 14×3.5 inches. ~6M SNP: ~2 min/trait. ~16M SNP: ~5 min/trait.

```python
# QQ: read ALL p-values from midresult (no sampling needed, fast for QQ)
pvals = []
with open(f"{trait}_midresult.csv") as fh:
    for row in csv.DictReader(fh):
        try: pv = float(row['pval'])
        except: continue
        if 0 < pv <= 1: pvals.append(pv)
n = len(pvals)
obs = -np.log10(np.sort(pvals))
exp = -np.log10((np.arange(1, n+1) - 0.5) / n)
ax.scatter(exp, obs, s=0.3, alpha=0.3, color='#333333', rasterized=True)
ax.plot([0, max(exp)], [0, max(exp)], 'r-', linewidth=0.8)
```

Scan genome in 500 kb windows. Count traits with ≥1 QTL per window.

```python
window_size = 500000
co_loc = {}
for qtl in all_qtl:
    w = (qtl['pos'] // window_size) * window_size
    co_loc.setdefault((qtl['chr'], w), set()).add(qtl['trait'])

hotspots = [(chr, w, len(traits)) for (chr, w), traits in co_loc.items() if len(traits) >= 3]
```

Outputs: `hotspots.tsv`, colocalization scatter plot, QTL density heatmap (trait × chromosome).

## 3c. Haplotype-Phenotype Boxplots

For top N QTL (by LOD), extract genotypes and plot phenotype by genotype class.

**CRITICAL: BIM SNP naming and plink2 extraction**

Many PLINK BIM files have `.` as SNP names (identified by chr:pos only). `--extract` by SNP name won't find anything. Always use positional extraction:

```bash
# CORRECT: position-based extraction (works even when BIM SNP name is ".")
plink2 --bfile genotype --chr $CHR --from-bp $POS --to-bp $POS --export A --out top_qtl_geno
```

```python
# IMPORTANT: consensus_loci trait names use "." separators (from rMVP CSV)
# but phenotype file uses "_" (safe names). Must normalize:
import re

# Consensus trait → safe name matching  
for orig_name in locus_traits:
    candidates = [
        orig_name,
        orig_name.replace(';', '_'),
        re.sub(r'\._|_\.|\.', '_', orig_name)  # "0._purple1" → "0_purple1"
    ]
    candidates = [re.sub(r'_+', '_', c) for c in candidates]  # collapse "__"→"_"
    for c in candidates:
        if c in pheno_trait_names:
            safe_name = c
            break
```

**Pitfall**: `Stem_color_greem0._purple1_2021` with naive `.replace('.','_')` → `Stem_color_greem0__purple1_2021` (double underscore). Must collapse with regex.

**Genotype column detection**: plink2 `--export A` produces columns like `._G` or `._T` (reference/alternate allele count). When using Python's `csv.DictReader` with `delimiter='\t'`, detect the genotype column by finding any column name not in `{'FID','IID','PAT','MAT','SEX','PHENOTYPE'}`:

```python
plink_meta = {'FID','IID','PAT','MAT','SEX','PHENOTYPE'}
geno_col = None
for row in csv.DictReader(f, delimiter='\t'):
    if not geno_col:
        for c in row:
            if c not in plink_meta:
                geno_col = c  # will be '._G' or '._T'
                break
    if geno_col:
        v = row.get(geno_col, '')
        if v and v.strip() and v != 'NA':
            geno_vals[row['IID']] = int(float(v))
```

## 3d. GO/KEGG Enrichment

Map QTL to nearest genes → eggNOG annotation → clusterProfiler enrichment.

**CRITICAL: Full genome background is required.** Without annotating the entire genome with eggNOG, enrichment is just descriptive frequency, not statistical enrichment.

**CRITICAL: Gene ID version suffix matching.** GFF gene IDs often lack version suffixes (e.g. `Pisat01G0001300`), but eggNOG-mapper output includes them (e.g. `Pisat01G0001300.1`). If candidate genes and term2gene use different ID formats, enrichment silently returns 0 results. Fix: add `.1` suffix to candidate genes before enrichment:

```python
# In extract_nearby_genes.py: add .1 suffix for emapper compatibility
if '.' not in gene:
    gene = gene + '.1'
```

Verify with: `comm -12 <(cut -f1 candidate_genes.tsv | tail -n+2 | sort) <(sort term2gene/gene_universe.tsv) | wc -l` — should be > 80% of candidate genes.

**GFF chromosome naming mismatch.** Consensus loci use numeric chromosomes (`1,2,...7`) but GFF may use prefixed names (`chr01, chr02,...`). Map before gene extraction:

```python
chr_map = {str(i): f"chr{i:02d}" for i in range(1, 8)}
gff_chr = chr_map.get(chr_val, chr_val)  # 1→chr01, etc.
```

```bash
# Extract all proteins from genome GFF
gffread -y all_proteins.faa -g genome.fa genome.gff3

# Run eggNOG on full genome (CNS2, ~40 min for 30k genes)
ssh CNS2 "export PATH=/media/nfs1/hermes/miniforge3/bin:\$PATH && \
  emapper.py --data_dir /media/nfs1/hermes/db/eggnog \
    -i all_proteins.faa -o emapper/full_genome \
    --cpu 40 --tax_scope 33090 --override --temp_dir /tmp"
```

Then build term2gene from full annotations → `clusterProfiler::enricher()` with QTL genes as foreground and all annotated genes as universe.

**Enrichment result reporting rules:**
- Report ALL terms with FDR < 0.05 in the results table — never delete based on subjective relevance
- Check for cross-annotation: if two terms share the same gene set, merge into the most specific one
- Flag marginal terms (p near 0.05, few genes) as "exploratory"
- Biological interpretation goes in discussion/footnotes, not by deleting results
**Phase 2 — Full-genome annotation** (2-3h for 30k genes, required for enrichment):
```bash
# Extract all proteins from GFF3 (use chrXXa-only if genome is haploid-a)
gffread -y full_genome_proteins.faa -g genome.a.fasta genome_a_only.gff3

emapper.py --data_dir /media/nfs1/hermes/db/eggnog \
    -i full_genome_proteins.faa -o emapper/full_genome --cpu 40 \
    --tax_scope 33090 --go_evidence non-electronic --override --temp_dir /tmp
```

Then build term2gene from full annotation → `clusterProfiler::enricher(gene=qtl_genes, universe=all_annotated_genes)`.

**Full-genome eggNOG notes**:
- Diamond runs `--sensitive --iterate` — hits file stays 0 bytes until ALL iterations complete (normal, not stuck)
- 30k queries × 40 threads ≈ 2-3 hours on CNS2
- Monitor with `ps -p PID -o etime,%cpu` (should show >2000% CPU for diamond)
- Run on CNS2 (more /tmp space for diamond temp files)

**Network pitfalls**:
- KEGG REST API (`rest.kegg.jp`) is BLOCKED → use pathway IDs as names
- Bioconductor is BLOCKED → install clusterProfiler from GitHub (`remotes::install_github("YuLab-SMU/clusterProfiler")`)
- go-basic.obo: download direct (no proxy), cache at `/media/nfs1/hermes/db/eggnog/`
- For proper enrichment, need full-genome eggNOG annotation as background (not just QTL genes)

---

# Part 4: Figure Standards

**ALL figures: 600 DPI minimum, BOTH PNG + PDF per figure.**

```python
DPI = 600
FIG_DIR = "04.gwas/figures/"

def save_both(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=DPI, bbox_inches='tight')
    fig.savefig(f'{FIG_DIR}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)
```

Required figures:
1. Per-trait Manhattan (all traits, all models) — **ALL SNPs must be plotted, no sampling. Point size 0.25 (larger than default). PNG only (no PDF needed for supplemental per-trait Manhattans due to file size).**
2. Per-trait QQ (all traits, all models)
3. Stacked multi-trait Manhattan (one per model)
4. Trait × chromosome QTL density heatmap
5. Colocalization hotspot scatter
6. Top QTL haplotype boxplots
7. GO/KEGG enrichment dotplot/bubble

---

# Part 5: HTML Report & Delivery

**Follow `bio-html-report-standards` skill for all report writing rules.** Key GWAS-specific additions:

## Report content checklist

- [ ] Species name verified (grep the report before delivering)
- [ ] Data sources: user-provided data marked "由用户提供", no fabricated pipeline details
- [ ] All 17 traits listed in a table with Chinese/English names and carbon chain features
- [ ] GWAS method explained in "what → why → how" three-step format
- [ ] Every figure has a "图 N（xxx）解读方法" block
- [ ] QQ plot included (one representative trait embedded, all 17 in bundle)
- [ ] Manhattan + QQ figures for ALL traits included in the bundle (even if only one embedded)
- [ ] Enrichment: ALL FDR<0.05 terms reported; cross-annotated redundancies merged in footnotes
- [ ] No subjective deletion of enrichment terms
- [ ] Discussion includes study limitations and future directions

## Delivery packaging

Package all figures and tables into a **single `report_bundle.tar.gz`**.

**Report naming convention**: Include dataset identifiers in filenames for quick identification:
- Pattern: `report_{N_samples}{id_type}_{N_traits}traits_{species}.html`
- Examples: `report_224srr_114traits_wandou.html`, `report_140num_29traits_wandou.html`, `report_146pv25_4traits_bean.html`
- Bundle follows same pattern: `report_224srr_114traits_wandou_bundle.tar.gz`
- Always provide absolute paths for both report and bundle.

```
XX.post_gwas/
├── report.html                    # self-contained HTML
└── report_bundle.tar.gz           # 附图/ + 附表/
    ├── 附图/
    │   ├── 图1_堆叠曼哈顿图.png/.pdf
    │   ├── 图2_QTL分布热图.png/.pdf
    │   ├── 图3_QQ图_代表性性状.png/.pdf
    │   ├── 图4_共定位热点分布.png/.pdf
    │   ├── 图5_单倍型箱线图.png/.pdf
    │   ├── 图6_GO_BP富集点图.png/.pdf
    │   ├── 图7_KEGG富集点图.png/.pdf
    │   ├── 补充_曼哈顿图_per_trait/
    │   ├── 补充_QQ图_per_trait/
    │   └── 补充_单倍型箱线图/
    └── 附表/
        ├── 表1_性状列表.txt
        ├── 表2_Top5共定位热点.tsv
        ├── 表2_完整_全N个共定位热点.tsv
        ├── 表3_Top5_QTL标记.txt
        ├── 表3_Top5_QTL标记_基因型矩阵.raw
        ├── 表4_富集_GO_BP.tsv
        ├── 表4_富集_KEGG.tsv
        └── 补充_QTL邻近基因*.tsv
```

**Packaging rules:**
- Figure naming: `图N_中文描述.png/.pdf` — N matches report figure number exactly
- Table naming: `表N_中文描述.tsv` — N matches report table number exactly
- All figures: 600 DPI PNG + vector PDF
- Supplement files in `补充_*/` subdirectories, never in root
- Absolute path for delivery: `/media/nfs1/hermes/project/<project>/XX.post_gwas/report_bundle.tar.gz`

## QQ plot generation

Fast3VmrMLM does NOT output QQ plots. Generate from midresult CSV.

**CRITICAL: Run on CNS2 for speed.** Reading 114 × 200MB midresult files over cross-node NFS from Python takes hours (0% CPU, stalled on I/O). R on CNS2 (where NFS is locally mounted) processes all 114 in ~5 minutes — **10× faster**.

```r
#!/usr/bin/env Rscript
# Batch QQ — run on CNS2 for fast local NFS access
suppressMessages(library(ggplot2))

mid_dir <- "fast3vmrmlm/results"
out_dir <- "figures"
files <- list.files(mid_dir, pattern = "_midresult\\.csv$", full.names = TRUE)

for (f in files) {
    trait <- gsub("_midresult\\.csv$", "", basename(f))
    dt <- read.csv(f, colClasses = c(pval = "numeric"))
    pvals <- dt$pval[is.finite(dt$pval) & dt$pval > 0 & dt$pval <= 1]
    if (length(pvals) < 100) next
    
    n <- length(pvals)
    obs <- -log10(sort(pvals))
    exp <- -log10((seq_len(n) - 0.5) / n)
    lam <- median(obs) / (-log10(0.5))
    
    df <- data.frame(Expected = exp, Observed = obs)
    p <- ggplot(df, aes(Expected, Observed)) +
        geom_point(size = 0.3, alpha = 0.3, color = "#333333") +
        geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.5) +
        labs(title = bquote(.(trait) ~ (lambda[GC] == .(round(lam, 4))))) +
        theme_bw(base_size = 11) + theme(panel.grid = element_blank())
    
    safe <- gsub("[/;: ]", "_", substr(trait, 1, 80))
    ggsave(file.path(out_dir, sprintf("qq_%s.png", safe)), p, width = 5, height = 5, dpi = 600)
    ggsave(file.path(out_dir, sprintf("qq_%s.pdf", safe)), p, width = 5, height = 5, dpi = 600)
}
```

Submit: `ssh CNS2 "nohup Rscript qq_all.R > qq.log 2>&1 &"`

**Performance comparison** (varies by dataset size):

| Dataset | Method | Speed |
|---|---|---|
| 114 traits × 200MB midresult (6.4M SNPs) | R on CNS2 | 19 QQ/min → ~5 min total |
| 4 traits × 370MB midresult (12M SNPs) | R on CNS2 | ~7 min/QQ → ~28 min total |
| Any size | Python cross-node NFS | **stalled** (0% CPU, I/O-bound) |

**Rule**: Always run QQ on CNS2 (local NFS). Per-trait time scales with midresult file size (~370MB = 7min, ~200MB = 3min).

---

# Part 6: Verification Checklist

After pipeline completion:
- [ ] rMVP: `N_traits × 3` CSV files in results/
- [ ] Fast3VmrMLM: `N_traits` midresult CSV + XLSX in results/
- [ ] λ_GC median near 1.0 for FarmCPU
- [ ] Hotspots found (≥1 hotspot with ≥3 traits) if traits are correlated
- [ ] All figures: 600 DPI PNG + PDF, English titles only
- [ ] QQ plots: generated for ALL traits from midresult CSV
- [ ] `report.html` opens standalone (all images base64-embedded)
- [ ] Correct species name throughout report
- [ ] All enrichment terms with FDR<0.05 reported in table
- [ ] Figure/table numbering matches between report and bundle
- [ ] `report_bundle.tar.gz` has 附图/ and 附表/ subdirectories
- [ ] No `__pycache__/` left in any directory

---

# Common Pitfalls (Cross-Tool)

| Issue | Tool | Fix |
|---|---|---|
| Using wrong VCF version | Pre-GWAS | Check all VCF in project vcf/, verify record counts, use only the FINAL filtered one (MAF+missing+GT) |
| IID mismatch (SRR IDs fine, numeric need zero-padding) | Both | Check if IDs are numeric; only pad if rMVP strips leading zeros |
| NFS read.table line count error | rMVP | Always use `sep=\"\\\\t\"` explicitly |
| NFS file.rename fails | rMVP | Use `outpath` in MVP() call — writes directly to target dir |
| NFS file.copy+file.remove loses files | rMVP | DO NOT use this pattern; use `outpath` instead |
| phe not all genotype ind | rMVP | Include ALL ind in phe_df, use NA for missing phenotypes |
| pheno not reordered to ind | rMVP | `idx <- match(ind, pheno$IID); pheno <- pheno[idx, ]` |
| BIM SNP names are \".\" | Both | Use `--chr X --from-bp Y --to-bp Y` instead of `--extract` |
| Consensus trait names use \".\" not \"_\" | Post-GWAS | Try `re.sub(r'\\.', '_', name)` + collapse underscores |
| Chinese chars in plots | Both | English-only plot titles |
| HTCondor not available | Both | ssh CNS2 + nohup |
| Rscript not found on SSH | Both | Full path `/media/nfs1/hermes/miniforge3/bin/Rscript` |
| rMVP built-in plots slow (>6M SNPs) | rMVP | Use pre-computed results if data unchanged; else accept 20+ min/trait |
| FarmCPU fails on binary | rMVP | Collapse to 0/1, drop small categories |
| TIFF too large for PIL | Fast3VmrMLM | `MAX_IMAGE_PIXELS = None` |
| QQ from midresult (Python cross-node) | Stalled at 0% CPU after minutes | Run R/ggplot2 on CNS2 directly — 10× faster (local NFS) |
| Enrichment: 0 results despite valid genes | Gene IDs lack `.1` suffix | Add `.1` to GFF gene IDs before enrichment |
| GFF chr names ≠ consensus chr names | 5 genes found near 1376 loci | Map `1→chr01, 2→chr02` before overlap |
| QQ from midresult over NFS | Reading 114×200MB files over NFS is infeasible | Run R on CNS2 (local NFS access) — 19 QQ/min |
| Fast3VmrMLM midresult col names | KeyError: 'Chromosome' or 'Position' in Python | Columns are `CHR,POS,MarkerID,Waldst,pval` (not `Chromosome,Position,RS#`). Always verify with `head -1 midresult.csv` before writing extraction scripts. |
| Manhattan only plotting sig SNPs | Figure looks empty — 4 dots for a trait with 4 sig SNPs | **Plot ALL SNPs** from midresult with `ax.plot(..., ',')` pixel marker; never plot only Bonf-sig subset |
| Manhattan mixed rendering (plot+scatter) | Sig SNPs look visually different from background | **ALL points must use the SAME rendering method** — `ax.plot`, not a mix of `plot` and `scatter`. Significance shown by position above threshold line, not by different marker style. |
| Manhattan sampling instead of all SNPs | User detects sparse/weird appearance | If user wants all points, use `ax.plot` with numpy arrays — 4M points renders fine with pixel marker + rasterized. |
| rm -rf deletes scripts too | After redo, R/Python scripts gone | Keep scripts in a `scripts/` dir at project root that survives cleanup. The `fast3vmrmlm/` and `post_gwas/` dirs get wiped on redo — scripts inside them vanish. |

## Reusing prior run results

When **input data is identical** (same genotype/phenotype/PCA), GWAS results are deterministic.
Copying prior results saves **days** of computation for large datasets (6M+ SNPs).

```bash
# Hardlink big.matrix and kinship (saves 20+ min)
ln prior_run/mvp.geno.{bin,desc,ind,map} new_run/rmvp/
ln prior_run/kinship.rds new_run/rmvp/

# Copy summary files
cp prior_run/summary/all_significant_snps.csv new_run/post_gwas/
cp prior_run/summary/cross_model_validated_loci.csv new_run/post_gwas/
cp prior_run/summary/method_trait_stats.csv new_run/post_gwas/summary.csv
```

The only non-deterministic step is Fast3VmrMLM (random seed), so always re-run that.

---

# Rejected Approaches (Don't re-add)

| Approach | Reason |
|---|---|
| HTCondor/SLURM | User explicitly disabled (2026-06-10) |
| mrMLM R package | Unusable at >10k markers, internal bugs |
| BLINK in rMVP 1.4.6 | Not available; only in rMVP ≥2.x |
| GAPIT | Not installed; adds unnecessary complexity |
| GEMMA | Different ecosystem; user prefers rMVP+Fast3VmrMLM |
| Per-chromosome QQ | λ_GC is genome-wide; per-chr not diagnostic |
| HapMap genotype format | Too slow for >10k markers; BEDMatrix only |
| pdf("/dev/null") to suppress rMVP plots | Plots still computed, just redirected; let rMVP plot normally then clean up devices |
| 4 large batches (28 traits each) | Graphics device leak kills traits after ~10/batch; use 10 smaller batches |


# NFS Pitfalls & Batch Strategy (6M+ SNPs)

1. **Always `sep="\t"`**: R's auto-sep-detection fails over NFS → "line 1 did not have N elements"
2. **Never `file.copy+remove` on NFS**: `file.remove` can execute before copy flush completes → data loss. Use `outpath` in MVP().
3. **10 batches × 12 traits**: Graphics device leak causes failures at ~10 traits per R session with larger batches
4. **`while(dev.list()>0) dev.off()`**: Clean up before/after every MVP() call

## Reusing prior results

When input data is **identical** (same genotype/phenotype/PCA), GWAS is deterministic. Reuse big.matrix + kinship + CSV summaries:

```bash
ln prior_run/mvp.geno.{bin,desc,ind,map} new_run/rmvp/
ln prior_run/kinship.rds new_run/rmvp/
cp prior_run/summary/*.csv new_run/post_gwas/
```

Fast3VmrMLM has random seed → always re-run.
| pdf("/dev/null") to suppress rMVP plots | Plots still computed, just redirected; let rMVP plot normally then clean up devices |

---

# NFS Pitfalls & Batch Strategy (Learned from Pea GWAS, 2026-06-12)

## NFS-specific issues

1. **Always use `sep="\t"` in read.table()**: R's auto-detection of separator fails silently when reading over NFS — works fine on local disk but returns "line 1 did not have N elements" over NFS.

2. **Never use `file.copy() + file.remove()` on NFS**: The `file.remove()` can execute before `file.copy()` finishes flushing, causing data loss. Always use `outpath=results_dir` in the MVP() call to write directly to the target directory.

3. **NFS `file.rename()` also fails silently**: Same root cause — use `outpath` instead.

## Batch sizing for 6M+ SNP datasets

| Batch strategy | Outcome |
|---|---|
| 4 batches × 28 traits | Too many traits per R session; graphics device leak causes failures after ~10 traits |
| 10 batches × 12 traits | ✅ Works reliably; each batch finishes within 60 min |
| 23 batches × 5 traits | Even safer but more SSH connections to manage |

For 6.4M SNP × 224 sample datasets: use **10 batches** (12 traits each). Each trait takes ~3-5 min (GLM+MLM+FarmCPU+built-in visualization).

## Graphics device management

rMVP v1.4.6 cannot suppress built-in plotting. After each trait, the graphics device may not clean up properly. Add before AND after every MVP() call:
```r
while (length(dev.list()) > 0) try(dev.off(), silent = TRUE)
```

## SSH batch submission

When submitting multiple `ssh CNS2 + nohup` jobs, SSH may timeout on later connections due to CNS2 load. Use a single SSH connection with a loop:
```bash
ssh CNS2 "
cd /path/rmvp
for i in 1 2 3 4 5 6 7 8 9 10; do
    nohup /media/nfs1/hermes/miniforge3/bin/Rscript run_rmvp.R \$i 10 > ../logs/batch\${i}.log 2>&1 &
    sleep 2
done
"
```

---

# Citations

## rMVP
Yin L, Zhang H, Tang Z, et al. **rMVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study.** *Genomics, Proteomics & Bioinformatics.* 2021;19(4):619-628. doi:10.1016/j.gpb.2020.10.007

## Fast3VmrMLM
Wang JT, Chen Y, Shu GP, et al. **Fast3VmrMLM: A fast algorithm that integrates genome-wide scanning with machine learning to accelerate gene mining and breeding by design for polygenic traits in large-scale GWAS datasets.** *Plant Communications* 2025;6(7):101385. doi:10.1016/j.xplc.2025.101385

## MLM
Yu J, Pressoir G, Briggs WH, et al. **A unified mixed-model method for association mapping that accounts for multiple levels of relatedness.** *Nature Genetics.* 2006;38(2):203-208. doi:10.1038/ng1702

## FarmCPU
Liu X, Huang M, Fan B, Buckler ES, Zhang Z. **Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies.** *PLoS Genetics.* 2016;12(2):e1005767. doi:10.1371/journal.pgen.1005767

## VanRaden K
VanRaden PM. **Efficient methods to compute genomic predictions.** *Journal of Dairy Science.* 2008;91(11):4414-4423. doi:10.3168/jds.2007-0980

## PLINK2
Chang CC, Chow CC, Tellier LCAM, et al. **Second-generation PLINK: rising to the challenge of larger and richer datasets.** *GigaScience.* 2015;4:7. doi:10.1186/s13742-015-0047-8

```bibtex
@article{Yin2021rMVP,
  title={{rMVP}: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study},
  author={Yin, Lilin and Zhang, Haohao and Tang, Zhenshuang and Xu, Jingya and Yin, Dong and Zhang, Zhiwu and Yuan, Xiaohui and Zhu, Mengjin and Zhao, Shuhong and Li, Xinyun and Liu, Xiaolei},
  journal={Genomics, Proteomics \& Bioinformatics}, volume={19}, number={4}, pages={619--628}, year={2021},
  doi={10.1016/j.gpb.2020.10.007}
}
@article{Wang2025Fast3VmrMLM,
  title={{Fast3VmrMLM}: A fast algorithm that integrates genome-wide scanning with machine learning to accelerate gene mining and breeding by design for polygenic traits in large-scale {GWAS} datasets},
  author={Wang, Jing-Tian and Chen, Ying and Shu, Guo-Ping and Zhao, Meng-Meng and Zheng, Ao and Chang, Xin-Yu and Li, Guo-Qing and Wang, Yu-Bo and Zhang, Yuan-Ming},
  journal={Plant Communications}, volume={6}, number={7}, pages={101385}, year={2025},
  doi={10.1016/j.xplc.2025.101385}
}
@article{Yu2006MLM,
  title={A unified mixed-model method for association mapping that accounts for multiple levels of relatedness},
  author={Yu, Jianming and Pressoir, Gael and Briggs, William H and others},
  journal={Nature Genetics}, volume={38}, number={2}, pages={203--208}, year={2006},
  doi={10.1038/ng1702}
}
@article{Liu2016FarmCPU,
  title={Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies},
  author={Liu, Xiaolei and Huang, Meng and Fan, Bin and Buckler, Edward S and Zhang, Zhiwu},
  journal={PLoS Genetics}, volume={12}, number={2}, pages={e1005767}, year={2016},
  doi={10.1371/journal.pgen.1005767}
}
@article{VanRaden2008K,
  title={Efficient methods to compute genomic predictions},
  author={VanRaden, P M},
  journal={Journal of Dairy Science}, volume={91}, number={11}, pages={4414--4423}, year={2008},
  doi={10.3168/jds.2007-0980}
}
@article{Chang2015PLINK2,
  title={Second-generation {PLINK}: rising to the challenge of larger and richer datasets},
  author={Chang, Christopher C and Chow, Carson C and Tellier, Laurent C A M and Vattikuti, Shashaank and Purcell, Shaun M and Lee, James J},
  journal={GigaScience}, volume={4}, pages={7}, year={2015},
  doi={10.1186/s13742-015-0047-8}
}
```

### Methods Template (Chinese)

> 全基因组关联分析使用 rMVP v1.4.6（Yin et al., 2021），整合了 GLM、MLM（Yu et al., 2006）和 FarmCPU（Liu et al., 2016）三种模型。亲缘关系矩阵基于 VanRaden（2008）方法计算。将 LD 修剪后基因型的前三个主成分作为固定效应协变量纳入模型。全基因组显著性阈值设为 Bonferroni 校正（α = 0.05 / 有效 SNP 数），建议性阈值设为 1/SNP 数。同时使用 Fast3VmrMLM v2.0（Wang et al., 2025）进行多位点 GWAS 分析作为补充，筛选标准为 LOD ≥ 3 且 P < 0.01。高置信度位点定义为至少在 2 种模型中均显著的 SNP 聚类结果。


# Part 7: Fast3VmrMLM Standalone Report

When the user wants a report focused ONLY on Fast3VmrMLM results (not rMVP), or when the Fast3VmrMLM XLSX files are empty:

## 7a. XLSX Emptiness Fallback

Fast3VmrMLM v2.0 default thresholds (LOD≥3, P<0.01, multi-locus EBayesEM) can be overly strict. For datasets where rMVP finds 1M+ Bonferroni-significant SNPs, Fast3VmrMLM XLSX files may ALL be empty (only "Trait ID" header, 0 QTL rows). **This is normal behavior for the multi-locus model with stringent defaults.**

**Fallback**: Extract significant SNPs directly from midresult CSVs using standard Bonferroni/suggestive thresholds:

```python
for f in midresult_files:
    with open(f) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pv = float(row['pval'])
            if pv < BONF:
                sig_snps.append({...})  # SIG
            elif pv < SUGGESTIVE:
                sug_snps.append({...})  # SUG
```

Then cluster SIG SNPs into loci (100kb), compute hotspots (500kb, ≥2 traits), and report λ_GC per trait.

## 7b. Standalone Report Structure

**Preferred approach: single Python script** that does QTL extraction → clustering → hotspots → figures → HTML report → bundle ALL in one go. This eliminates hours of iterative debugging from running separate scripts. Key sections:

**Preferred approach: single Python script** that does QTL extraction → clustering → hotspots → figures → HTML report → bundle ALL in one go. This eliminates hours of iterative debugging from running separate scripts. Key sections:

1. **分析概览** — Stats grid: total SIG SNPs, loci count, traits with QTL, hotspots, median λ_GC
2. **分析方法** — Fast3VmrMLM "是什么-为什么-怎么做" three-step, parameter table
3. **分析结果**:
   - 3.1 Summary table (top 15 traits by SIG count + λ_GC)
   - 3.2 Top QTL loci (chr, range, n_traits, n_snps, peak p-value)
   - 3.3 Hotspot table (chr, window, n_traits)
   - 3.4 Embedded figures (Manhattan top3, QTL heatmap, QQ, hotspot scatter)
   - 3.5 **Per-trait interpretation** — automated generation using signal strength, λ_GC evaluation, top locus coordinates, and candidate gene hints
4. **附图附表清单** — Match report figure/table numbers to bundle files
5. **软件与参考文献**
   - 性状遗传网络：从 loci 性状共现矩阵 → 热图 + 边列表
   - QTL 效应量排名：从 midresult Waldst 统计量 → 最大效应 SNP
5. **附图附表清单** — Match report figure/table numbers to bundle files
6. **软件与参考文献**

## 7d. Per-Trait Figures — User Preferences

### Manhattan Plots
**Plot ALL SNPs, never downsample.** User explicitly rejects sampling. Use:
- R/ggplot2 on CNS2 (local NFS for speed)
- `geom_point(size = 0.25, alpha = 0.4, stroke = 0, shape = 16)`
- **PNG only** (no PDF needed for per-trait Manhattan — files too large)
- `data.table::fread` for fast CSV reading
- 600 DPI, width=14, height=3.5

```r
dt <- data.table::fread(f, select = c("CHR","POS","pval"))
ggplot(dt, aes(x = cum_pos, y = neglogp, color = CHR)) +
    geom_point(size = 0.25, alpha = 0.4, stroke = 0, shape = 16) +
    scale_color_manual(values = setNames(colors, chr_order), guide = "none") +
    ggsave(out_png, width = 14, height = 3.5, dpi = 600)
```

Performance: 6.4M SNPs ≈ 2 min/trait, 16M SNPs ≈ 5 min/trait.

### QQ Plots
Same approach — all p-values, R on CNS2, PNG+PDF, size=0.3.

### Naming
Output: `manhattan_full_{trait}.png` for full Manhattan, placed in `附图/补充_全量曼哈顿图_per_trait/` in bundle.

**Candidate gene annotation**: Use eggNOG annotations + GFF to find genes near QTL peaks.

### Selection & Diversity for Selfing Species (Post-GWAS)

For selfing crops (pea, common bean, soybean), some standard population genetics analyses have limitations:

| Analysis | Viable? | Note |
|---|---|---|
| π (nucleotide diversity) | ✅ vcftools `--site-pi` | VCF-based π is inflated (~750× for pea, missing invariant sites). Useful for relative comparison between populations/regions but NOT absolute values. |
| π from BAM (ANGSD) | ✅ | Correct absolute π, but very slow (365 BAM × 4.5Gb genome). Submit to NFS and let run for days. |
| Ho/He | ❌ | Selfing species are homozygous — Ho≈0 everywhere, He≈0.5 everywhere. No signal. |
| Tajima's D | ⚠️ | Selfing bottleneck inflates negative D. Interpret cautiously. |
| Fst | ✅ vcftools `--weir-fst-pop` | Cross-population differentiation works. Need group labels. |
| π ratio (log2(π_A/π_B)) | ✅ | Reveals population-specific selective sweeps even with VCF-inflated π. |
| LD decay | ✅ PopLDdecay | Compare decay rates between populations — faster decay = higher diversity. |
| ROH | ❌ | Selfing makes long homozygous tracts the default, no selection signal. |
| iHS / XP-EHH | ⚠️ | Needs phased haplotypes (Beagle). Marginal value for selfing species. |

**Trait-trait network**: Count shared loci between all trait pairs → adjacency matrix → heatmap.
```python
for l in loci:
    for t1, t2 in combinations(l['traits'], 2):
        trait_pairs[tuple(sorted([t1,t2]))] += 1
```

**PVE ranking**: Use Waldst from midresult as effect-size proxy. Find max-Waldst SNP per trait.

**LD decay** (skip on 6M+ SNP datasets): plink2 `--r2` can produce empty vcor files on large datasets. Skip unless specifically needed.

**Figures required**:
- Stacked Manhattan (top 3 traits) from midresult CSVs (sample every 10th SNP for speed)
- QTL density heatmap: trait × chr matrix of SIG SNP counts
- Top QTL haplotype boxplots (plink2 chr:pos extraction → merge phenotype → matplotlib boxplot)
- Hotspot scatter plot (cumulative chr position vs N traits per window)

**Tables required**:
- 表1: trait summary (SIG/SUG counts, λ_GC)
- 表2: QTL loci list (peak SNP, range, n_traits)
- 表3: hotspots (chr, window, n_traits)

### Single-Script Pipeline

For new GWAS projects, consolidate ALL downstream analysis into **one Python script** that runs on CNS2:

```
build_pipeline.py:
  ├── Extract QTLs from midresult CSVs (streaming)
  ├── Cluster into loci (100kb), compute hotspots (500kb)
  ├── Generate all figures (Manhattan, heatmap, QQ, hotspot scatter)
  ├── Build self-contained HTML report (base64-embedded figures)
  ├── Write all TSV tables
  └── Package report_bundle.tar.gz
```

This eliminates hours of iterative debugging from running separate scripts. Just wait for Fast3VmrMLM to finish, then run one script.

### 16M+ SNP Dataset Notes

For datasets with >10M SNPs:
- First trait per batch takes 8-10 min (kinship computation + scanning)
- Subsequent traits ~4-6 min each
- Midresult files ~400-500MB each (vs ~200MB for 6M SNPs)
- Manhattan plot drawing often fails ("Drawing manplot failed") — benign, process continues
- Use **3 batches** for <30 traits, **6 batches** for 100+ traits
- Monitor with: `ls results/*_midresult.csv | wc -l`
- **Manhattan generation**: Use R on CNS2 with `data.table::fread`, sample to 800K points for plotting speed. ~45s per trait for 16M SNPs.
- **Processes stuck on "Drawing manplot"**: If midresult CSV exists, the computation is done — kill and proceed with results.

### Automated Per-Trait Interpretation

Generate concise biological interpretations for each trait based on GWAS statistics. This pattern produces meaningful text without manual curation:

```python
def interpret(trait, sig_count, lam, top_wald, top_chr, top_pos, n_loci):
    parts = []
    
    # Signal strength classification
    if sig_count == 0:
        parts.append('未检出Bonferroni水平显著关联')
    elif sig_count < 100:
        parts.append(f'检出{sig_count}个显著SNP，为寡基因遗传特征')
    elif sig_count < 1000:
        parts.append(f'检出{sig_count}个显著SNP，为中等复杂度遗传')
    else:
        parts.append(f'检出{sig_count}个显著SNP，提示多基因遗传或强LD')
    
    # λ_GC calibration
    if lam < 0.9:
        parts.append(f'λ_GC={lam:.3f}偏低，可能由群体结构过度校正导致')
    elif lam < 1.1:
        parts.append(f'λ_GC={lam:.3f}，模型校准良好')
    else:
        parts.append(f'λ_GC={lam:.3f}偏高，提示存在真实多基因信号或残留群体分层')
    
    # Top QTL location
    if top_wald > 0:
        parts.append(f'最强信号位于Chr{top_chr}:{top_pos:,}（Waldst={top_wald:.1f}）')
    
    if n_loci > 0:
        parts.append(f'聚类为{n_loci}个QTL位点')
    
    return '；'.join(parts)
```

Trait categories for classification: 抗逆性 (disease/stress), 形态与发育 (morphological), 产量相关 (yield), 物候期 (phenology), 其他.

### Single-Script Pipeline

For new GWAS projects, consolidate ALL downstream analysis into **one Python script** that runs on CNS2:

```
build_pipeline.py:
  ├── Extract QTLs from midresult CSVs (streaming)
  ├── Cluster into loci (100kb), compute hotspots (500kb)
  ├── Generate all figures (Manhattan, heatmap, QQ, hotspot scatter)
  ├── Build self-contained HTML report (base64-embedded figures)
  ├── Write all TSV tables
  └── Package report_bundle.tar.gz
```

This eliminates hours of iterative debugging from running separate scripts. Just wait for Fast3VmrMLM to finish, then run one script.

### 16M+ SNP Dataset Notes

For datasets with >10M SNPs:
- First trait per batch takes 8-10 min (kinship computation + scanning)
- Subsequent traits ~4-6 min each
- Midresult files ~400-500MB each (vs ~200MB for 6M SNPs)
- Manhattan plot drawing often fails ("Drawing manplot failed") — benign, process continues
- Use **3 batches** for <30 traits, **6 batches** for 100+ traits
- Monitor with: `ls results/*_midresult.csv | wc -l`
- **Manhattan generation**: Use R on CNS2 with `data.table::fread`, sample to 800K points for plotting speed. ~45s per trait for 16M SNPs.
- **Processes stuck on "Drawing manplot"**: If midresult CSV exists, the computation is done — kill and proceed with results.
