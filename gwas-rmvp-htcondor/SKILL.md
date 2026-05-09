---
name: gwas-rmvp-htcondor
description: Run multi-trait multi-model GWAS (GLM + MLM + FarmCPU) with rMVP on a HTCondor cluster. Handles PLINK BED → rMVP big.matrix conversion, VanRaden kinship, PC covariates, per-trait batching, adaptive λ_GC-driven PC tuning (sweep PC ∈ {0,1,2,3,5,7,10} and pick winner per (trait,model)), hardlink-based final_results consolidation to save hundreds of GB, and Python post-processing (Manhattan + QQ + top-hits summary) with English-only titles to avoid CJK font issues.
---

# GWAS with rMVP + HTCondor

Use this skill when the user wants GWAS on a SLURM/HTCondor cluster with many traits and wants all three of GLM, MLM, FarmCPU in one sweep, with Manhattan/QQ plots and top-hits tables.

## When to use

Triggers:
- "Run GWAS with rMVP"
- "GWAS 多性状多模型"
- User has a PLINK BED + phenotype TSV + PC covariate TSV
- Cluster has HTCondor (CNS1 Schedd in the wandou project) or SLURM (minor adjustment)

Do NOT use for:
- Small datasets that fit on a laptop (just run rMVP directly in R)
- Only one model needed and speed-critical → use GEMMA (see `gwas-gemma-slurm`)
- Binary case/control with rare variants → REGENIE

## Project layout (canonical)

Directory naming convention:
- Keep numeric prefixes (`01.ld_filter`, `02.pca`, ..., `05.gwas`) only for the **top-level pipeline stage directories** that define execution order
- Inside a single stage (e.g. `05.gwas/`) use plain names (`inputs`, `scripts`, `report`, `logs`) — no inner `00.` / `01.` prefixes
- Reports for a stage live inside that stage, e.g. `05.gwas/report/`, not in a project-global `00.report/`

```
<project>/05.gwas/
├── inputs/                           # source data, not written during run
│   ├── pheno_filtered.tsv            # FID IID trait_01 trait_02 ...  (IID zero-padded to 3 digits if fam uses that)
│   ├── sample_names.tsv
│   └── pheno_category_map.json       # optional: {original_to_safe: {"中文名":"trait_01", ...}, category_encodings: {...}}
├── scripts/                          # one-off prep scripts (not HTCondor entrypoints)
│   ├── prepare_pheno.py
│   └── 01_build_bed_pc.sh
├── logs/                             # top-level orchestration logs
├── report/                           # markdown summaries — LIVE IN THE STAGE, not a global 00.report/
│   ├── gwas_summary.md
│   └── scripts_usage.md
├── genotype_full_bed/
│   └── full_snps_NNN.{bed,bim,fam}   # full SNP set, NO LD pruning for GWAS
├── regenie_covariates/
│   └── quant_pc3.tsv                 # FID IID PC1 PC2 PC3  (computed from LD-pruned PCA)
├── rMVP/
│   ├── 01_prepare_data.R             # BED → big.matrix (one-shot)
│   ├── 02_run_rMVP.R                 # main GWAS driver
│   ├── run_batch.sh                  # HTCondor wrapper: Rscript 02_run_rMVP.R start end
│   ├── rmvp.condor                   # HTCondor submit file (vanilla universe)
│   ├── mvp.geno.{bin,desc,ind,map}   # produced by 01_prepare_data.R
│   ├── kinship.rds                   # VanRaden K, cached on first run
│   ├── logs/                         # HTCondor .log .out .err
│   ├── tune/                         # (optional) adaptive PC tuning artifacts
│   │   ├── runs/                     # trait_XX_pcNN/ per-combination CSVs — DELETE after final_results built (hardlink preserves data)
│   │   ├── state/                    # λ_GC per (trait, model, PC) — keep as audit trail
│   │   ├── compute_lambda.py
│   │   └── orchestrator.py
│   └── final_results/                # 78 hardlinks pointing into tune/runs/ — THE output for downstream
│       ├── trait_XX.{GLM,MLM,FarmCPU}.csv
│       ├── best_pc_table.tsv
│       └── _summary.tsv
└── post_gwas/
    ├── plot_rmvp_all.py              # Manhattan/QQ + summary tables
    ├── extract_all_bonf.py           # Bonferroni / suggestive hit extraction
    ├── plots/                        # trait_XX.{model}.{manhattan,qq}.png  (ENGLISH titles only, see pitfall)
    └── summary/
        ├── per_trait_model_summary.tsv    # n_p, min_p, lambda_gc, n_bonf, n_sugg
        ├── top10_per_trait_model.tsv
        ├── bonferroni_significant_hits.tsv
        ├── suggestive_hits.tsv
        └── model_overlap_by_trait.tsv
```

Symlink bridges (for scripts with hardcoded relative paths like `../pheno_filtered.tsv` in rMVP .R scripts):
```
05.gwas/pheno_filtered.tsv    -> inputs/pheno_filtered.tsv
05.gwas/sample_names.tsv      -> inputs/sample_names.tsv
05.gwas/pheno_category_map.json -> inputs/pheno_category_map.json
```
Cost zero space, keep scripts working without edits.

## Steps

### 1. Verify inputs

```bash
# BED samples match PC samples match pheno samples (IID column)
awk '{print $2}' <proj>/05.gwas/genotype_full_bed/full_snps_*.fam | sort > /tmp/fam.ids
tail -n +2 <proj>/05.gwas/pheno_filtered.tsv | awk '{print $2}' | sort > /tmp/pheno.ids
tail -n +2 <proj>/05.gwas/regenie_covariates/quant_pc3.tsv | awk '{print $2}' | sort > /tmp/pc.ids
comm -3 /tmp/fam.ids /tmp/pheno.ids   # expect empty
comm -3 /tmp/fam.ids /tmp/pc.ids      # expect empty
```

If IIDs don't match, zero-pad / strip leading zeros in the source until they do.

### 2. Convert BED → rMVP big.matrix (one-shot)

`01_prepare_data.R`:
```r
suppressPackageStartupMessages(library(rMVP))
bed_prefix <- "../genotype_full_bed/full_snps_NNN"
dir.create("data", showWarnings=FALSE)
MVP.Data(fileBed=bed_prefix, fileOut="mvp", verbose=TRUE)
```

Run:
```bash
cd <proj>/05.gwas/rMVP
/media/nfs1/hermes/miniforge3/bin/Rscript 01_prepare_data.R
```

Produces `mvp.geno.{bin,desc,ind,map}`. Takes 5–20 min for ~16M SNPs × 140 samples.

### 3. Main GWAS driver `02_run_rMVP.R`

Key parts (full template in `templates/02_run_rMVP.R`):

```r
args <- commandArgs(trailingOnly=TRUE)
trait_start <- as.integer(args[1]); trait_end <- as.integer(args[2])

geno <- attach.big.matrix("mvp.geno.desc")
map  <- read.table("mvp.geno.map", header=TRUE, stringsAsFactors=FALSE)
ind_raw <- readLines("mvp.geno.ind")
# CRITICAL: rMVP strips leading zeros → restore zero-padding to match fam/pheno
ind <- sprintf("%03d", as.integer(ind_raw))

pheno_all <- read.delim("../pheno_filtered.tsv", stringsAsFactors=FALSE,
                        colClasses=c("FID"="character","IID"="character"))
pheno_all <- pheno_all[match(ind, pheno_all$IID), ]
trait_names <- colnames(pheno_all)[-(1:2)]

pc_data <- read.delim("../regenie_covariates/quant_pc3.tsv", stringsAsFactors=FALSE,
                      colClasses=c("FID"="character","IID"="character"))
pc_data <- pc_data[match(ind, pc_data$IID), ]
CV <- as.matrix(pc_data[, c("PC1","PC2","PC3")])

if (file.exists("kinship.rds")) {
    K <- readRDS("kinship.rds")
} else {
    K <- MVP.K.VanRaden(geno, verbose=TRUE)
    saveRDS(K, "kinship.rds")
}

bonf <- 0.05 / nrow(map)
for (i in trait_start:min(trait_end, length(trait_names))) {
    tname <- trait_names[i]
    y <- as.numeric(pheno_all[[tname]])
    if (sum(!is.na(y)) < 50) next        # too few samples
    pheno_df <- data.frame(IID=ind, y); colnames(pheno_df)[2] <- tname
    tryCatch({
        MVP(
            phe=pheno_df, geno=geno, map=map, K=K,
            CV.MLM=CV, CV.FarmCPU=CV, nPC.GLM=3,      # PC covariates: PC1-3
            maxLoop=10, method.bin="static",
            threshold=bonf,
            method=c("GLM","MLM","FarmCPU"),
            ncpus=as.integer(Sys.getenv("OMP_NUM_THREADS","16")),
            file.output=TRUE, file.type="csv", outpath="results",
            verbose=FALSE
        )
    }, error=function(e) cat("ERR:", conditionMessage(e), "\n"))
}
```

### 4. HTCondor submission

`run_batch.sh`:
```bash
#!/bin/bash
set -eo pipefail
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
cd /abs/path/to/05.gwas/rMVP
/media/nfs1/hermes/miniforge3/bin/Rscript 02_run_rMVP.R "$1" "$2"
```

`rmvp.condor`:
```
universe       = vanilla
executable     = run_batch.sh
initialdir     = /abs/path/to/05.gwas/rMVP
getenv         = True
should_transfer_files = NO
request_cpus   = 16
request_memory = 64GB
request_disk   = 10GB
log            = logs/rmvp_$(Cluster)_$(Process).log
output         = logs/rmvp_$(Process).out
error          = logs/rmvp_$(Process).err
# batch traits by 5 (adjust for total N traits)
arguments = 1 5
queue
arguments = 6 10
queue
# ...
```

Submit:
```bash
cd <proj>/05.gwas/rMVP
mkdir -p logs results
condor_submit rmvp.condor
condor_q -submitter $USER
```

### 5. Monitor

```bash
condor_q
condor_tail <cluster>.<process>              # live stderr
tail -f logs/rmvp_0.err                      # batch 0 stderr
ls results/ | wc -l                          # should grow to N_traits × 3 (minus skipped)
```

Expected throughput: ~2–5 min per trait × model on 16 cores (15M SNPs, 140 samples). 29 traits × 3 models in 6 batches ≈ 1–3h wall-clock.

### 6. Post-processing (Python)

Use `templates/plot_rmvp_all.py` — reads every `results/trait_XX.MODEL.csv`, makes Manhattan + QQ, aggregates top hits. Run:

```bash
cd <proj>/05.gwas/post_gwas
mkdir -p plots summary
nohup /media/nfs1/hermes/miniforge3/bin/python plot_rmvp_all.py > plot_rmvp_all.log 2>&1 &
```

CSV p-value column is named `trait_XX.MODEL` (last column), NOT `p-value` or `P`. Read by position.

Requires: numpy, pandas, matplotlib, scipy (all available in miniforge3 base env).

## Parameters rationale

| Param | Value | Why |
|---|---|---|
| PC count | PC1–PC3 (default), tune up to PC10 if λ misbehaves | Matches 05.srr precedent. Use adaptive tuning (next section) when λ_GC > 1.15 |
| K method | VanRaden | rMVP default, robust for additive model |
| GLM PC | `nPC.GLM=3` | rMVP computes internal PCs from K; using 3 here matches external PCs |
| MLM PC | `CV.MLM=CV` | External PCs as fixed effects + K as random |
| FarmCPU PC | `CV.FarmCPU=CV` | Initial covariates; algorithm adds pseudo-QTNs iteratively |
| maxLoop | 10 | FarmCPU default, usually converges in 3–5 |
| method.bin | "static" | Faster than "FaST-LMM", adequate for ≤ 200 samples |
| Bonferroni | 0.05 / nSNP | Standard. For 16M SNPs ≈ 3.13e-09 |
| Suggestive | 1 / nSNP | For 16M SNPs ≈ 6.25e-08 |
| skip threshold | n_valid < 50 | Too few samples to fit reliably |
| IID padding | `sprintf("%03d", as.integer(x))` | rMVP strips leading zeros from fam; project uses 3-digit IDs |
| ncpus | 16 | Match request_cpus in condor file |
| request_memory | 64GB | 15M SNPs × 140 samples fits comfortably; raise to 96GB if >200 samples |
| λ_GC pass band | [0.85, 1.15] | Considered well-calibrated; outside this triggers PC tuning |

## Adaptive PC tuning (when λ_GC misbehaves)

When the initial PC=3 run produces λ_GC outside [0.85, 1.15] for many (trait, model) cells, sweep PC counts and pick the best per cell. Cuts dozens of inflated λ values without hand-tuning.

### Sweep design

- **PC test points**: `{0, 1, 2, 3, 5, 7, 10}` (skip 4, 6, 8, 9 — diminishing returns)
- **Selection rule per (trait, model)**: pick the PC count that brings λ_GC closest to 1.0 while staying inside [0.85, 1.15]. If no PC value lands in band, pick whichever has λ closest to 1.0 and flag with ⚠ in the report.
- **Output**: `final_results/best_pc_table.tsv` records chosen PC per cell.

### Hardlink trick (saves hundreds of GB)

Each (trait × PC) sweep produces 3 model CSVs in `tune/runs/trait_XX_pcNN/`. After picking the winner per cell, **hardlink** the chosen CSV into `final_results/`:

```bash
ln results/trait_XX_pc05/trait_XX.MLM.csv final_results/trait_XX.MLM.csv
```

Verify link count is 2 before deleting `tune/runs/`:
```bash
stat -c '%h %n' final_results/*.csv | awk '$1<2'   # any output = unsafe to delete
```

Once verified, `rm -rf tune/runs/` is safe — the inode survives via final_results/. On the wandou project this saved 332GB.

### Orchestrator pattern

`tune/orchestrator.py` should:
1. List (trait, PC) combinations remaining (skip done by checking `tune/state/`)
2. Submit HTCondor batches with env var `RMVP_NPC=<n>` and `RMVP_TRAIT_IDX=<i>`
3. After each batch, `compute_lambda.py` reads CSVs and writes `tune/state/trait_XX_pcNN.json` with λ per model
4. When all PC counts done for a trait, pick winner per model and hardlink into `final_results/`

`02_run_rMVP.R` reads `RMVP_NPC` to pick column count from the PC matrix:
```r
npc <- as.integer(Sys.getenv("RMVP_NPC", "3"))
CV  <- if (npc == 0) NULL else as.matrix(pc_data[, paste0("PC", 1:npc)])
```

## Interpreting outputs

`per_trait_model_summary.tsv` columns: trait, trait_zh, model, n_p, min_p, lambda_gc, n_bonferroni, n_suggestive.

Quality heuristics:
- λ_GC ≈ 1.0 ± 0.05 → model well-calibrated
- λ_GC > 1.2 → inflation (population structure not fully controlled)
- λ_GC < 0.9 → over-conservative (often MLM on small n)

Typical pattern on this project:
- GLM: fastest, may be inflated (λ up to 2.7 on polygenic traits)
- MLM: over-conservative, often 0 Bonferroni hits
- FarmCPU: most balanced → report FarmCPU as primary, use MLM/GLM for cross-check

Significant loci should ideally be supported by FarmCPU AND (MLM OR GLM with reasonable λ).

## Common pitfalls

### IID mismatch (most frequent)

Symptom: `MVP()` errors "sample names do not match" or all-NA output.
Cause: fam uses 3-digit zero-padded IIDs (`001`, `002`, ...), pheno uses the same, but rMVP silently strips leading zeros.
Fix: `ind <- sprintf("%03d", as.integer(readLines("mvp.geno.ind")))` — see step 3.

### CSV column name

Symptom: `header.index("p-value")` raises ValueError.
Cause: rMVP writes `trait_XX.MODEL` as the p-value column header.
Fix: use the last column (`df.columns[-1]` or `idx_p = len(header) - 1`).

### HTCondor evict

Symptom: batch shows "Normal termination" but started multiple times.
Cause: preemption on shared cluster. rMVP writes to disk atomically per trait, so partial progress is preserved and Condor auto-resubmits. Final results are correct.
Fix: none needed; just verify all `results/trait_XX.*.csv` are present after completion.

### Memory OOM on LTR-rich genomes

Symptom: Job goes on hold with `RequestMemory` exceeded.
Fix: `condor_qedit <JobId> RequestMemory 100000` (MB), `condor_release <JobId>`. For genomes >1Gb consider 120GB.

### Plots library

Symptom: `plot.type = NULL` but you still want plots.
Design: rMVP's built-in plotting is fine but inflexible. We use Python post-processing (`plot_rmvp_all.py`) for consistent styling and parallel rendering (168 PNGs in ~5 min on 16 cores).

### GLM inflated

Symptom: `λ_GC > 1.5` on many traits.
Diagnosis: residual population structure.
Options (in order of preference):
1. Run the adaptive PC tuning sweep (previous section) — lets each (trait, model) pick its own PC count
2. Raise to PC1–PC5 globally (edit `CV` matrix and `nPC.GLM=5`)
3. Check for cryptic relatedness in K heatmap
4. Drop GLM, keep MLM/FarmCPU only

Expected outcome after PC tuning: MLM ~95% pass rate, FarmCPU ~70%, GLM ~45% (GLM is fundamentally hard to calibrate on structured populations).

### Chinese characters in plot titles render as boxes/tofu

Symptom: matplotlib renders Chinese trait names as `□□□` or similar.
Root cause: default matplotlib font has no CJK glyphs.
**Preferred fix (Option A, recommended)**: strip Chinese from `title`, keep English-only:
```python
title = f'{trait_id} - {model}'   # e.g. "trait_01 - MLM"
```
Keep the Chinese name in the output TSV (`trait_zh` column) so downstream readers can cross-reference.

Why option A over "configure a CJK font": zero font dependency, plot files portable, no risk of font missing on reviewer machines. The Chinese↔English mapping lives in `pheno_category_map.json`.

### Disk explosion on adaptive tuning

Symptom: `05.gwas/` grows to 400+GB.
Cause: `tune/runs/` keeps every (trait × PC) sweep CSV (78 × 7 PC counts = up to 1638 CSVs × ~150MB each).
Fix: after `final_results/` is built via hardlinks, verify link count ≥ 2 on every final CSV, then `rm -rf tune/runs/`. State JSON in `tune/state/` is tiny — keep as audit trail.

### FarmCPU zero hits

Check `_signals.csv` too — FarmCPU's internal threshold may differ. Also inspect MAF distribution (rare variants filtered out often).

## Verification checklist

After the pipeline finishes:
- [ ] `ls final_results/*.csv | wc -l` == `N_traits_analyzed × 3` (exclude skipped traits)
- [ ] If using adaptive PC tuning: `stat -c '%h %n' final_results/*.csv | awk '$1<2'` returns nothing (all hardlinked)
- [ ] `per_trait_model_summary.tsv` has one row per (trait, model)
- [ ] `λ_GC` median near 1.0 for FarmCPU and MLM
- [ ] Plot folder has `2 × N × 3` PNG files, all with English-only titles (`grep -l '[\x{4e00}-\x{9fff}]' plots/*.png` returns nothing — but PNGs are binary, so just spot-check 2-3 with vision)
- [ ] HTCondor logs show no `ERR:` lines (re-run any error traits individually)
- [ ] Reports written to `<proj>/05.gwas/report/gwas_summary.md` (NOT `<proj>/00.report/`)
- [ ] `tune/runs/` deleted if adaptive tuning was used and hardlinks verified
- [ ] No `__pycache__/` left in scripts/post_gwas dirs (`find . -name __pycache__ -exec rm -rf {} +`)

## See also

- `gwas-gemma-slurm` — alternative single-model GWAS
- `population-genetics-vcf-analysis` — PCA / structure upstream of GWAS
- `slurm-cross-node-bioinformatics` — cluster troubleshooting
