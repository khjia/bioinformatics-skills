---
name: gwas-rmvp-htcondor
description: Run multi-trait multi-model GWAS (GLM + MLM + FarmCPU + BLINK) with rMVP on a HTCondor cluster. Handles PLINK BED → rMVP big.matrix conversion, VanRaden kinship with genotype-fingerprint cache validation, PC covariates, per-trait batching with automatic retry, adaptive λ_GC-driven PC tuning (sweep PC ∈ {0,1,2,3,5,7,10} and pick winner per (trait,model)), hardlink-based final_results consolidation, multi-model consensus SNP integration (high-confidence loci), Python post-processing (per-trait Manhattan + QQ, multi-trait stacked Manhattan, trait×chrom hotspot heatmaps), and a single-file embedded HTML report. English-only plot titles to avoid CJK font issues.
---

# GWAS with rMVP + HTCondor

Use this skill when the user wants GWAS on a SLURM/HTCondor cluster with many traits and wants all four of GLM, MLM, FarmCPU, BLINK in one sweep, with Manhattan/QQ plots, multi-model consensus loci, and an embedded HTML report.

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
│   ├── kinship.fp                    # fingerprint "<size>_<mtime>_<md5>" of mvp.geno.bin — invalidates K on geno change
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
    ├── extract_all_bonf.py           # Bonferroni / suggestive hit extraction (run 1st)
    ├── extract_high_confidence.py    # multi-model consensus SNPs + locus clustering (run 2nd)
    ├── plot_rmvp_all.py              # per-trait Manhattan/QQ (run 3rd) — MODELS=[GLM,MLM,FarmCPU,BLINK]
    ├── plot_multitrait_summary.py    # stacked multi-trait Manhattan + trait×chrom heatmaps (run 4th)
    ├── generate_html_report.py       # single-file embedded HTML report (run 5th/last)
    ├── plots/                        # trait_XX.{model}.{manhattan,qq}.png  (ENGLISH titles only, see pitfall)
    ├── report.html                   # portable single-file report (base64-embedded PNGs)
    └── summary/
        ├── per_trait_model_summary.tsv    # n_p, min_p, lambda_gc, n_bonf, n_sugg (one row per trait×model)
        ├── top10_per_trait_model.tsv
        ├── bonferroni_significant_hits.tsv
        ├── suggestive_hits.tsv
        ├── model_overlap_by_trait.tsv
        ├── high_confidence_snps.tsv       # SNPs significant in ≥N models (default 2)
        ├── high_confidence_loci.tsv       # windowed cluster of consensus SNPs
        └── model_agreement_matrix.tsv     # pairwise model overlap across traits
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

Use `templates/02_run_rMVP.R` directly — it already bakes in three critical features that the naive driver lacks:

**A. Kinship fingerprint cache validation.** Naive `if (file.exists(\"kinship.rds\")) readRDS(...)` will happily load a stale K when the genotype file has changed underneath. This driver computes a fingerprint `"<size>_<mtime>_<md5>"` of `mvp.geno.bin` and stores it in `kinship.fp`. On next run, if the fingerprint mismatches, K is recomputed; otherwise the cached `kinship.rds` is reused. Invalidation is automatic and free on cache-hit.

**B. Trait-level retry (inside R, not shell).** Each trait is wrapped in a `tryCatch` + retry loop (default 3 attempts, 30s sleep + `gc()` between). Controlled via env `RMVP_RETRY`. **Why trait-level, not shell-level**: a shell-level `&& break` would rerun the whole batch on any single trait failure — wasting hours on the 20 traits that already succeeded. Retries at trait granularity preserve completed work; unrecoverable failures log `[FAIL]` and `next` to the following trait.

**C. Four-model sweep.** `method = c("GLM","MLM","FarmCPU","BLINK")`. BLINK is rMVP-built-in, costs ~10–20% extra wall-clock but adds an independent method for consensus filtering. No external deps.

Env vars honored by the driver:

| Var | Default | Purpose |
|---|---|---|
| `OMP_NUM_THREADS` | 16 | → `ncpus` in MVP() call |
| `RMVP_NPC` | 3 | Number of PC columns used as covariates (0 = none). Adaptive tuning sets this per-run. |
| `RMVP_RETRY` | 3 | Max attempts per trait before giving up |
| `RMVP_OUTDIR` | `results` | Output directory for CSVs (tuning sets this to `results/trait_XX_pcNN`) |

IID zero-padding (`sprintf("%03d", as.integer(ind_raw))`) is project-specific — adjust if your fam uses a different convention.

### 4. HTCondor submission

`run_batch.sh`:
```bash
#!/bin/bash
set -eo pipefail
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export RMVP_RETRY=${RMVP_RETRY:-3}         # trait-level retry count
export RMVP_NPC=${RMVP_NPC:-3}             # PC covariate count (adaptive tuning overrides)
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

Expected throughput: ~2–5 min per trait × model on 16 cores (15M SNPs, 140 samples). 4 models (GLM/MLM/FarmCPU/BLINK) × 29 traits in 6 batches ≈ 2–4h wall-clock.

### 6. Post-processing (Python) — 5-stage pipeline

Run all 5 scripts in order from `<proj>/05.gwas/post_gwas/`. All produce outputs under `summary/`, `plots/`, and `report.html`.

```bash
cd <proj>/05.gwas/post_gwas
mkdir -p plots summary
PY=/media/nfs1/hermes/miniforge3/bin/python
BASE=/abs/path/to/05.gwas                  # or the parent dir of final_results/

# 1. Bonferroni + suggestive hit extraction → summary/bonferroni_significant_hits.tsv etc.
$PY extract_all_bonf.py --base "$BASE"

# 2. Multi-model consensus (high-confidence) SNPs and loci
#    SNPs hit by ≥ --min-models get logged; neighbors within --window bp collapse to loci
$PY extract_high_confidence.py --base "$BASE" --min-models 2 --window 100000 --cutoff bonferroni

# 3. Per-trait Manhattan + QQ (one PNG pair per trait × model; 4 models now including BLINK)
$PY plot_rmvp_all.py

# 4. Multi-trait visualizations: stacked Manhattan + trait×chrom hotspot heatmap
#    --model picks which model's p-values to use for the stacked plot (FarmCPU recommended)
$PY plot_multitrait_summary.py --base "$BASE" --model FarmCPU --threshold suggestive

# 5. Single-file HTML report (base64-embeds all PNGs → emailable, no external assets)
$PY generate_html_report.py --base "$BASE"
# → produces report.html (typically 1–5 MB for ~30 traits)
```

Stages 1–4 are independent after stage 1 writes the summary TSVs. Stage 5 consumes everything. The HTML report is pure stdlib + pandas — no Jinja/Rmarkdown — portable to any browser.

CSV p-value column is named `trait_XX.MODEL` (last column), NOT `p-value` or `P`. Read by position.

Requires: numpy, pandas, matplotlib, scipy (all available in miniforge3 base env).

## Parameters rationale

| Param | Value | Why |
|---|---|---|
| Models | `c("GLM","MLM","FarmCPU","BLINK")` | 4-model sweep enables consensus filtering (see Multi-model consensus section). BLINK adds ~10–20% wall-clock. |
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
- FarmCPU: most balanced → report FarmCPU as primary, use MLM/GLM/BLINK for cross-check
- BLINK: similar spirit to FarmCPU but different pseudo-QTN selection; useful as an independent voter for consensus filtering

Significant loci should ideally be supported by FarmCPU AND at least one of {MLM, GLM with reasonable λ, BLINK}. See the high-confidence section below for automation.

## Multi-model consensus (high-confidence loci)

`extract_high_confidence.py` integrates the four model outputs to produce publication-grade consensus calls.

**Algorithm**:
1. For each (trait, model), collect SNPs passing `--cutoff` (default Bonferroni; `suggestive` or a raw float like `1e-6` also accepted).
2. Per SNP, count how many models called it; keep SNPs with count ≥ `--min-models` (default 2).
3. Cluster consecutive consensus SNPs within `--window` bp (default 100kb) on the same chromosome into loci. Each locus records lead SNP, supporting models union, and participating traits.

**Outputs** (in `summary/`):
- `high_confidence_snps.tsv`: one row per (trait, SNP) with `models_hit` (comma-joined), `n_models`, and min p across models
- `high_confidence_loci.tsv`: one row per locus (trait, chrom, start, end, lead_snp, lead_p, models_union, n_snps)
- `model_agreement_matrix.tsv`: pairwise model overlap counts across all traits — a quick health check (if GLM∩MLM is very low while GLM∩FarmCPU is high, GLM likely inflated)

**Tuning guidance**:
- For strict papers: `--min-models 3 --cutoff bonferroni`
- For discovery / QTL prioritization: `--min-models 2 --cutoff suggestive`
- `--window` should match your LD decay scale — 100kb is a reasonable default for selfing plants; widen to 500kb for outbreeders.

## Rejected optimizations (don't re-add these)

When reviewing this skill, the following suggestions were evaluated and **rejected with reason**. Don't re-add without new evidence.

| Suggestion | Reason rejected |
|---|---|
| `BiocParallel::MulticoreParam` in R | rMVP uses its own `ncpus` parameter; BiocParallel is not plumbed through and adds zero benefit. |
| `bigmemory::attach.big.matrix(..., backingcache="mmap")` | rMVP already memory-maps the geno via big.matrix backing file. Adding `backingcache="mmap"` doesn't change behavior. |
| Shell-level `while retry; do Rscript ... && break; done` | Reruns the ENTIRE batch on any single trait failure, wasting hours. Replaced with R-level `tryCatch` per trait. |
| Per-chromosome QQ plots | λ_GC is a genome-wide statistic; per-chrom QQ is not diagnostic and clutters the report. Genome-wide QQ is sufficient. |
| FASTLMM | Not in rMVP; separate tool. If needed, spin up a new skill. |
| `FarmCPU++`, `SUPERMODEL` | Not real tools (at time of writing) — likely hallucinated names. |
| Auto GFF annotation of lead SNPs | Kept as a future extension — too project-specific (species, GFF version, nomenclature) to bake into a general skill. Annotate downstream with `bedtools intersect` against your local GFF. |
| LD block visualization (Haploview-style) | Heavy dep (LDBlockShow / PLINK --blocks) for marginal gain in an HTML report. Can be added per-project if needed. |

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

### HTCondor memory hold

Symptom: Job goes on hold with `RequestMemory` exceeded.
Diagnosis for GWAS specifically:
- ~150 samples × 16M SNPs typically peaks at 30–50GB
- Kinship matrix is n×n (small for n<500), not the bottleneck
- Multiple parallel batches each mmap the big.matrix backing file — count copies × per-process working set
Fix: `condor_qedit <JobId> RequestMemory 80000` (MB), `condor_release <JobId>`. For >300 samples or running ≥3 parallel batches per node, request 96GB up front.

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
- [ ] `ls final_results/*.csv | wc -l` == `N_traits_analyzed × 4` (GLM/MLM/FarmCPU/BLINK; exclude skipped traits)
- [ ] If using adaptive PC tuning: `stat -c '%h %n' final_results/*.csv | awk '$1<2'` returns nothing (all hardlinked)
- [ ] `per_trait_model_summary.tsv` has one row per (trait, model) including BLINK
- [ ] `λ_GC` median near 1.0 for FarmCPU, MLM, and BLINK
- [ ] `kinship.fp` exists next to `kinship.rds` (fingerprint cache active)
- [ ] Plot folder has `2 × N × 4` PNG files (Manhattan + QQ × 4 models), all with English-only titles
- [ ] `summary/high_confidence_loci.tsv` exists and is non-empty on traits with known QTL
- [ ] `report.html` opens standalone in a browser (all PNGs base64-embedded)
- [ ] HTCondor logs show no lingering `[FAIL]` lines (traits that exhausted `RMVP_RETRY`). Re-run those individually with higher memory.
- [ ] Reports written to `<proj>/05.gwas/report/gwas_summary.md` (NOT `<proj>/00.report/`)
- [ ] `tune/runs/` deleted if adaptive tuning was used and hardlinks verified
- [ ] No `__pycache__/` left in scripts/post_gwas dirs (`find . -name __pycache__ -exec rm -rf {} +`)

## Citations (mandatory for any paper using this skill)

If you publish results produced by this pipeline, you **must** cite rMVP itself and the original method paper for every model you report. Reviewers in plant/animal genomics routinely check this.

**Must cite (the tool)**:

- Yin L, Zhang H, Tang Z, Xu J, Yin D, Zhang Z, Yuan X, Zhu M, Zhao S, Li X, Liu X. **rMVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study.** *Genomics, Proteomics & Bioinformatics.* 2021;19(4):619–628. doi:10.1016/j.gpb.2020.10.007

**Must cite (per model you report)**:

- **MLM** (Yu et al. 2006) — Yu J, Pressoir G, Briggs WH, et al. *A unified mixed-model method for association mapping that accounts for multiple levels of relatedness.* Nature Genetics. 2006;38(2):203–208. doi:10.1038/ng1702
- **FarmCPU** (Liu et al. 2016) — Liu X, Huang M, Fan B, Buckler ES, Zhang Z. *Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies.* PLoS Genetics. 2016;12(2):e1005767. doi:10.1371/journal.pgen.1005767
- **BLINK** (Huang et al. 2019) — Huang M, Liu X, Zhou Y, Summers RM, Zhang Z. *BLINK: a package for the next level of genome-wide association studies with both individuals and markers in the millions.* GigaScience. 2019;8(2):giy154. doi:10.1093/gigascience/giy154
- **GLM** — no single canonical reference; linear regression for GWAS is textbook. Cite the rMVP paper for the implementation.

**Must cite (kinship)**:

- **VanRaden K** (VanRaden 2008) — VanRaden PM. *Efficient methods to compute genomic predictions.* Journal of Dairy Science. 2008;91(11):4414–4423. doi:10.3168/jds.2007-0980

**Recommended (upstream / infrastructure)**:

- **PLINK 2** (if used for BED prep / LD pruning) — Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. *Second-generation PLINK: rising to the challenge of larger and richer datasets.* GigaScience. 2015;4:7. doi:10.1186/s13742-015-0047-8
- **HTCondor** (if reporting compute environment) — Thain D, Tannenbaum T, Livny M. *Distributed computing in practice: the Condor experience.* Concurrency and Computation: Practice and Experience. 2005;17(2–4):323–356. doi:10.1002/cpe.938

### BibTeX (copy-paste)

```bibtex
@article{Yin2021rMVP,
  title   = {{rMVP}: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study},
  author  = {Yin, Lilin and Zhang, Haohao and Tang, Zhenshuang and Xu, Jingya and Yin, Dong and Zhang, Zhiwu and Yuan, Xiaohui and Zhu, Mengjin and Zhao, Shuhong and Li, Xinyun and Liu, Xiaolei},
  journal = {Genomics, Proteomics \& Bioinformatics},
  volume  = {19}, number = {4}, pages = {619--628}, year = {2021},
  doi     = {10.1016/j.gpb.2020.10.007}
}
@article{Yu2006MLM,
  title   = {A unified mixed-model method for association mapping that accounts for multiple levels of relatedness},
  author  = {Yu, Jianming and Pressoir, Gael and Briggs, William H and others},
  journal = {Nature Genetics}, volume = {38}, number = {2}, pages = {203--208}, year = {2006},
  doi     = {10.1038/ng1702}
}
@article{Liu2016FarmCPU,
  title   = {Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies},
  author  = {Liu, Xiaolei and Huang, Meng and Fan, Bin and Buckler, Edward S and Zhang, Zhiwu},
  journal = {PLoS Genetics}, volume = {12}, number = {2}, pages = {e1005767}, year = {2016},
  doi     = {10.1371/journal.pgen.1005767}
}
@article{Huang2019BLINK,
  title   = {{BLINK}: a package for the next level of genome-wide association studies with both individuals and markers in the millions},
  author  = {Huang, Meng and Liu, Xiaolei and Zhou, Yao and Summers, Ryan M and Zhang, Zhiwu},
  journal = {GigaScience}, volume = {8}, number = {2}, pages = {giy154}, year = {2019},
  doi     = {10.1093/gigascience/giy154}
}
@article{VanRaden2008K,
  title   = {Efficient methods to compute genomic predictions},
  author  = {VanRaden, P M},
  journal = {Journal of Dairy Science}, volume = {91}, number = {11}, pages = {4414--4423}, year = {2008},
  doi     = {10.3168/jds.2007-0980}
}
@article{Chang2015PLINK2,
  title   = {Second-generation {PLINK}: rising to the challenge of larger and richer datasets},
  author  = {Chang, Christopher C and Chow, Carson C and Tellier, Laurent C A M and Vattikuti, Shashaank and Purcell, Shaun M and Lee, James J},
  journal = {GigaScience}, volume = {4}, pages = {7}, year = {2015},
  doi     = {10.1186/s13742-015-0047-8}
}
@article{Thain2005Condor,
  title   = {Distributed computing in practice: the {Condor} experience},
  author  = {Thain, Douglas and Tannenbaum, Todd and Livny, Miron},
  journal = {Concurrency and Computation: Practice and Experience}, volume = {17}, number = {2--4}, pages = {323--356}, year = {2005},
  doi     = {10.1002/cpe.938}
}
```

### Suggested Methods-section sentence

> Genome-wide association analysis was performed using rMVP v1.x (Yin et al., 2021), which implements GLM, MLM (Yu et al., 2006), FarmCPU (Liu et al., 2016), and BLINK (Huang et al., 2019). The genomic relationship matrix was computed following VanRaden (2008). The first three principal components from LD-pruned genotypes were included as fixed-effect covariates. Genome-wide significance was set at the Bonferroni threshold (α = 0.05 / N_SNP) and a suggestive threshold of 1/N_SNP. High-confidence loci were defined as SNPs significant in ≥ 2 of the 4 models, clustered within 100 kb windows.

## See also

- `gwas-gemma-slurm` — alternative single-model GWAS
- `population-genetics-vcf-analysis` — PCA / structure upstream of GWAS
- `slurm-cross-node-bioinformatics` — cluster troubleshooting
