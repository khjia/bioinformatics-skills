---
name: vcf-to-gwas
description: End-to-end pipeline from a VCF to a Fast3VmrMLM GWAS and per-trait HTML reports — PLINK QC, LD pruning, PCA, Fast3VmrMLM input files, multi-locus GWAS, nearest-gene QTN annotation, Manhattan/QQ, main-peak locus zoom, the full ±1 Mb gene list, regional QTN signals, haplotype boxplots, optional SuSiE fine-mapping with LD diagnostics, and a plant known-locus literature check via Europe PMC. Use when the user asks 从 VCF 做 GWAS, 跑 Fast3VmrMLM, 做 PCA/群体结构, 注释 QTN, 单性状 GWAS 报告, 主峰基因列表, 批量性状报告, 精细定位/SuSiE, or 查已知位点.
---

# VCF → GWAS → per-trait reports

Two phases in one skill:

- **A. VCF → GWAS**: PLINK QC, LD pruning, PCA, Fast3VmrMLM inputs, multi-locus
  GWAS, QTN annotation, gene coordinates.
- **B. Per-trait reports**: one self-contained HTML per trait with 7 figures
  and 12 tables.

Reference dataset: Yuanbaofeng, 200 samples, 2,494,579 SNPs on 13 chromosomes.

## Requirements

- `plink2` ≥ 2.0 (tested v2.0.0-a.6.9LM)
- R ≥ 4.2: `data.table`, `Fast3VmrMLM` (v2.0), `ggplot2`, `patchwork`,
  `ggrepel`, `knitr`, `rmarkdown`
- Python 3: `pandas`, `numpy`, `openpyxl`
- Pandoc (for HTML rendering), a CJK font (`Noto Sans CJK SC`)
- gene GFF3 (`final.gene.mRNA.longest.gff3`)

## Layout produced

```
<PROJECT>/
  01.phenotype/pheno_all.tsv
  02.ld_filter/beds/full_snps{,_renamed}.{bed,bim,fam}
  02.ld_filter/beds/ld_pruned.prune.in
  03.pca/ld_pruned_snps.{bed,bim,fam}
  03.pca/pca.eigenvec
  06.Fast3VmrMLM/inputs/{genotype.*, trait.csv, PCA.csv}
  06.Fast3VmrMLM/results/<trait>_{result.xlsx,midresult.csv,Manhattan_plot.tiff}
  06.Fast3VmrMLM/report/qtl_annotated.tsv
  09.postgwas2/
    data/gene_coords.rds
    scripts/{00_config.R, 20_trait_report_figs.R, 22_render_all.R,
             23_finemap_susie.R, 26_known_loci.py}
    report/{21_trait_report.Rmd, *.html}
    results/<prefix>_<trait>_report/*.tsv
    results/susie_<trait>/{credible_sets,pip,qtl_overlap}.tsv
    figures/<prefix>_<trait>_report/*.png
    figures/susie_<trait>/susie_pip.png
```

## Quick start

```bash
cp templates/pipeline.env.example pipeline.env     # edit PROJECT, filters, tools
set -a; source pipeline.env; set +a

# Phase A
bash templates/run_all.sh

# Phase B (after creating 09.postgwas2/scripts/00_config.R)
RUN_REPORTS=1 bash templates/run_all.sh            # or run only 22_render_all.R

# Phase C: fine-mapping (individual-level SuSiE, in-sample LD)
Rscript templates/23_finemap_susie.R fa fa_C24_1

# Phase D: known-locus lookup (Europe PMC)
python3 templates/26_known_loci.py --gene AT5G43760 --out known_kcs.tsv
```

Individual steps read the same env vars:

```bash
bash    templates/01_vcf2plink.sh
Rscript templates/02_make_inputs.R
Rscript templates/03_run_fast3vmrmlm.R
python3 templates/04_annotate_qtl.py
Rscript templates/05_gene_coords.R
```

# Phase A — VCF → GWAS

## 1. VCF → PLINK + PCA (`01_vcf2plink.sh`)

```bash
plink2 --vcf $VCF --allow-extra-chr --geno 0.2 --maf 0.05 --make-bed --out full_snps
plink2 --bfile full_snps --allow-extra-chr --set-all-var-ids '@:#:$r:$a' --make-bed --out full_snps_renamed
plink2 --bfile full_snps --allow-extra-chr --set-all-var-ids '@:#:$r:$a' \
       --indep-pairwise 50 10 0.2 --out ld_pruned
plink2 --bfile full_snps_renamed --allow-extra-chr \
       --extract ld_pruned.prune.in --make-bed --out 03.pca/ld_pruned_snps
plink2 --bfile 03.pca/ld_pruned_snps --allow-extra-chr --pca 10 --out 03.pca/pca
```

Reference run: 2,494,579 SNPs → 197,382 pruned SNPs for PCA.

## 2. Fast3VmrMLM inputs (`02_make_inputs.R`)

- `trait.csv`: first column header literally `<Phenotype>`, rows **in `.fam`
  order**, one column per trait.
- `PCA.csv`: first column header `<pca>`, then `pc_1..pc_NPC` (3 used).
- Rows must align with the genotype; Fast3VmrMLM matches by row, not ID.

## 3. GWAS (`03_run_fast3vmrmlm.R`)

```r
Fast3VmrMLM(fileGen=".../inputs/genotype", filePhe=".../trait.csv",
            filePS=".../PCA.csv", PopStrType="PC", fileOut=".../results/",
            genoType="SNP", trait=1:N, svrad=20, svpal=0.01, svmlod=3,
            nThreads=40, DrawPlot=TRUE, Plotformat="*.tiff")
```

Bonferroni: p = 0.05 / n_SNP ⇔ **LOD = 7.7** for 2,494,579 SNPs.
Per trait: `_result.xlsx` (LOD, add/dom, r2, P-value, significance),
`_midresult.csv` (single-marker scan), `_Manhattan_plot.tiff`.

## 4. QTN annotation (`04_annotate_qtl.py`)

Nearest gene within 50 kb (0 if inside a gene); farther QTN are kept as
`intergenic`. Output `qtl_annotated.tsv`:
`trait, marker, chr, pos, lod, significance, gene_id, gene_name, dist`.

## 5. Gene coordinates (`05_gene_coords.R`)

Builds `gene_coords.rds` (`gene_id, gene_name, chr, start, end, strand`) for
gene tracks and ±1 kb haplotype extraction.

# Phase B — Per-trait HTML reports

## Install and configure

1. Copy `templates/00_config.example.R` to
   `<PROJECT>/09.postgwas2/scripts/00_config.R` and edit paths, trait lists,
   Chinese names, `BONF_N`.
2. Copy the report assets:
   - `templates/20_trait_report_figs.R` → `<PROJECT>/09.postgwas2/scripts/`
   - `templates/21_trait_report.Rmd` → `<PROJECT>/09.postgwas2/report/`
   - `templates/22_render_all.R` → `<PROJECT>/09.postgwas2/scripts/`

## Run

```bash
export GWAS_WDIR=/path/to/09.postgwas2
export RSCRIPT=/path/to/Rscript

# one trait (figures, then render)
$RSCRIPT scripts/20_trait_report_figs.R fa fa_C24_1
$RSCRIPT -e 'rmarkdown::render("report/21_trait_report.Rmd",
  output_file="元宝枫_神经酸C24-1（fa_C24_1）_GWAS单性状报告_2026-10-26.html",
  output_dir="report",
  params=list(prefix="fa", trait="fa_C24_1", cn="神经酸C24:1"))'

# batch (2–3 min per trait)
setsid nohup $RSCRIPT scripts/22_render_all.R > logs_batch_report.log 2>&1 &
grep -c "OK$" logs_batch_report.log
```

## Report contents

| # | Figure | # | Table |
|---|--------|---|-------|
| 1 | Phenotype distribution | 1 | Dataset & phenotype |
| 2 | Manhattan (−log₁₀P, QTN marked) | 2 | Software versions |
| 3 | QQ (λ_GC) | 3 | Key parameters |
| 4 | QTN per chromosome | 4 | Overall statistics |
| 5 | Main-peak locus zoom + gene track | 5 | QTN per chromosome |
| 6 | Regional QTN LOD distribution | 6 | Top 20 QTN |
| 7 | Haplotype boxplots (top 3 genes) | 7 | Candidate genes |
|   |   | 8 | All genes in peak ±1 Mb |
|   |   | 9 | Regional QTN signals (50 kb pruning) |
|   |   | 10 | Haplotype results |
|   |   | 11 | Deliverables |
|   |   | 12 | All QTN |

# Phase C — Fine-mapping (optional)

Use after Phase A when a peak needs credible sets. Individual-level SuSiE with
in-sample LD is the default; summary-statistics mode needs an LD reference.

```bash
export GWAS_WDIR=/path/to/09.postgwas2
Rscript templates/23_finemap_susie.R fa fa_C24_1
# outputs: results/susie_fa_C24_1/{credible_sets,pip,qtl_overlap,summary}
#          figures/susie_fa_C24_1/susie_pip.png
```

- Region = anchor QTN ±1 Mb (`SUSIE_HALF_MB`), `L = 10` (`SUSIE_L`).
- Report `coverage = 0.9` at low N and say so; PIPs are rankings, not
  calibrated probabilities.
- Re-run with `L = 5` and `L = 15` to check stability.
- For summary statistics, run `susie_rss()` and `estimate_s_rss()` first;
  `lambda > 0.05` means the LD reference does not match the cohort.
- Full guidance and pitfalls: `references/fine-mapping.md`.

# Phase D — Known-locus lookup (optional, plant literature)

Before claiming novelty, map the candidate to its Arabidopsis/rice ortholog and
search the literature for that gene and the trait. No curated plant-wide
association catalog exists for this purpose.

```bash
python3 templates/26_known_loci.py --gene AT5G43760 \
    --trait "fatty acid seed oil" --out known.tsv --meta known_meta.tsv
```

The script queries Europe PMC, sorts by citation count, and writes a paper
table plus a metadata row (gene, trait, hit count). Treat hits as
hypothesis-generating evidence from model species, not as validation. See
`references/known-loci.md` for AraPheno, 1001 Genomes and Ensembl Plants.

## Method notes (keep consistent)

- **Anchor** = highest-LOD QTN with a gene annotation (fallbacks: highest-LOD
  QTN → strongest single marker).
- **Peak window** = anchor ±1 Mb → locus zoom, gene table, candidate genes.
- **Regional signals** = Fast3VmrMLM QTN in the window, greedily pruned at
  50 kb. Do *not* re-fit a single-marker conditional model: it is underpowered
  for most traits and returns zero signals.
- **SIG/SUG** = LOD ≥ 7.7 / 3 < LOD < 7.7.
- **Haplotypes** = top genes by QTN count → max LOD → distance; gene ±1 kb,
  tags p < 0.05, r² < 0.8, ≤4 tags, frequency ≥ 5%; KW + ANOVA.
- Never plot QTN LOD on a −log₁₀P axis. Figure 2 places QTN at their own
  single-marker −log₁₀P and says so in the subtitle.

## Verification

- Every HTML: 7 `<img>`, 12 `<table>`, zero `缺失图片`, no `hermes@`.
- Filenames: sanitise `:` → `-`, no spaces.
- Archive finished reports to `C:\Users\khjia\Pictures\` with project + date.

## Pitfalls

See `references/pitfalls.md`: encoding/BOM, Windows filename corruption from
tar/zip, `sprintf` and `%`, bare `+`/`-` table cells becoming markdown lists,
`results='asis'` for scroll tables, the `--geno` fallback, TBB/`.so`, and
ggrepel label placement.
