---
name: vcf-to-gwas
description: End-to-end pipeline from a VCF to a Fast3VmrMLM GWAS and per-trait HTML reports — PLINK QC, LD pruning, PCA, Fast3VmrMLM input files, multi-locus GWAS, nearest-gene QTN annotation, then Manhattan/QQ, main-peak locus zoom, the full ±1 Mb gene list, regional QTN signals, and haplotype boxplots. Use when the user asks 从 VCF 做 GWAS, 跑 Fast3VmrMLM, 做 PCA/群体结构, 注释 QTN, 单性状 GWAS 报告, 主峰基因列表, or 批量性状报告.
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
    scripts/{00_config.R, 20_trait_report_figs.R, 22_render_all.R}
    report/{21_trait_report.Rmd, *.html}
    results/<prefix>_<trait>_report/*.tsv
    figures/<prefix>_<trait>_report/*.png
```

## Quick start

```bash
cp templates/pipeline.env.example pipeline.env     # edit PROJECT, filters, tools
set -a; source pipeline.env; set +a

# Phase A
bash templates/run_all.sh

# Phase B (after creating 09.postgwas2/scripts/00_config.R)
RUN_REPORTS=1 bash templates/run_all.sh            # or run only 22_render_all.R
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

## Citations

**Fast3VmrMLM (GWAS engine)**
Wang J, Chen Y, Shu G, Zhao M, Zheng A, Chang X, Li G, Wang Y, Zhang Y-M.
Fast3VmrMLM: a fast algorithm that integrates genome-wide scanning with
machine learning to accelerate gene mining and breeding by design for
polygenic traits in large-scale GWAS datasets. *Plant Communications*.
2025;6(7):101385. doi:10.1016/j.xplc.2025.101385

**3VmrMLM (underlying multi-locus model)**
Li M, Zhang Y-W, Zhang Z, et al. A compressed variance component mixed model
for detecting QTNs and QTN-by-environment and QTN-by-QTN interactions in
genome-wide association studies. *Molecular Plant*. 2022;15(4):630–650.
doi:10.1016/j.molp.2022.02.012

**PLINK**
Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ.
Second-generation PLINK: rising to the challenge of larger and richer
datasets. *GigaScience*. 2015;4:7. doi:10.1186/s13742-015-0047-8

**R packages**
Wickham H. *ggplot2: Elegant Graphics for Data Analysis*. 2nd ed. Springer; 2016.
Dowle M, Srinivasan A. *data.table: Extension of data.frame*. R package.
Xie Y. *Dynamic Documents with R and knitr*. 2nd ed. Chapman & Hall/CRC; 2015.
Xie Y, Allaire JJ, Grolemund G. *R Markdown: The Definitive Guide*. Chapman & Hall/CRC; 2018.

```bibtex
@article{wang2025fast3vmrmlm,
  title   = {Fast3VmrMLM: a fast algorithm that integrates genome-wide scanning with machine learning to accelerate gene mining and breeding by design for polygenic traits in large-scale GWAS datasets},
  author  = {Wang, Jingtian and Chen, Ying and Shu, Guoping and Zhao, Miaomiao and Zheng, Ao and Chang, Xiaoyu and Li, Guiqi and Wang, Yibo and Zhang, Yuan-Ming},
  journal = {Plant Communications},
  volume  = {6}, number = {7}, pages = {101385}, year = {2025},
  doi     = {10.1016/j.xplc.2025.101385}
}
@article{li2022compressed,
  title   = {A compressed variance component mixed model for detecting QTNs and QTN-by-environment and QTN-by-QTN interactions in genome-wide association studies},
  author  = {Li, Meiyue and Zhang, Ya-Wen and Zhang, Zhen and others},
  journal = {Molecular Plant},
  volume  = {15}, number = {4}, pages = {630--650}, year = {2022},
  doi     = {10.1016/j.molp.2022.02.012}
}
@article{chang2015plink,
  title   = {Second-generation PLINK: rising to the challenge of larger and richer datasets},
  author  = {Chang, Christopher C and Chow, Carson C and Tellier, Laurent CAM and Vattikuti, Shashaank and Purcell, Shaun M and Lee, James J},
  journal = {GigaScience},
  volume  = {4}, pages = {7}, year = {2015},
  doi     = {10.1186/s13742-015-0047-8}
}
@book{wickham2016ggplot2,
  title     = {ggplot2: Elegant Graphics for Data Analysis},
  author    = {Wickham, Hadley},
  edition   = {2nd},
  publisher = {Springer},
  year      = {2016}
}
```

**Methods template**

"Genome-wide association analysis was performed with Fast3VmrMLM v2.0
(Wang et al., 2025) on 2,494,579 biallelic SNPs, using the first three
principal components as covariates. Variants with LOD ≥ 7.7 (Bonferroni
p = 2.0 × 10⁻⁸) were declared significantly associated; 3 < LOD < 7.7 were
considered suggestive. Regional signals were summarised from the multi-locus
QTN output after 50-kb pruning, and haplotype association used tag SNPs
(p < 0.05, r² < 0.8) tested by Kruskal-Wallis and ANOVA."
