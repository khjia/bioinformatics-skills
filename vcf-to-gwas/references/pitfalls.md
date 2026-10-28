# Pitfalls — VCF → GWAS → reports

## Chromosomes and variant IDs

- Non-model genomes often use non-standard contig names (`chr01a`…`chr13a`).
  Always pass `--allow-extra-chr`; otherwise PLINK silently drops them.
- Rename variants once and reuse the key everywhere:
  `--set-all-var-ids '@:#:$r:$a'` → `chr04a:8780794:T:A`. Mixed ID styles are
  the most common cause of empty joins between GWAS output and annotation.
- Apply `--set-all-var-ids` before pruning/export when downstream code matches
  on IDs.

## Sample alignment

- Fast3VmrMLM reads `trait.csv` and `PCA.csv` **by row**, not by ID. Reorder
  both to `.fam` order and assert identity.
- Phenotype tables may use a different ID range than the genotype; match on
  IID, never on position.
- Header names matter: `<Phenotype>` for traits, `<pca>` for PCs.

## QC choices

- Reference run: `--geno 0.2 --maf 0.05` at VCF→bed.
- Regional extraction later uses `--maf 0.05` with a **looser** `--geno`
  fallback (0.1 → 0.25 → 0.5): some regions have 10–17% missing genotypes.
- LD pruning is two-stage: prune the full set, then `--extract prune.in
  --make-bed` for the PCA set. Do not prune again on the extracted set.

## Fast3VmrMLM

- The pre-compiled `.so` may need an older TBB. If `library(Fast3VmrMLM)`
  fails, run with `LD_PRELOAD=/path/to/libtbb.so.2`.
- `_result.xlsx` column names vary (`Position (bp)` vs `Position`,
  `significance` vs `Significance`). The annotation script matches by keyword;
  watch its skip messages.
- `_midresult.csv` is the single-marker scan; QTN LOD is multi-locus. They are
  not interchangeable.

## Annotation

- 50 kb window; inside-gene QTN get `dist = 0`; farther QTN stay
  `intergenic` so counts match the GWAS output.
- GFF `ID`/`Name` are parsed by regex; `Name` falls back to `ID`.

## R Markdown / knitr

- A `kable()` table wrapped in a scroll `<div>` needs a chunk with
  `results='asis'` and `cat(kable(...))`; otherwise the HTML is escaped.
- A table cell whose entire content is `+` or `-` is parsed by pandoc as a
  markdown list: it destroys the table and leaves `<div>` unclosed
  (`[WARNING] Div ... unclosed`). Use `正链`/`负链` or `—`.
- `sprintf("... ≥ 5% ...", x)` throws `invalid format '% 5'`. Use `paste0()`
  or `%%`.
- Inline R works in YAML (`date: "`r Sys.Date()`"`).

## Plotting

- **LOD ≠ −log₁₀P.** Fast3VmrMLM LOD comes from a multi-locus model; a QTN can
  have LOD 190 while its single-marker p is 1e-4. Never mix axes.
- 2.5M SNPs: bin at 1 kb (max −log₁₀P per bin) → ~416k points; point size
  1.3–1.4, alpha 0.7.
- Constant genotype columns make `cor()` warn; impute to the mean and set
  `NA` r² to 0.
- `ggrepel` with `direction = "y"` ignores `nudge_x`; lift labels with
  `nudge_y` (≈6).
- Do not re-fit a single-marker conditional model for regional signals: it is
  underpowered and returns 0 signals. Use the QTN list pruned at 50 kb.

## Batch rendering

- 22 traits take ~60–90 min. Use
  `setsid nohup Rscript scripts/22_render_all.R > logs_batch_report.log 2>&1 &`
  and poll with `grep -c "OK$"`.
- Check `attr(r, "status")` from `system2()`; a missing `.raw` aborts one trait
  but the loop continues — read the `FIG FAIL` lines.
- Verify 7 images / 12 tables / no `缺失图片` before delivery.

## Encoding and filenames (Windows ↔ Linux)

- Sources must be UTF-8 **without BOM**; a BOM gives
  `unexpected invalid token in "﻿"`.
- PowerShell `Set-Content -Encoding UTF8` adds a BOM; `Get-Content -Raw`
  misreads UTF-8 as ANSI. Use
  `[System.IO.File]::ReadAllText/WriteAllText(..., UTF8Encoding($false))`.
- Upload large files as base64 via `scp`, not inline SSH arguments.
- Do **not** extract Chinese-named tar/zip archives on a GBK Windows system:
  names become mojibake (`鍏冨疂…`, U+934F instead of `元` U+5143). Prefer
  per-file `scp`, or rename from file content afterwards.
- Sanitise `:` in trait names before building filenames (`C16:0` → `C16-0`).

## Delivery

- Archive reports to `C:\Users\khjia\Pictures\` with project + date in the name.
- Re-check filenames (character codes) and sizes after transfer.
