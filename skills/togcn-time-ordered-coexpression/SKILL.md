---
name: togcn-time-ordered-coexpression
description: Time-Ordered Gene Coexpression Network (TO-GCN) construction from time-series RNA-seq (root/tissue development). The SAME toolchain covers BOTH single-group (Single Time Series, one condition) AND double-group (Comparative, two conditions — often the user copies one condition into two to reuse their familiar peanut workflow). Covers iTAK TF annotation, TPM matrix prep (replicate-merged OR replicate-inclusive), single-seed BFS level assignment, PCC threshold selection for target network size (8k-10k genes), user's Mean/MedianExpHeatmap standard (raw/log2/row-scale/row-scale+numbers), and the accepted matplotlib two-circle concentric-ring network layout (inner=TF, outer=structural). Critical bug fixes (stack smashing from long gene IDs + PCC NaN), and the C++ BFS performance workaround (reproduce with R/igraph O(V+E)).
triggers:
  - TO-GCN
  - time-ordered gene coexpression network
  - 时间有序共表达网络
  - developmental time series coexpression
  - root development GCN
  - seed genes for coexpression
---

# TO-GCN: Time-Ordered Gene Coexpression Network

Construct a time-ordered gene coexpression network from **single-condition time-series** RNA-seq (e.g. root development root1幼→root6老). The pipeline assigns a time-order (BFS level) to each gene in a TF×gene coexpression network, revealing hierarchical regulatory cascades over development.

## When to use
- User has a **time-series** (multiple sample points over time), NOT a tissue/vs-control comparison
- Wants to infer regulatory hierarchy (which TF upstream/downstream) from coexpression over development
- Has single condition (no two-treatment comparative design) → use **Single Time Series** mode

## Reference project
- 珊瑚菜 (Glehnia littoralis) root development: `/media/nfs2/hermes/yanda/shanhucai/TO-GCN/`
- Data: 5 time points (root1-4 × 3 reps, root6 × 2 reps), TPM matrix 23,233 genes

## Two modes (choose correctly)
1. **Comparative** (main dir `Cutoff/GCN/TO-GCN`): needs two conditions × time (n1+n2 columns). 8 coexpression types.
2. **Single Time Series** (`Single_Time-series_data/source_code/` `*_single`): one condition only. C1+ positive coexpression only. **Use this for single-condition development.**

## USER GUIDANCE: BOTH single-group AND double-group use THIS SAME toolchain
The user explicitly stated (2026-08): **"单组还是双组流程都是用这套流程"** — single-group and double-group TO-GCN both run through the same Cutoff → GCN → TO-GCN(level) → heatmap → network pipeline. The difference is ONLY:
- **Single**: `Cutoff_single/GCN_single/TO-GCN_single`, input = N time-point columns (or N×reps), one condition.
- **Comparative (double)**: `Cutoff/GCN/TO-GCN` main binaries, input = n1+n2 columns. When the user has only ONE condition but prefers the Comparative syntax, **duplicate the data** (copy columns twice) → `C1=C2` → `C1+C2+` ≡ single-condition positive coexpression. This is the user's preferred route (familiar peanut workflow).
So there is ONE consolidated skill; do not create separate single vs double skills.

## KEY: Comparative mode with COPIED data = single-condition (user's preferred approach)
A user familiar with Comparative mode may ask to "copy the data" and run `GCN n1 n2` where both conditions are the same data:
- **C1=C2 (same columns) → R_LD == R_TD exactly** → C1+C2+ = single-condition positive coexpression (equivalent to Single mode C1+). Other 6 types (C1+C20 etc.) come out EMPTY (67-byte headers) — expected, don't worry.
- This reproduces the user's familiar peanut workflow. Data prep: make input files with n1+n2 columns (copy the n1 columns twice). Command: `GCN n1 n2 TF.txt gene.txt pos1 pos2 neg1 neg2`.
- **Sample-count pitfall**: user's "18" in peanut = 6 timepoints × 3 reps (replicates NOT merged, each rep a column). For single-condition dev data, use ALL raw replicates as columns (e.g. root1-4×3 + root6×2 = 14) → `GCN 14 14`. Replicates being uneven (root6 only 2) is NOT a problem for PCC — each rep is just a data point.
- Cutoff on 14 samples gives LOWER thresholds (~0.74/-0.63) than on 5 merged timepoints (~0.90/-0.86) because more samples = more df = lower critical PCC. ALWAYS re-run Cutoff after changing n1/n2.

## Network size control (target 8k-10k genes)
Genome-wide coexpression is dense; gene count drops SLOWLY with PCC threshold (hubs persist). 14-sample 0.93 → 8,630 genes / 83K edges. Scan `awk -F',' 'NR>1 && $4>=T'` over C1+C2+.csv to pick threshold hitting the target gene count. A single high threshold (e.g. 0.93) works better than trying to filter post-hoc.

## Seed selection: use ONE seed gene (not many)
TO-GCN official + peanut workflow use a **single representative seed** per BFS run. Multiple seeds flatten the level structure (139 seeds → 19 levels, low resolution). Pick a seed that is: (1) a TF, (2) high network degree (hub, covers most nodes), (3) matches the seed expression criterion (first-timepoint-high, declining). E.g. Glh1_09G003566 (AP2/ERF-ERF, degree 174, root1=114→root6=9). One seed → clean ~20 levels over 8k genes.

## Heatmap standard (user's MeanExpHeatmap.r / MedianExpHeatmap.r)
Follow the user's proven heatmap script style: for EACH of mean and median, produce ONE PDF with 4 panels:
1. raw data (scale="none", title "none centered&scaled")
2. log2(data+1) (scale="none")
3. row-scaled z-score (scale="row", "row centered&scaled")
4. row-scaled + display_numbers
Same palette (blue-white-red), cluster_rows=F, cluster_cols=F, angle_col=45.
Data matrix = level × time-point (levels from Node_level, columns = time-point mean TPM).

## Source + binaries
Repo: `https://github.com/petitmingchang/TO-GCN` (clone via proxy: `git -c http.proxy=http://127.0.0.1:10811 clone --depth 1 ...`)
Linux binaries NOT shipped for Single mode (only Windows/MacOS) → **compile from source**:
```bash
g++ Cutoff_single.cpp -o Cutoff_single
g++ GCN_single.cpp -o GCN_single
g++ TO-GCN_single.cpp -o TO-GCN_single
```

## Input formats (all TSV, NO header, tab-separated)
- Gene expression file: each line `geneID<TAB>val1<TAB>val2<...>valN` where N = # time points
- TF genes file: same format (subset of all genes with expression values)
- All genes file: same format (all genes)
- seeds.txt: one gene ID per line

## Data preparation (critical)
1. **Collapse replicates to time-point means**: root1-4 avg of 3 reps, root6 avg of 2 reps.
2. **UNIFY GENE IDs — THE #1 CRASH CAUSE.** TO-GCN uses `char GID[20]` buffers in `fscanf("%s", GID)`. Any gene ID >19 chars causes **stack smashing detected: terminated**. Plant annotation files often mix `Glh1_01G000001` and `gene.Glh1_01G000456.mRNA1` (25 chars). **Clean ALL IDs to the short form**: strip `gene.` prefix and `.mRNA\d+` suffix. Verify `max(len(ID)) <= 19` before running.
3. **TPM is fine** — user confirmed TPM works directly (no need to convert to RPKM).
4. iTAK TF ids use clean IDs (strip `.mRNA1` from proteins.fa headers).

## Pipeline steps
```
Step0 Data prep: 00.data/ All_genes_expression.tsv, TF_genes_expression.tsv, seeds.txt
Step1 Cutoff_single  N_features All_genes_expression.tsv  → pos+neg cutoffs (p<0.05)
Step2 GCN_single     N_features TF_genes_expression.tsv All_genes_expression.tsv <pos_cutoff>
Step3 TO-GCN_single  N_features All_genes_expression.tsv seeds.txt <pos_cutoff>
```

### Step 1: Cutoff_single
```bash
Cutoff_single 5 All_genes_expression.tsv   # outputs suggested pos/neg cutoffs
```
Use ALL genes file (gives distribution of all TF×gene PCCs). On 珊瑚菜 5 time points → pos≈0.91, neg≈-0.86.

### Step 2: GCN_single
```bash
GCN_single 5 TF_genes_expression.tsv All_genes_expression.tsv 0.91
```
TF × all-gene pairs, PCC ≥ cutoff. Outputs `C1+.csv` (TF gene ID, coexpression type, gene ID, PCC). **Computationally heavy**: 1869 TF × 23233 genes ≈ 43M PCC calcs → run background. 珊瑚菜 produced 1.67M edges at 0.91, **798K edges at 0.95**. **Use a cutoff that keeps the network tractable for BFS: prefer ~0.95 (≈800K edges) over 0.91 (1.67M).** At 0.95, node coverage barely drops (23,172 vs 23,233 nodes).

### Step 3: TO-GCN_single
```bash
TO-GCN_single 5 All_genes_expression.tsv seeds.txt 0.91
```
BFS assigns time levels. Outputs `Node_level.csv` (TF gene ID, level in GCN) + `Node_relation.csv` (node_1, node_2, PCC). **Can be slow on large (1.6M+ edge) networks** — run background, no timeout. If too slow, raise the PCC cutoff to shrink the network.

### CRITICAL: C++ BFS is O(levels×edges×TFs) — reproduce with R/igraph instead
On genome-wide networks (e.g. 798K edges, 1869 TFs) the C++ TO-GCN_single BFS is **impractical** — ran 7-20+ min with no output (RSS grows to ~560MB but never finishes, `set_neighbor_pos` scans ALL edges per level AND for each match scans all TFs). **Do NOT wait for it.** Reproduce the identical BFS in R/igraph using native C (O(V+E), seconds):

```r
library(igraph)
g <- graph_from_data_frame(edges[,c("tf","gene")], directed=FALSE)  # edges = C1+.csv (tf,gene,pcc)
d <- distances(g, v=which(V(g)$name %in% seeds), to=V(g), mode="all")  # shortest dist from any seed
min_d <- apply(d, 2, min); min_d[is.infinite(min_d)] <- -1
level <- ifelse(min_d >= 0, min_d + 1, NA)
level[! (V(g)$name %in% tf_set)] <- NA     # ONLY TF genes get a level (TO-GCN behavior)
# Write Node_level.csv: gene, level_in_GCN (==1 kept as 1, else level-1 per TO-GCN convention)
```

This gives the **exact same result** as the C++ BFS (both = BFS shortest-path distance from seeds), in seconds. Full working script: `reproduce_togcn_bfs.R` in the reference project's `scripts/`. Also, seeds that are NOT TF genes still anchor BFS (their neighbors get levels), but only TF nodes receive recorded levels.

### Making ALL genes get levels (TF→gene regulatory network)
By default TO-GCN only assigns levels to genes in the FIRST input file (the "TF list" / `TF_exp_table`). To build a **TF→gene regulatory network** where every gene gets a level, **change the input so all genes are nodes**: in TO-GCN_single pass `gene.txt gene.txt` (all-genes expression file) instead of `TF.txt`. In the R/igraph reproduction, simply **drop the `level[!tf] <- NA` line** so ALL genes receive levels:
```r
level <- ifelse(min_d >= 0, min_d + 1, NA)   # ALL genes get a level (not just TF)
```
珊瑚菜 result (all-genes): 23,172 genes, 11-12 levels, still shows the same early→late biology (L1-L4 early root6/root1<1, L6-L11 late root6/root1 up to 3.5).

Per-level gene/TF counts (all-genes BFS, 珊瑚菜, 0.95 cutoff): L1:3148/260TF, L2:4657/244, L3:3147/197, L4:2321/168, L5:1928/164, L6:1988/178, L7:1941/193, L8:2601/286, L9:1397/176, L10:43/3, L11:1/0. Compute via `awk -F',' 'NR>1{print $3}' Node_level_allgenes.csv | sort -n | uniq -c` (level in col 3 after gene,col2 check).

## Comparative mode — duplicating data does NOT change network
Users familiar with the double-sample (Comparative) workflow may ask to duplicate the single-condition data into two copies (C1=C2). **Key facts:**
- Comparative `GCN 18 18 TF.txt gene.txt` computes PCC under C1 (first n1 cols) and C2 (last n2 cols).
- If C1=C2 (same data duplicated), then R_LD == R_TD exactly, so **`C1+C2+` reduces to single-condition positive coexpression** — identical to Single mode `C1+`.
- The other 7 types (C1+C20, C1-C2- etc.) either can't occur (a correlation can't be positive AND uncorrelated on identical data) or collapse onto the same classes.
- **IMPORTANT: duplicating data does NOT change network size** — you still get ~23K genes / 80万 edges. It only lets you use the familiar Comparative command syntax (`Cutoff/GCN/TO-GCN` main binaries with `n1 n2` args). The input file must have n1+n2 columns (duplicate the 5 time points → 10 columns).
- Use Comparative mode mainly if the user's workflow/tools (e.g. MFSelector seed selection, specific downstream scripts) expect it; it's NOT a way to shrink the network.

## Replicates-inclusive input (DO NOT merge replicates) — 花生-style
Some users (esp. following a prior plant project like 花生/peanut) input **every biological replicate as an independent sample**, NOT time-point means. E.g. peanut `GCN 18 18` = 6 time points × 3 replicates, so 18 samples per condition. **Always ask the user whether to merge replicates or keep them** — this is a decision point, do not assume.
- For our 珊瑚菜 (14 raw samples: root1-4×3 + root6×2), replicate-inclusive means `GCN 14 14`, input = gene + 14(C1) + 14(C2) = 29 columns.
- **Cutoff behavior differs drastically from merged data.** With 14 samples, the p<0.05 PCC threshold is much LOWER (自由度增加 → 临界PCC下降): 珊瑚菜 14-sample Cutoff gave **pos≈0.74, neg≈-0.63** vs 0.91/-0.86 for 5 merged time points. Do NOT reuse the merged-data cutoff.
- **Gene-count-vs-threshold differs too.** 14-sample 珊瑚菜: 0.74→22,308 genes, 0.90→12,986, **0.925→9,464, 0.93→8,630**, 0.95→5,184. To hit 8,000-10,000 genes use ~0.925-0.93.
- root6 having only 2 replicates is statistically FINE — PCC uses all 14 samples as data points; a 2-replicate time point just contributes 2 points to the correlation, doesn't break it.
- Comparative binaries (`Cutoff/GCN/TO-GCN` main dir) read a single file with n1+n2 columns; `argv[3]`=TF file (regulators), `argv[4]`=all-gene file (targets) — same positional convention as peanut.

## Controlling network size (gene count) — PCC threshold is NOT enough
**Target size & user reminder:** the user's target is **8,000-10,000 genes** in the network. When targeting this, **ALWAYS remind the user** (do not silently pick a threshold): the chosen PCC threshold determines gene count, and raising it to hit 8k-10k drops many low-PCC (but possibly biologically relevant) edges/genes. State the tradeoff explicitly and let the user confirm the threshold before committing (they care about which genes get dropped). Example confirmation to give: "0.93 → 8,630 genes; 0.925 → 9,464; do you confirm 0.93?"
If the user wants a smaller network (e.g. 8,000-10,000 genes), **raising the PCC threshold alone is ineffective in the 0.90-0.99 range**: 0.95→23,172 genes, 0.96→23,059, 0.97→22,689, 0.98→21,620, 0.99→18,092. Coexpression networks are scale-free — hub genes stay connected even at very high PCC, so gene count drops slowly while edges collapse. **However, in the 0.995+ range the drop accelerates sharply** (empirical, 珊瑚菜 5-time-point data): 0.995→12,732 genes, 0.996→11,006, **0.9965→10,002, 0.997→8,882, 0.9975→7,705, 0.998→6,382**. So a single high threshold (0.9965-0.997) can hit the 8,000-10,000 target. Options:
- **Degree-based pruning**: keep genes with degree ≥ threshold (e.g. 度≥50→6,460 genes; pick threshold to hit target count). Loses low-degree genes that may be biologically relevant.
- **TF + top-N targets**: keep all TFs (1,869) + each TF's top-N high-PCC targets → best for TF→gene regulatory network (guarantees all TFs present, controls size via N).
- **Combination**: raise PCC to ~0.97-0.98 THEN degree-prune.
- Confirm the user's target and pruning preference before running — they care about which genes are dropped.

## Seed gene selection (official method)
TO-GCN README: "select genes that highly express in the first time point and lowly express in following time points."
Practical criteria (珊瑚菜): root1 TPM ≥ 20, strict monotonic decline root1>root2>root3>root4>root6, root6/root1 < 0.5. Yielded 139 seeds. **Ensure seeds.txt also uses cleaned short IDs.**

### CRITICAL: use ONE representative seed gene (not many seeds)
TO-GCN official examples and established plant workflows (e.g. 花生/peanut) use **a single seed gene** as the BFS root — NOT many seeds. Using many seeds (e.g. all 139) flattens the BFS level structure: 珊瑚菜 with 60 seeds-in-network gave a murky 19 levels. With a single well-chosen seed the cascade is much clearer (e.g. 20 levels with a clean early→late gradient).

**How to choose the single seed (empirical, 珊瑚菜):**
1. Compute each candidate seed's **network degree** (connectivity in the coexpression network).
2. Pick the **highest-degree seed that is also a TF with the right expression pattern** (early-high, monotonic decline). 珊瑚菜 best: `Glh1_09G003566` (AP2/ERF-ERF TF, degree=174 in 0.93 network, root1=114.7 → root6=9.7, ratio 0.084) — far above the next-best seed (degree 52).
3. A high-degree TF seed covers the largest network subgraph AND is biologically a plausible upstream regulator.
4. Optionally verify it's a TF and check its family (e.g. via iTAK `tf_classification.txt`).

```bash
# network degree of each seed (from C1+C2+ edges: col1=TF, col3=gene)
awk -F',' 'NR>1{print $1"\n"$3}' C1+C2+.csv | sort | uniq -c | sort -rn | head
# cross-reference with seeds list to find highest-degree seed in the network
```

**Decision point**: the user may want one seed, all seeds, or per-seed runs (peanut loops over seeds with `mkdir $i.c1+c2+`). ALWAYS ask which seed strategy before running BFS — this is a decision the user cares about, do not assume.

**Best single seed + expression trend plot (user requirement):** after choosing THE best single seed (highest network degree + TF + right pattern), **plot its expression trend across the time points** so the user can visually confirm it matches the seed criterion (high at first time point, declining). Save this trend plot. This is a required deliverable — the user explicitly asked to "画个表达趋势图" for the seed gene.

**Seed trend plot with 95% CI (user-approved final style):** the accepted trend plot shows the seed gene's expression mean ± **95% confidence interval** at each time point, computed from the RAW replicate data (NOT time-point means). For each time point compute `mean`, `SE = sd/sqrt(n)`, `CI = mean ± t(0.975, n-1)*SE`; plot with ggplot `geom_ribbon(alpha=0.25)` + `geom_errorbar` + points + mean value labels + line. Caveat: root6 has only 2 replicates so its CI is wider (t with 1 df) — fine, it reflects the fewer replicates. Script: `scripts/plot_seed_trend_CI.R` in the reference project. This is the final approved version — use CI version, not the plain line plot.

## iTAK TF annotation (for TF gene list)
- iTAK.py needs HMM DBs (Tfam_domain.hmm etc) + hmmscan (in base miniforge3 PATH).
- DB location: `~/miniforge3/share/itak/database/` (create + copy from a full iTAK install, e.g. `/media/nfs1/lll/softwares/iTAK/database/iTAK-db-v1/database/`).
- Run: `python iTAK.py proteins.fa -p 16 -o <outdir>` → outputs `tf_classification.txt` (3,060 TF from 40,881 proteins), `shiu_classification.txt` (protein kinases).
- Protein input: strip `.mRNA\d+` suffix from headers so IDs match expression matrix.
- TF expression file = subset of All_genes_expression.tsv for iTAK-identified TF ids.

## CRITICAL bug fix — stack smashing (REQUIRED for real plant data)
TO-GCN `*_single.cpp` (Cutoff/GCN/TO-GCN) crashes on real plant datasets with "*** stack smashing detected ***". Two root causes:
1. **Gene IDs >19 chars** overflow `char GID[20]` in `fscanf("%s", GID)` (Read_Time_Course_Data_TFs / Read_seed_genes). Fix data (clean IDs), not code. Diagnose with `g++ -g -fsanitize=address` to see exact line.
2. **PCC NaN**: zero-variance genes (all time points equal) make the PCC denominator 0 → NaN. In `histogram_calculation`, `(int)((NaN+1)*100)` → out-of-range index → stack smash. Fix r_calculator + histogram_calculation:
```cpp
// r_calculator: guard denominator
double denom = (SUM_X2-(SUM_X)*(SUM_X)/N_LD)*(SUM_Y2-(SUM_Y)*(SUM_Y)/N_LD);
if (denom <= 0) R = 0;
else { R = (SUM_XY-(SUM_X)*(SUM_Y)/N_LD)/sqrt(denom);
       if (R != R) R = 0; if (R > 1.0) R = 1.0; if (R < -1.0) R = -1.0; }
// histogram_calculation: clamp y to [0,200]
if (x != x || x > 1.0) x = 1.0; if (x < -1.0) x = -1.0;
int y = (int)((x+1)*100); if (y<0)y=0; if(y>200)y=200;
```
Apply the r_calculator fix to all three `*_single.cpp` files.

## Outputs & visualization
- `C1+.csv` (GCN_single): edges for Cytoscape
- `Node_level.csv` + `Node_relation.csv` (TO-GCN_single): levels + edges → Cytoscape (levels as node color/layout by level) AND R(igraph) layered plot (y=level, x=TF family or arbitrary). Both requested by user typically.
- **Level × time-point heatmap (validates biological meaning)**: for each TO-GCN level, compute mean & median TPM across time points → heatmap (log2 TPM). 珊瑚菜 result confirmed the core hypothesis: **low levels (L1-L3) = early-expressed (幼根高, root6/root1<1), high levels (L6-L9) = late-expressed (成熟根高, root6/root1 up to 4.7)**. This is the key sanity check that the time-ordering is biologically real.
- **Seed gene plots**: heatmap (z-score) + mean line + top-30 individual lines — shows seeds decline from root1→root6 as designed.
- **Peak-timepoint distribution per level** (stacked bar) shows the early→late progression across levels.
- **Double-ring network (user's preferred layout)**: the user's preferred way to visualize the TF-gene regulatory network is **concentric double-rings per level**, arranged left→right (L1→L13). Each level = a big ring (outer) + small ring (inner): **inner ring = TF nodes (red), outer ring = structural genes (blue)**. Edge count is large (83K) so render edges as very light grey (`edge.color=rgb(0.5,0.5,0.5,0.12)`, `edge.width=0.3`) and small nodes. Custom layout in R/igraph:
```r
# x = level center (lv*spacing), y = ring positions
for (lv in levels) {
  tf <- idx[tf_of[idx]]  # inner ring radius tf_r
  ge <- idx[!tf_of[idx]] # outer ring radius gene_r
  ang <- seq(0,2*pi,length.out=length(tf)+1)[1:length(tf)]
  xy[tf,1] <- lv*spacing + tf_r*cos(ang); xy[tf,2] <- tf_r*sin(ang)
  ang <- seq(0,2*pi,length.out=length(ge)+1)[1:length(ge)]
  xy[ge,1] <- lv*spacing + gene_r*cos(ang); xy[ge,2] <- gene_r*sin(ang)
}
```
**CRITICAL: rings come out squashed/elliptical unless you set `asp=1` in `plot()` AND match the PDF canvas physical W:H to the coordinate X:Y ratio.** Without asp=1 the igraph plot stretches coordinates so circles become ellipses. Compute `ratio=diff(xr)/diff(yr)` from the layout ranges, set `pdf(width=W, height=W/ratio)`, and pass `asp=1` to plot. Result is a wide horizontal figure (13 rings) with true circular rings, matching user's "左右稍微长一些" preference.

**CRITICAL experiential finding — fixed ring radius breaks on uneven levels:** each level has wildly different gene counts (珊瑚菜 L8=1253 structural genes vs L1=150). With a FIXED outer radius, high-count levels crowd badly (L8 ≈ 71 nodes/unit of circumference → nodes pile into an unreadable blob) while small levels look fine. **Make ring radius ADAPTIVE to node count:** `outer_r <- sqrt(level_gene_count) * scale` (scale so the largest level's outer radius ≈ 6.0). This makes big levels have bigger rings. Also scale the inner (TF) ring by `sqrt(TF_count)`. This is essential — a fixed-radius ring plot of real data will look broken.

**CRITICAL geometric tradeoff — 13 true circles in a row is inherently very wide:** after `asp=1`, arranging 13 non-overlapping rings left→right forces the figure to be extremely wide/narrow (珊瑚菜: X range ~13.5→175.5, Y range ~-6→6, i.e. aspect ~13:1). This is NOT a code bug — it's the geometric consequence of "one row of 13 true circles, non-overlapping." No parameter tweak (spacing, ring size, canvas size) removes the thinness. **This is a decision point: present the option to the user** — accept the wide single-row layout, split into 2 rows, or use a circular/spiral arrangement (Circos-style). Do NOT silently pick one; the user has strong layout preferences. When the user says "the figure looks wrong but your description is right", check ring-radius adaptivity and asp=1 first, then explain the inherent aspect-ratio tradeoff and ask their preference.

**FINAL PREFERRED LAYOUT — matplotlib TWO-CIRCLE concentric rings (one row, left→right):** after trying several options, the user's accepted layout is a **one-row left→right arrangement of 13 concentric double-ring circles** (matplotlib, `scripts/viz_two_circle_network.py`). Key properties the user explicitly approved:
- 13 levels left→right, each a concentric double-ring: **inner dashed ring = TF (orange #D95F02), outer solid ring = structural genes (green #1B9E77)**. Draw rings explicitly with `matplotlib.patches.Circle` (fill=False, inner dashed, outer solid).
- **NO gene ID labels** (`SHOW_LABEL=False`) — user asked to remove them (too many nodes, unreadable). Keep only L1-L13 labels below each circle + "TF:n Gene:n" count above.
- **Compact horizontal spacing**: LEVEL_SPACING ≈ 2.3 (outer ring radius 0.90, net gap between outer rings ≈0.5). User asked to make it progressively more compact (3.2 → 2.6 → 2.3); keep reducing spacing if they ask.
- Circle centers at `((lv-1)*SPACING, 0)`; inner TF nodes on radius R_TF=0.45, outer structural on R_STRUCTURAL=0.90 (angles evenly spread).
- `ax.set_aspect("equal")` + `ax.set_xlim/ylim` keeps rings true circles. Wide figure `FIGSIZE=(28,8)`.
- **CRITICAL performance: do NOT draw 80K+ edges with `FancyArrowPatch` in a loop** — 7+ min, never finishes. Use `matplotlib.collections.LineCollection` (segments list + one `add_collection`) → 83K edges in seconds. Node scatter in 2 calls for legend.
- Inputs: `C1+C2+_0.93.csv` (rename cols to `TF_gene_ID,type,coexpress_gene_ID,PCC_C1,PCC_C2`), `tf_gene_ids.txt`, `Node_level_*_L13.csv` (gene, level_in_GCN).
- Alternative arc/fan layout (`viz_arc_network.py`) also works (13 levels along an arc, near-square figure) if the user later prefers a more compact non-row shape; it was tried first but the two-circle one-row layout was the final accepted version.

## Level merging (collapse tail levels)

## Key biological interpretation
TO-GCN levels reflect regulatory cascade order: upstream TFs (low level) are expressed early and drive later (high-level) TFs. Validate by correlating level with expression timing — a monotonic early→late shift with increasing level confirms the method captured real developmental ordering.

Single-seed BFS (珊瑚菜, seed Glh1_09G003566, 0.93 cutoff, 14-sample) gave the cleanest gradient:
- **L1-L8 early** (root6/root1 0.26-0.83), **L9-L10 flat** (~1), **L11-L19 late** (root6/root1 up to 5.74). 8,243 genes, 1,016 TFs, 20 levels. This is the recommended presentation when the user wants a clear cascade.

## Single-seed output files (珊瑚菜 reference)
- `04.togcn/Node_level_seed_<SEED>.csv` — gene, level_in_GCN
- `05.visualize/singleseed_level_mean_heatmap.pdf/png`, `_median_heatmap` — level×time heatmap
- `05.visualize/regulatory_TF2TF_singleseed.csv` — TF→TF cascade edges with reg/tgt levels (for readable network plot)
- `05.visualize/cytoscape_singleseed_node_table.csv`, `_edge_table.csv` — Cytoscape import
- `05.visualize/singleseed_TF2TF_network.pdf/png` — TF cascade network colored by level
- R heatmap bug note: when merging levels with an expression df, the df MUST have a `gene` column (not just rownames) or `merge(by="gene")` fails with "'by' must specify a uniquely valid column".

## Level merging (collapse tail levels)
Users often want to **merge trailing sparsely-populated levels** into the last meaningful one for a cleaner cascade (e.g. L15-L20 → L14). Do it by remapping `level >= K → K` in the Node_level file, then regenerate the level×time mean/median matrices and heatmaps. 珊瑚菜 single-seed: L1-L14 after merging L15-20 into L14 (L14 = 380 genes: 223 + 157 from L15-20). Always regenerate the per-level expression data AND heatmaps from the merged level file so counts and plots stay consistent.
```python
# merge level>=15 into 14, write new Node_level_*_merged.csv
for row in reader:
    lv = int(row['level_in_GCN'])
    if lv >= 15: lv = 14
    writer.writerow([row['gene'], lv])
```

## Level adjustment criteria (user's workflow)
The user tunes the level structure based on **TWO things together**, iteratively:
1. **Gene count per level** — merge sparsely-populated tail levels so no level is trivially small (e.g. L14=380, L20=1 → merge tail up). Target a sensible distribution where each retained level has a meaningful number of genes.
2. **The expression heatmap's diagonal pattern** — the desired outcome is a **top-left → bottom-right developmental gradient**: low levels (top rows) show HIGH expression at the EARLY time point (left column), high levels (bottom rows) show HIGH expression at the LATE time point (right column). Inspect the row-scaled heatmap; if the gradient is broken (e.g. a late level's peak is in the middle, or an early level doesn't decline), **merge/adjust levels until a clean left-top→right-bottom diagonal emerges**. A diagonal pattern confirms the time-ordered cascade is biologically coherent.
So: pick the level set by balancing (1) per-level gene counts and (2) heatmap diagonal coherence — iterate merge until both look clean, then regenerate heatmaps. 珊瑚菜 final: merged to L1-L13 (after merging L14-20 into L13) because that gave the cleanest diagonal.

## Stage-wise enrichment + merged multi-group bubble plot (user requirement)
After levels are merged and validated, the user often wants to divide levels into **developmental stages** (前期/early, 中期/middle, 后期/late) and run GO/KEGG enrichment on each stage's genes, then **plot the three groups' enrichment bubbles in ONE merged figure**.

**Stage division (user-confirmed for 珊瑚菜 L1-L13):** `L1-7 early / L8-10 middle / L11-13 late`. This follows the level-expression gradient: L1-8 peak early (root6/root1<1), L9-10 flat (~1), L11-13 late (root6/root1>1). Ask the user to confirm the exact cut before running (they tune it).

**Data prep:** extract per-stage gene lists from `Node_level_*.csv` (gene, level_in_GCN) by level range; background = ALL expressed genes (the same gene_14/gene matrix used for the network, e.g. 23,232 genes — NOT just the network's level genes). All level genes must be present in background.

**Enrichment (reuse existing eggNOG annotation):** the RNA-seq project already has `*.emapper.annotations` (e.g. `RNA-seq/5.DEG/annotation/gli_shanhucai.emapper.annotations`). Build term2gene with clusterProfiler `enricher()` exactly as `RNA-seq/6.Enrichment/run_enrichment.R` does:
- **GO**: parse `GOs` column (comma-split) → TERM2GENE (columns `TERM`,`GENE` uppercase). Term2name from GO.db `AnnotationDbi::select(GO.db, keys=go_ids, columns="TERM", keytype="GOID")` (TERM column = description). **Plot ONLY GO BP** per project standard (filter results by `ONTOLOGY=="BP"` via GO.db).
- **KEGG: use `KEGG_Pathway` column (ko##### pathway level), NOT `KEGG_ko` (K-number level).** KO-level gives few (5-18) unreadable `K###` labels; Pathway-level shows readable names ("Photosynthesis", "Ribosome", "MAPK signaling pathway - plant"). Parse `KEGG_Pathway` col, keep only `^ko[0-9]{5}$` (drop `map#####` dups), apply **plant filter** (exclude prefixes `ko05/ko047/ko049/ko046/ko043` — animal pathways). Term2name from the global DB `/media/nfs1/hermes/db/kegg_pathway_names.tsv` (col1=TERM, col2=NAME; 1170 rows, reusable across projects). Pass `TERM2NAME` to enricher so Description is readable.
- Background for enrichment = the expressed-gene universe intersected with each term2gene's genes (enricher `universe=`).
- `minGSSize=5, maxGSSize=500, pvalueCutoff=0.05, qvalueCutoff=0.2`, BH adjustment.

**Merged bubble plot — THE ONLY final style (user decided 2026-08): multi-group dotplot with x = stage.** The user's final, approved enrichment figure for comparing multiple groups is the **multi-group dotplot** where **x-axis = the group (stage), y-axis = term**, dot size = Count, color = `-log10(p.adjust)` (blue→red `#4575B4`→`#D73027`). One row per unique term, each term appears as a dot in each group column it is enriched in — this lets you compare the SAME term's significance across stages at a glance. Save PDF + PNG (600 dpi), figure ~8×7.
```r
# x = Group, y = Description; top-N significant terms per group, merged & deduped
plot_df <- df %>%                                   # df = all groups' sig results + Group col
  group_by(Group) %>% arrange(p.adjust, .by_group=TRUE) %>%
  slice_head(n=10) %>% ungroup() %>% pull(Description) %>% unique()
df2 <- df %>% filter(Description %in% show_terms)
term_order <- df2 %>% group_by(Description) %>%
  summarise(min_p=min(p.adjust, na.rm=TRUE)) %>% arrange(min_p) %>% pull(Description)
df2$Description <- factor(df2$Description, levels=rev(term_order))   # top-most = most significant
ggplot(df2, aes(x=Group, y=Description)) +
  geom_point(aes(size=Count, color=logP), alpha=0.85) +
  scale_color_gradient(low="#4575B4", high="#D73027",
                       name=expression(-log[10](adjusted~italic(P)))) +
  scale_size_continuous(range=c(2.5,9), name="Count") + theme_bw(base_size=12) +
  theme(panel.grid.major=element_line(color="grey90", linewidth=0.3), panel.grid.minor=element_blank())
```
Do NOT use the facet-grid (`. ~ Group`) or per-stage single plots as the final deliverable — the user wants the x=Group multi-group dotplot for enrichment comparison. (Other styles may still be generated during exploration, but the x=Group one is the one to keep.)

**CRITICAL — Chinese labels break PDF rendering:** facet strip labels like `前期_L1-7` cause `conversion failure on '期' in 'mbcsToSbcs'` in the PDF device (PNG is fine). **Use English stage labels** `Early (L1-7) / Middle (L8-10) / Late (L11-13)` so PDF renders cleanly. Keep Chinese only in filenames/CSV if desired. Also avoid the em-dash `—` in ggplot titles (PDF substitutes `-`); use plain ASCII.

**珊瑚菜 reference result (validates biology):**
- GO BP: Early 87 / Middle 188 / Late 104 significant terms.
  - Early: chloroplast thylakoid / photosynthesis — juvenile root metabolism.
  - Middle: ribosome biogenesis, mRNA binding, ribosomal subunit — cell proliferation/protein synthesis peak.
  - Late: response to chitin, defense response to fungus, response to bacterium — mature root defense/secondary metabolism.
- KEGG Pathway: Early 12 / Middle 4 / Late 9.
  - Early: DNA replication, Photosynthesis, Proteasome.
  - Middle: Ribosome, Spliceosome, Oxidative phosphorylation.
  - Late: MAPK signaling (plant), Plant hormone signal transduction, alpha-Linolenic acid metabolism (JA precursor), Linoleic acid.
- Reference dir: `TO-GCN/07.enrich/{gene_lists,results,figures,scripts}`. Script: `run_stage_enrich.R`.

## Final HTML report (user-approved pattern, 2026-08)
After the pipeline + enrichment are done, produce a **standalone self-contained HTML report** for the TO-GCN analysis, packaged with a `report_bundle.tar.gz`. The 珊瑚菜 report (deemed "很好" by the user) is the reference template: `TO-GCN/06.report/report_14sample_5stage_togcn_shanhucai.html` + `_bundle.tar.gz`.

**Writing style: treat the reader as a PRIMARY-SCHOOL STUDENT.** Explain EVERY concept, computation step, and number in plain language with analogies. This is mandatory per the user's global report standard (`bio-html-report-standards`). Every report must:
- Explain each term on first use with a "打个比方" analogy (PCC = 步调是否一致; 共表达 = 同涨同落; 富集 = 抓球颜色; BFS = 水面涟漪扩散; TPM = 标准化活跃度; 转录因子 = 工头; 热图 = 颜色表数值; 种子基因 = 投石入水的起点).
- Walk through the whole pipeline step-by-step (数据来源 → 数据准备 → 建网 → 种子 → 分层 → 热图 → 网络图 → 富集 → 结论), with real verified numbers at each step.
- For every figure, give a "解读方法" box: what axes mean, what colors/sizes mean, what a GOOD result looks like, then the SPECIFIC observation from THIS data.

**Report structure (11 sections, user-approved):**
1. 项目背景 — research aim + data source + 打比方 overview + stat cards
2. 数据准备 — TPM matrix, gene-ID unification (strip `gene.`/`.mRNA\d+`), iTAK TF annotation (40,881 proteins → 3,060 TF → 1,869 expressed)
3. 共表达网络构建 — PCC concept, Cutoff (14-sample pos 0.74 / neg -0.63), GCN (1,869 TF × 23,233 genes), 0.93 subnetwork (8,630 genes / 83,344 edges)
4. 种子基因选择 — criteria, chosen seed (Glh1_09G003566, AP2/ERF-ERF, degree 174, root1=114.7→root6=9.7), expression-trend figure with 95% CI
5. 时间顺序分层 — BFS concept, 8,243 genes / 1,016 TF, level merging to L1–L13, per-level gene/TF table
6. 表达热图验证 — 4-panel MeanExpHeatmap, diagonal pattern (L1-7 early / L8-10 middle / L11-13 late)
7. 网络可视化 — two-circle ring network figure + reading guide
8. 三阶段富集分析 — stage split table, GO (BP) + KEGG (Pathway) enrichment methods, multigroup dotplots + biology
9. 总结与结论 — bullet summary
10. 软件与参考文献 — tool/version table + citations
11. 附图附表清单 — figure/table index matching the bundle

**Embedded figures (5, base64 ≤1500px, each <700KB):** (1) seed trend CI, (2) stitched 4-panel heatmap, (3) two-circle network, (4) GO multigroup dotplot, (5) KEGG multigroup dotplot. Resize with PIL before base64. Stitch heatmap panels horizontally with PIL.

**Bundle structure** (`report_<Nsamples>_<Nstage>_togcn_<species>_bundle.tar.gz`): `附图/` (图N_描述.png + 图N_描述.pdf, 600dpi+vector) and `附表/` (表N_描述.tsv). For TO-GCN: 表1=level_distribution.tsv, 表2=GO_BP_by_stage.tsv, 表3=KEGG_by_stage.tsv. Report filename: `report_14sample_5stage_togcn_shanhucai.html`.

**Generation rule (from `bio-html-report-standards`):** write the HTML via Python `open(path,'w',encoding='utf-8')`, NEVER via the write_file tool (breaks CSS on Windows). Verify: 0 external links, 0 `<script>`, 0 `<link stylesheet>`, 0 "Hermes/Agent", 0 `<code>`, correct species name, all numbers verified against actual files.

**珊瑚菜 report key numbers to reuse (verify against files before writing):** root1-4×3 + root6×2 = 14 samples; 23,233 genes; 40,881 proteins → 3,060 TF / 1,493 PK / 1,869 expressed TF; Cutoff 0.74/-0.63; 0.93 → 8,630 genes / 83,344 edges; seed Glh1_09G003566 (AP2/ERF-ERF, deg 174, root1=114.68→root6=9.69); 8,243 leveled genes / 1,016 TF; L1-L13 counts (L1:163, L2:292, L3:459, L4:466, L5:519, L6:797, L7:742, L8:1351, L9:508, L10:522, L11:470, L12:1125, L13:829); stages Early 3,438 / Middle 2,381 / Late 2,424; GO BP 87/188/104; KEGG Pathway 12/4/9.

## Related skills
- `wgcna-tf-network` — iTAK TF annotation + WGCNA module network (complementary, static modules vs time-ordered)
- `rnaseq-hisat2-featurecounts-deseq2` — upstream RNA-seq pipeline producing the TPM/counts matrix
- `bio-gene-regulatory-networks-coexpression-networks` — WGCNA/CEMiTool coexpression general reference

## Citations

### Tools & methods
- **TO-GCN**: Cui et al. "A Comparative Transcriptomics Method to Infer Time-ordered Gene Coexpression Networks". Software repo: https://github.com/petitmingchang/TO-GCN
- **iTAK**: Zheng Y et al. iTAK: A Program for Genome-wide Prediction and Classification of Plant Transcription Factors, Transcriptional Regulators, and Protein Kinases. Mol Plant 2016, 9(12):1667-1670. doi:10.1016/j.molp.2016.09.014
- **clusterProfiler**: Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2021, 2(3):100141. doi:10.1016/j.xinn.2021.100141
- **igraph**: Csárdi G, Nepusz T. The igraph software package for complex network research. InterJournal, Complex Systems, 2006, 1695(5):1-9
- **ggplot2**: Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. ISBN 978-3-319-24277-4

### BibTeX
```bibtex
@software{togcn,
  author = {Cui, J. and others},
  title = {TO-GCN: A Comparative Transcriptomics Method to Infer Time-ordered Gene Coexpression Networks},
  url = {https://github.com/petitmingchang/TO-GCN},
  year = {2019}
}
@article{itak2016,
  author = {Zheng, Y. and Jiao, C. and Sun, H. and Rosli, H. G. and Pombo, M. A. and Zhang, P. and Banf, M. and Dai, X. and Martin, G. B. and Giovannoni, J. J. and Zhao, P. X. and Rhee, S. Y. and Fei, Z.},
  title = {iTAK: A Program for Genome-wide Prediction and Classification of Plant Transcription Factors, Transcriptional Regulators, and Protein Kinases},
  journal = {Molecular Plant},
  year = {2016},
  volume = {9},
  number = {12},
  pages = {1667-1670},
  doi = {10.1016/j.molp.2016.09.014}
}
@article{clusterprofiler2021,
  author = {Wu, T. and Hu, E. and Xu, S. and Chen, M. and Guo, P. and Dai, Z. and Feng, T. and Zhou, L. and Tang, W. and Zhan, L. and Fu, X. and Liu, S. and Bo, X. and Yu, G.},
  title = {clusterProfiler 4.0: A universal enrichment tool for interpreting omics data},
  journal = {The Innovation},
  year = {2021},
  volume = {2},
  number = {3},
  pages = {100141},
  doi = {10.1016/j.xinn.2021.100141}
}
@article{igraph2006,
  author = {Csárdi, G. and Nepusz, T.},
  title = {The igraph software package for complex network research},
  journal = {InterJournal, Complex Systems},
  year = {2006},
  volume = {1695},
  number = {5},
  pages = {1-9}
}
```

### Methods template sentence
> Transcription factors (TFs) were annotated using iTAK (Zheng et al., 2016). A time-ordered gene coexpression network (TO-GCN) was constructed from time-series RNA-seq data following Cui et al.; genes were assigned to temporal levels by breadth-first search from a single seed gene. Level-specific gene sets were tested for GO and KEGG enrichment using clusterProfiler (Wu et al., 2021). Network analysis was performed with igraph (Csárdi & Nepusz, 2006) and visualization with ggplot2 (Wickham, 2016).