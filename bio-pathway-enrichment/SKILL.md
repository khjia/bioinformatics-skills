---
name: bio-pathway-enrichment
description: "Complete pathway enrichment analysis: GO (ORA + GSEA), KEGG, Reactome, WikiPathways enrichment with clusterProfiler. Includes visualization (dot plots, bar plots, cnetplot, emapplot), gene ID conversion, and background universe selection."
---


## When to Use ORA vs GSEA


| Scenario | Method | Why |
|----------|--------|-----|
| Clear DE gene list with arbitrary cutoff (padj + FC) | ORA, but consider GSEA instead | ORA discards magnitude; GSEA uses all genes ranked by statistic |
| Genes from co-expression module, GWAS loci, screen hits | ORA | No ranking available; ORA is appropriate |
| All genes with DE statistics available | GSEA (gseGO) | Avoids arbitrary cutoff; detects subtle coordinated changes |
| Very few DE genes (< 20) | GSEA | ORA has no power with small lists |
| RNA-seq with known length bias | GOseq (goseq package) | Standard ORA ignores length bias; longer genes are more likely DE |

ORA converts continuous measures into binary (significant/not), losing information. When in doubt, run both ORA and GSEA and compare.


## Core Pattern


**Goal:** Identify enriched Gene Ontology terms in a gene list from differential expression or similar analyses.

**Approach:** Test for over-representation of GO terms using the hypergeometric test via clusterProfiler enrichGO.

**"Run GO enrichment on my gene list"** → Test whether biological process, molecular function, or cellular component terms are over-represented among significant genes.

```r
library(clusterProfiler)
library(org.Hs.eg.db)  # Human - change for other organisms

ego <- enrichGO(
    gene = gene_list,           # Character vector of gene IDs
    OrgDb = org.Hs.eg.db,       # Organism annotation database
    keyType = 'ENTREZID',       # ID type: ENSEMBL, SYMBOL, ENTREZID, etc.
    ont = 'BP',                 # BP, MF, CC, or ALL
    pAdjustMethod = 'BH',       # p-value adjustment method
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
```


## Prepare Gene List from DE Results


**Goal:** Extract significant gene IDs from differential expression results and convert to the format required by enrichGO.

**Approach:** Filter DE results by adjusted p-value and fold change, then convert gene symbols to Entrez IDs using bitr.

```r
library(dplyr)

de_results <- read.csv('de_results.csv')

sig_genes <- de_results %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(gene_id)

# If using gene symbols, convert to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```


## ID Conversion with bitr


**Goal:** Convert between gene identifier types (Ensembl, Symbol, Entrez) for compatibility with enrichment tools.

**Approach:** Use clusterProfiler bitr to map between ID types using organism annotation databases.

```r
# Check available key types
keytypes(org.Hs.eg.db)

# Convert between ID types
converted <- bitr(genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Multiple output types
converted <- bitr(genes, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'), OrgDb = org.Hs.eg.db)
```


## Background Universe (Critical)


**Goal:** Improve enrichment specificity by restricting the background to genes actually tested in the experiment.

**Approach:** Pass all expressed genes (not just significant ones) as the universe parameter to enrichGO.

The background must be genes that *could have* appeared in the list. Getting this wrong is the single most common ORA error (95% of published analyses fail to specify an appropriate background). Using the whole genome (~20,000 genes) when only 12,000 were expressed inflates significance for tissue-specific pathways.

| Experiment Type | Correct Background |
|----------------|-------------------|
| RNA-seq | All genes with detectable expression (e.g., > 1 CPM in >= N samples) |
| Microarray | All probes on the array (mapped to genes) |
| Proteomics | All detected proteins |
| Targeted panel | Only genes on the panel |

```r
# Background = all genes that were tested (NOT the full genome)
# For DESeq2: genes with non-NA pvalue survived independent filtering
all_tested <- de_results$gene_id[!is.na(de_results$pvalue)]
universe_ids <- bitr(all_tested, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego <- enrichGO(
    gene = gene_list,
    universe = universe_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)
```

**Warning:** clusterProfiler silently drops unannotated genes from the background. To prevent this: `options(enrichment_force_universe = TRUE)` before running enrichGO.


## All Three Ontologies


```r
# Run all ontologies at once
ego_all <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'ALL',  # BP, MF, and CC combined
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)

# Results include ONTOLOGY column
head(as.data.frame(ego_all))
```


## Make Results Readable


```r
# Convert Entrez IDs to gene symbols in results
ego_readable <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# Or use readable = TRUE directly (only works with ENTREZID input)
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    readable = TRUE  # Converts to symbols
)
```


## Extract and Export Results


```r
# View top results
head(ego)

# Convert to data frame
results_df <- as.data.frame(ego)

# Key columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count

# Export to CSV
write.csv(results_df, 'go_enrichment_results.csv', row.names = FALSE)

# Filter for specific criteria
sig_terms <- results_df[results_df$p.adjust < 0.01 & results_df$Count >= 5, ]
```


## Simplify Redundant Terms


**Goal:** Remove highly similar GO terms to reduce redundancy in enrichment results.

**Approach:** Cluster GO terms by semantic similarity and retain representative terms using the simplify function.

GO terms form a DAG (directed acyclic graph), not a flat list. If "mitotic cell cycle" is enriched, parent terms ("cell cycle", "cell cycle process") will also be enriched because they contain supersets of the same genes. Always simplify before interpretation.

```r
# Remove redundant GO terms (keeps representative terms)
ego_simplified <- simplify(ego, cutoff = 0.7, by = 'p.adjust', select_fun = min)

# measure options: 'Wang' (default, graph-based, stable across releases),
# 'Resnik', 'Lin', 'Jiang', 'Rel' (IC-based, depend on annotation version)
ego_simplified <- simplify(ego, cutoff = 0.7, measure = 'Wang')
```

**Limitations:** `simplify()` does NOT work with `ont='ALL'` -- run BP, MF, CC separately. Cutoff 0.7 is a reasonable default; lower retains more terms, higher is more aggressive.


## Different Organisms


```r
# Mouse
library(org.Mm.eg.db)
ego_mouse <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = 'BP')

# Zebrafish
library(org.Dr.eg.db)
ego_zfish <- enrichGO(gene = genes, OrgDb = org.Dr.eg.db, ont = 'BP')

# Yeast
library(org.Sc.sgd.db)
ego_yeast <- enrichGO(gene = genes, OrgDb = org.Sc.sgd.db, ont = 'BP', keyType = 'ORF')
```


## Group GO Terms by Ancestor


**Goal:** Classify genes by broad GO slim categories for a high-level functional overview.

**Approach:** Use groupGO to assign genes to GO terms at a specific hierarchy level.

```r
# Classify genes by GO slim categories
ggo <- groupGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    level = 3,  # GO hierarchy level
    readable = TRUE
)
```


## Key Parameters


| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| OrgDb | required | Organism database |
| keyType | ENTREZID | Input ID type |
| ont | BP | BP, MF, CC, or ALL |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value (FDR) threshold |
| pAdjustMethod | BH | BH, bonferroni, etc. |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per term |
| maxGSSize | 500 | Max genes per term |
| readable | FALSE | Convert to symbols |


## Interpreting Results


Always examine effect size alongside p-values. A pathway with 500 genes can achieve p < 1e-15 with a modest 1.2x fold enrichment, while a 10-gene pathway with 4x enrichment at p = 0.01 is biologically more interesting.

- **Fold enrichment** = GeneRatio / BgRatio. Values > 2 suggest strong enrichment.
- **Count**: number of query genes in the term. Very large counts (> 50) may indicate overly broad terms.
- `minGSSize=10, maxGSSize=500` filters out uninformative extremes.


## Gene ID Mapping Pitfalls


- **Many-to-many mappings**: one Ensembl gene can map to multiple Entrez IDs. Deduplicate after `bitr()` to avoid counting genes multiple times.
- **Lost genes**: if > 15% of genes fail to convert, results may be unreliable. Always report the conversion rate.
- **Best practice**: use the same ID type throughout the pipeline. Convert at the last step if possible.


## RNA-seq Gene Length Bias


In RNA-seq, longer transcripts produce more fragments, increasing statistical power to detect DE. This systematically biases ORA toward pathways enriched in long genes (extracellular matrix, cell adhesion) and against short-gene pathways (ribosomal, mitochondrial). Standard normalization (RPKM, TMM) does NOT fix this.

For length-corrected GO enrichment, use GOseq:
```r
library(goseq)
pwf <- nullp(de_vector, 'hg38', 'ensGene', bias.data = gene_lengths)
goseq_results <- goseq(pwf, 'hg38', 'ensGene', method = 'Wallenius')
```


## Related Skills


- kegg-pathways - KEGG pathway enrichment
- gsea - Gene Set Enrichment Analysis for GO
- enrichment-visualization - Visualize enrichment results
- differential-expression/de-results - Generate input gene lists
- **sci-pathway-enrichment** — Python-based GSEA + multi-database ORA (gseapy, g:Profiler, MSigDB, Reactome, WikiPathways). Use when you need GSEA or databases beyond GO, or prefer a Python toolchain.


## Interpretation Pitfalls


1. **Significance ≠ relevance**: A 500-gene term can reach p < 1e-15 with modest 1.2× enrichment, while a 10-gene term with 4× enrichment at p = 0.01 is often biologically more interesting. Always check fold enrichment and gene count alongside p-values.
2. **Tiny gene sets**: GO terms with <5 genes can reach nominal significance by chance. Use `minGSSize=10` minimum.
3. **Multiple-testing across databases**: FDR is computed *within* each enrichment run. Running GO + KEGG + Reactome multiplies tests — report per-database FDR and stay conservative.
4. **Reproducibility**: GO annotations and enrichment libraries are versioned and drift over time. Record `OrgDb` version, `clusterProfiler` version, and analysis date in Methods.


## Prepare Gene List


**Goal:** Extract significant Entrez gene IDs from DE results in the format required by enrichKEGG.

**Approach:** Filter by significance thresholds and convert gene symbols to Entrez IDs (KEGG requires NCBI Entrez).

```r
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

# KEGG requires NCBI Entrez gene IDs (kegg, ncbi-geneid)
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```


## KEGG ID Conversion


**Goal:** Convert between KEGG-specific identifiers and other gene ID formats.

**Approach:** Use bitr_kegg to map between kegg, ncbi-geneid, ncbi-proteinid, and uniprot ID types.

```r
# Convert between KEGG and other IDs
kegg_ids <- bitr_kegg(gene_list, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')

# Available types: kegg, ncbi-geneid, ncbi-proteinid, uniprot
```


## Run KEGG Pathway Enrichment


**Goal:** Perform KEGG pathway over-representation analysis with customizable parameters.

**Approach:** Run enrichKEGG with specified organism, ID type, and statistical thresholds.

```r
kk <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',
    keyType = 'ncbi-geneid',    # or 'kegg'
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = 500
)

# View results
head(kk)
results <- as.data.frame(kk)
```


## KEGG Module Enrichment


**Goal:** Test for enrichment of KEGG modules (smaller functional units than pathways).

**Approach:** Use enrichMKEGG which tests against KEGG module definitions rather than full pathways.

```r
# KEGG modules are smaller functional units than pathways
mkk <- enrichMKEGG(
    gene = gene_list,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```


## Common Organism Codes


| Code | Organism | Notes |
|------|----------|-------|
| hsa | Human (Homo sapiens) | |
| mmu | Mouse (Mus musculus) | |
| rno | Rat (Rattus norvegicus) | |
| dre | Zebrafish (Danio rerio) | |
| dme | Fruit fly (Drosophila) | |
| cel | Worm (C. elegans) | |
| sce | Yeast (S. cerevisiae) | |
| ath | Arabidopsis thaliana | |
| eco | E. coli K-12 | Bacterial |
| pae | P. aeruginosa PAO1 | Bacterial |
| bsu | B. subtilis 168 | Bacterial |
| sau | S. aureus N315 | Bacterial |
| mtc | M. tuberculosis H37Rv | Bacterial |
| ko | KEGG Orthology | Cross-species, use with KO IDs |

KEGG covers 8,000+ organisms. Always verify the code for the specific strain:
```r
search_kegg_organism('Pseudomonas', by = 'scientific_name')
search_kegg_organism('aeruginosa', by = 'scientific_name')
```


## Browse KEGG Pathways


**Goal:** Visualize enriched genes overlaid on KEGG pathway diagrams.

**Approach:** Use browseKEGG for interactive browser view or pathview to generate annotated pathway images.

```r
# View pathway in browser (opens KEGG website)
browseKEGG(kk, 'hsa04110')

# Download pathway image
library(pathview)
pathview(gene.data = gene_list, pathway.id = 'hsa04110', species = 'hsa')
```


## Compare Multiple Gene Lists


**Goal:** Compare KEGG pathway enrichment across multiple gene lists (e.g., upregulated vs downregulated).

**Approach:** Use compareCluster with enrichKEGG to run enrichment per group and visualize with dotplot.

```r
# Compare KEGG enrichment across groups
gene_lists <- list(
    up = up_genes,
    down = down_genes
)

ck <- compareCluster(
    geneClusters = gene_lists,
    fun = 'enrichKEGG',
    organism = 'hsa'
)

dotplot(ck)
```


## Prokaryotic / Non-Model Organism KEGG


Bacteria and non-model organisms do NOT use org.*.eg.db packages or bitr(). Bacterial genes use locus tags (e.g., PA0001 for P. aeruginosa, b0001 for E. coli) that map directly as KEGG gene IDs.

```r
# Bacterial KEGG ORA -- no bitr() or OrgDb needed
# Gene IDs should be locus tags matching the KEGG genome
kegg_bac <- enrichKEGG(
    gene = sig_locus_tags,       # e.g., c('PA0001', 'PA0612', 'PA3476')
    organism = 'pae',            # P. aeruginosa PAO1
    keyType = 'kegg',            # use locus tags directly
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)

# Note: setReadable() requires an OrgDb which does not exist for most bacteria
# Instead, map gene IDs manually or use KEGG gene names from the result
```

For organisms without KEGG strain-specific annotation, use KEGG Orthology (KO) with organism = 'ko'. Map genes to KO IDs via eggNOG-mapper or BlastKOALA first.


## Multi-Condition Comparison


**Goal:** Find shared and condition-specific enriched pathways across experimental conditions.

**Approach:** Run enrichKEGG per condition, then use set operations on significant pathway IDs. Do NOT compare p-values across conditions (they depend on sample size and DE gene count).

```r
# Run enrichment per condition
kk_A <- enrichKEGG(gene = sig_genes_A, organism = 'hsa', pvalueCutoff = 0.05)
kk_B <- enrichKEGG(gene = sig_genes_B, organism = 'hsa', pvalueCutoff = 0.05)

# Set operations on enriched pathway IDs
paths_A <- as.data.frame(kk_A)$ID
paths_B <- as.data.frame(kk_B)$ID
shared <- intersect(paths_A, paths_B)
only_A <- setdiff(paths_A, paths_B)
only_B <- setdiff(paths_B, paths_A)

# Or use compareCluster for side-by-side visualization
gene_clusters <- list(ConditionA = sig_genes_A, ConditionB = sig_genes_B)
ck <- compareCluster(geneClusters = gene_clusters, fun = 'enrichKEGG', organism = 'hsa')
dotplot(ck, showCategory = 10)
```

For proper multi-contrast enrichment that avoids p-value comparison pitfalls, use the mitch package (rank-MANOVA approach).


## Notes


- **No readable parameter** - use `setReadable()` with OrgDb (eukaryotes only)
- **Requires internet** - queries KEGG database online
- **use_internal_data** - set TRUE to use cached KEGG data (may be outdated)
- **Pathway IDs** - format is organism code + 5 digits (e.g., hsa04110)
- **Licensing** - KEGG data is free for academic web browsing but bulk downloads and commercial use require a license; for reproducibility-critical work, consider Reactome or WikiPathways (fully open)
- **Background universe** - always specify; default uses all KEGG-annotated genes which inflates significance


## Core Concept


GSEA uses **all genes ranked by a statistic** (log2FC, signed p-value) rather than a subset of significant genes. It finds gene sets where members are enriched at the top or bottom of the ranked list.


## When to Use GSEA vs ORA


| Scenario | Preferred | Why |
|----------|-----------|-----|
| Have ranked DE results for all genes | GSEA | Uses full information; no arbitrary cutoff |
| Biological signal involves many modest but coordinated changes | GSEA | Core strength -- detects "distributed enrichment" ORA misses |
| Gene list NOT from ranking (co-expression module, GWAS hits) | ORA | No meaningful ranking exists |
| Few total measured genes, cannot construct meaningful ranking | ORA | GSEA needs large ranked lists to be powerful |

In benchmarks, GSEA-family methods outperform ORA by ~35% higher F1 score on simulated data. GSEA is strictly preferred for DE-derived analyses.


## Prepare Ranked Gene List


**Goal:** Create a sorted named vector of gene-level statistics suitable for GSEA input.

**Approach:** Extract fold changes (or other statistics) from DE results, name by gene ID, and sort in decreasing order.

**"Run GSEA on my differential expression results"** → Rank all genes by expression statistic and test whether predefined gene sets cluster toward the extremes of the ranked list.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')

# Create named vector: values = statistic, names = gene IDs
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$gene_id

# Sort in decreasing order (REQUIRED)
gene_list <- sort(gene_list, decreasing = TRUE)
```


## Convert Gene IDs for GSEA


**Goal:** Map gene symbols to Entrez IDs while preserving the ranked statistic values.

**Approach:** Use bitr for ID conversion, then rebuild the named sorted vector with Entrez IDs as names.

```r
# Convert symbols to Entrez IDs
gene_ids <- bitr(names(gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Create ranked list with Entrez IDs
gene_list_entrez <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID[match(names(gene_list_entrez), gene_ids$SYMBOL)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
```


## Ranking Metric Selection


**Goal:** Choose a ranking metric that balances magnitude and significance for GSEA.

**Approach:** The ranking metric choice matters enormously. Match the metric to the DE tool used.

| DE Tool | Recommended Metric | Column | Why |
|---------|-------------------|--------|-----|
| DESeq2 | Wald statistic | `stat` | Combines effect size + variance; best overall for RNA-seq |
| DESeq2 (shrunk) | Shrunken log2FC | `log2FoldChange` | Use `type='apeglm'` or `type='ashr'`; NOT `type='normal'` (deprecated) |
| limma/voom | Moderated t-statistic | `t` | Borrows strength across genes |
| edgeR | Signed p-value | `sign(logFC) * -log10(PValue)` | edgeR has no Wald-equivalent column |

```r
# DESeq2 Wald statistic (default recommendation)
gene_list <- de_results$stat
names(gene_list) <- de_results$gene_id
gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)

# Signed p-value (for edgeR or when Wald stat unavailable)
# Replace p=0 with small value to avoid Inf
pvals <- pmax(de_results$pvalue, 1e-300)
gene_list <- -log10(pvals) * sign(de_results$log2FoldChange)
names(gene_list) <- de_results$gene_id
gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)
```

**Never use:** shrunken log2FC from `lfcShrink(type='normal')` -- the prior distorts rankings. Also: `lfcShrink()` with type='apeglm'/'ashr' drops the `stat` column, so pull stat from unshrunk `results(dds)` if needed.


## GSEA with GO


**Goal:** Detect coordinated expression changes across GO gene sets without requiring a significance cutoff.

**Approach:** Run gseGO on a ranked gene list, testing whether GO term members are enriched at the top or bottom of the list.

```r
gse_go <- gseGO(
    geneList = gene_list_entrez,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',                     # BP, MF, CC, or ALL
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    pAdjustMethod = 'BH'
)

# Make readable
gse_go <- setReadable(gse_go, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```


## GSEA with KEGG


**Goal:** Identify KEGG pathways with coordinated expression changes across all genes.

**Approach:** Run gseKEGG on the ranked gene list using KEGG pathway definitions.

```r
gse_kegg <- gseKEGG(
    geneList = gene_list_entrez,
    organism = 'hsa',
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

# Make readable
gse_kegg <- setReadable(gse_kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```


## GSEA with Custom Gene Sets


**Goal:** Run GSEA against user-provided or non-standard gene set collections.

**Approach:** Load a GMT file and use the generic GSEA function with TERM2GENE mapping.

```r
# Read GMT file (Gene Matrix Transposed)
gene_sets <- read.gmt('msigdb_hallmarks.gmt')

gse_custom <- GSEA(
    geneList = gene_list_entrez,
    TERM2GENE = gene_sets,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05
)
```


## MSigDB Gene Sets


**Goal:** Run GSEA using curated gene set collections from the Molecular Signatures Database.

**Approach:** Retrieve gene sets via msigdbr, format as TERM2GENE data frame, and run GSEA.

```r
# Use msigdbr package for MSigDB gene sets
library(msigdbr)

# Hallmark gene sets
hallmarks <- msigdbr(species = 'Homo sapiens', category = 'H')
hallmarks_t2g <- hallmarks[, c('gs_name', 'entrez_gene')]

gse_hallmark <- GSEA(
    geneList = gene_list_entrez,
    TERM2GENE = hallmarks_t2g,
    pvalueCutoff = 0.05
)

# Other categories: C1 (positional), C2 (curated), C3 (motif), C5 (GO), C6 (oncogenic), C7 (immunologic)
```


## Understanding Results


```r
# View results
head(gse_go)
results <- as.data.frame(gse_go)

# Key columns:
# - NES: Normalized Enrichment Score (positive = upregulated, negative = downregulated)
# - pvalue: Nominal p-value
# - p.adjust: FDR-adjusted p-value
# - core_enrichment: Leading edge genes
```


## Interpreting NES (Normalized Enrichment Score)


| NES | Interpretation |
|-----|----------------|
| Positive (> 0) | Gene set enriched in upregulated genes |
| Negative (< 0) | Gene set enriched in downregulated genes |
| |NES| > 1.5 | Strong enrichment (but see caveats below) |

**Correct interpretation order:**
1. Check FDR first. Use FDR < 0.25 (Broad Institute recommendation) or FDR < 0.05 (common in publications). High |NES| with non-significant FDR is meaningless.
2. Use NES for prioritization among significant results.
3. Examine the leading edge genes to understand what drives the signal.

**NES caveats:** Very large gene sets (> 500 genes) can achieve high |NES| even randomly. Very small sets (< 10 genes) can be driven by a single outlier. Always cross-check with minGSSize/maxGSSize filtering.


## Leading Edge Interpretation


The `core_enrichment` column contains the "leading edge" genes -- those driving the enrichment signal. These appear before the enrichment peak in the ranked list.

- **High leading edge count, concentrated at the extreme of the ranked list:** Strong, trustworthy enrichment. The pathway's genes are coordinated at one end.
- **Low leading edge count:** Enrichment may be driven by 1-2 extreme outlier genes, not coordinated pathway regulation. Inspect the individual genes.
- The leading edge genes are the most biologically actionable output of GSEA -- use them for downstream analysis (pathway visualization, network analysis).


## Export Results


**Goal:** Save GSEA results and extract leading edge genes for downstream analysis.

**Approach:** Convert enrichment object to data frame, export to CSV, and parse core_enrichment for driving genes.

```r
results_df <- as.data.frame(gse_go)
write.csv(results_df, 'gsea_go_results.csv', row.names = FALSE)

# Get leading edge genes for a term
leading_edge <- strsplit(results_df$core_enrichment[1], '/')[[1]]
```


## Duplicate Gene Handling


Duplicate gene IDs in the ranked list will bias enrichment scores. After ID conversion, some genes may map to multiple IDs. Always deduplicate:

```r
# Remove duplicates -- keep the entry with the largest absolute value
gene_list <- gene_list[!duplicated(names(gene_list))]

# Or more carefully, keep the most extreme signal per gene:
gene_df <- data.frame(id = names(gene_list), val = gene_list)
gene_df <- gene_df[order(-abs(gene_df$val)), ]
gene_df <- gene_df[!duplicated(gene_df$id), ]
gene_list <- setNames(gene_df$val, gene_df$id)
gene_list <- sort(gene_list, decreasing = TRUE)
```


## When to Use Reactome


| Scenario | Reactome? | Alternative |
|----------|-----------|-------------|
| Signaling pathway detail (reaction-level) | Yes -- best choice | KEGG (pathway-level only) |
| Metabolic pathway focus | Supplement | KEGG has stronger metabolic coverage |
| Reproducibility / open license required | Yes (CC0) | WikiPathways (CC0) |
| Non-model organism (bacteria, plants) | No (7 species only) | KEGG (8,000+ species) |
| Non-human model organism (mouse, rat, fly) | Caution | Annotations are computationally inferred via orthology from human; may contain errors |

Reactome pathways are curated by PhD-level biologists and externally peer-reviewed, making them the highest-quality curated pathway database. Human is the primary species; all others are computationally inferred.


## Core Pattern - Over-Representation Analysis


**Goal:** Identify Reactome pathways over-represented in a gene list from differential expression or other analyses.

**Approach:** Test for enrichment using the hypergeometric test via ReactomePA enrichPathway against curated peer-reviewed pathways.

**"Run pathway enrichment against Reactome"** → Test whether genes in curated Reactome pathways are over-represented among significant genes.

```r
library(ReactomePA)
library(org.Hs.eg.db)

pathway_result <- enrichPathway(
    gene = entrez_ids,         # Character vector of Entrez IDs
    organism = 'human',        # human, rat, mouse, celegans, yeast, zebrafish, fly
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    readable = TRUE            # Convert to gene symbols
)

head(as.data.frame(pathway_result))
```


## GSEA on Reactome Pathways


**Goal:** Detect coordinated expression changes in Reactome pathways using all genes ranked by a statistic.

**Approach:** Create a sorted named vector from DE results and run gsePathway for rank-based enrichment.

```r
# Create ranked gene list (named vector sorted by statistic)
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$entrez_id
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_result <- gsePathway(
    geneList = gene_list,
    organism = 'human',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    verbose = FALSE
)

head(as.data.frame(gsea_result))
```


## With Background Universe


**Goal:** Restrict enrichment testing to only genes that were actually measured in the experiment.

**Approach:** Pass all tested gene IDs as the universe parameter to enrichPathway.

```r
all_genes <- de_results$entrez_id  # All tested genes

pathway_result <- enrichPathway(
    gene = entrez_ids,
    universe = all_genes,      # Background gene set
    organism = 'human',
    pvalueCutoff = 0.05,
    readable = TRUE
)
```


## Visualization


**Goal:** Create publication-quality plots of Reactome enrichment results.

**Approach:** Use enrichplot functions (dotplot, barplot, emapplot, cnetplot, gseaplot2) on enrichment result objects.

```r
library(enrichplot)

# Dot plot
dotplot(pathway_result, showCategory = 15)

# Bar plot
barplot(pathway_result, showCategory = 15)

# Enrichment map (requires pairwise_termsim first)
pathway_result <- pairwise_termsim(pathway_result)
emapplot(pathway_result)

# Gene-concept network
cnetplot(pathway_result, categorySize = 'pvalue')

# GSEA plot
gseaplot2(gsea_result, geneSetID = 1:3)
```


## View Pathway in Browser


```r
# Open pathway in Reactome browser
viewPathway('R-HSA-109582', organism = 'human')  # Uses pathway ID

# Get pathway ID from results
top_pathway_id <- pathway_result@result$ID[1]
viewPathway(top_pathway_id, organism = 'human')
```


## Compare Clusters


**Goal:** Compare Reactome pathway enrichment across multiple gene lists (e.g., upregulated vs downregulated).

**Approach:** Use compareCluster with enrichPathway to run enrichment per group and visualize side by side.

```r
# Compare pathways across multiple gene lists
gene_clusters <- list(
    upregulated = up_genes,
    downregulated = down_genes
)

compare_result <- compareCluster(
    geneClusters = gene_clusters,
    fun = 'enrichPathway',
    organism = 'human',
    pvalueCutoff = 0.05
)

dotplot(compare_result)
```


## Supported Organisms


| Organism | Name | OrgDb |
|----------|------|-------|
| Human | human | org.Hs.eg.db |
| Mouse | mouse | org.Mm.eg.db |
| Rat | rat | org.Rn.eg.db |
| Zebrafish | zebrafish | org.Dr.eg.db |
| Fly | fly | org.Dm.eg.db |
| C. elegans | celegans | org.Ce.eg.db |
| Yeast | yeast | org.Sc.sgd.db |


## Interpretation Notes


- Reactome is very granular -- some pathways contain only 2-3 genes. Use `minGSSize = 10` to filter these out.
- The deep hierarchy means parent pathways will often appear alongside child pathways. Look for the most specific (deepest) enriched pathway.
- Always specify a background universe (all tested genes) to avoid inflated significance.
- Examine fold enrichment (GeneRatio / BgRatio), not just p-values.
- For non-human species, note that annotations are orthology-inferred and may not capture species-specific pathway biology.


## When to Use WikiPathways


WikiPathways is community-curated (wiki model), not expert or peer-reviewed like KEGG/Reactome. This means:

- **Strengths**: disease-specific and drug-related pathways not found in KEGG/Reactome; fully open (CC0 license); newer pathways contributed by active researchers; 30+ species
- **Limitations**: quality varies by pathway -- some are meticulously curated by domain experts, others may be incomplete or contributed by non-specialists
- **Best use**: complement to KEGG/Reactome, not a standalone primary database. Run WikiPathways alongside KEGG or Reactome to catch pathways unique to the WikiPathways collection.

Check the "Last edited" date and contributor for specific pathways before relying on them for key conclusions.


## GSEA on WikiPathways


**Goal:** Detect coordinated expression changes in WikiPathways using a ranked gene list.

**Approach:** Sort genes by fold change and run gseWP for rank-based enrichment testing.

```r
# Create ranked gene list
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$entrez_id
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_wp <- gseWP(
    geneList = gene_list,
    organism = 'Homo sapiens',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)

head(as.data.frame(gsea_wp))
```


## Using rWikiPathways Directly


**Goal:** Query the WikiPathways database directly for pathway metadata, gene lists, and GMT files.

**Approach:** Use rWikiPathways API functions to list organisms, retrieve pathway info, and download gene set definitions.

```r
library(rWikiPathways)

# List available organisms
listOrganisms()

# Get all pathways for an organism
human_pathways <- listPathways('Homo sapiens')

# Get pathway info
pathway_info <- getPathwayInfo('WP554')  # ACE Inhibitor Pathway

# Get genes in a pathway
pathway_genes <- getXrefList('WP554', 'H')  # HGNC symbols
pathway_entrez <- getXrefList('WP554', 'L')  # Entrez IDs

# Download pathway as GMT for custom analysis
downloadPathwayArchive(organism = 'Homo sapiens', format = 'gmt')
```


## Custom GMT-Based Analysis


**Goal:** Run enrichment using a downloaded WikiPathways GMT file for offline or custom analysis.

**Approach:** Download the GMT archive via rWikiPathways, read it with read.gmt, and run enricher.

```r
# Download WikiPathways GMT
library(rWikiPathways)
downloadPathwayArchive(organism = 'Homo sapiens', format = 'gmt', destpath = '.')

# Read GMT and run enrichment
wp_gmt <- read.gmt('wikipathways-Homo_sapiens.gmt')

wp_custom <- enricher(
    gene = entrez_ids,
    TERM2GENE = wp_gmt,
    pvalueCutoff = 0.05
)
```


## Common Organisms


| Common Name | Scientific Name |
|-------------|-----------------|
| Human | Homo sapiens |
| Mouse | Mus musculus |
| Rat | Rattus norvegicus |
| Zebrafish | Danio rerio |
| Fruit fly | Drosophila melanogaster |
| C. elegans | Caenorhabditis elegans |
| Arabidopsis | Arabidopsis thaliana |
| Yeast | Saccharomyces cerevisiae |


## WikiPathways vs Other Databases


| Feature | WikiPathways | KEGG | Reactome |
|---------|--------------|------|----------|
| Curation | Community | Expert | Peer-reviewed |
| License | Open (CC0) | Commercial | Open |
| Species | 30+ | 4000+ | 7 |
| Focus | Disease, drug | Metabolic | Signaling |
| Updates | Continuous | Ongoing | Quarterly |


## Setup


**Goal:** Load required packages for visualizing enrichment analysis results.

**Approach:** Import clusterProfiler, enrichplot, and ggplot2 which provide the plotting functions for enrichment objects.

```r
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Assume ego (enrichGO result), kk (enrichKEGG result), or gse (GSEA result) exists
```


## Dot Plot


**Goal:** Summarize enrichment results showing gene ratio, count, and significance in a single figure.

**Approach:** Use enrichplot dotplot which maps gene ratio to x-axis, term to y-axis, dot size to count, and color to p-value.

Most common visualization - shows gene ratio, count, and significance.

```r
dotplot(ego, showCategory = 20)

# Customize
dotplot(ego, showCategory = 15, font.size = 10, title = 'GO Enrichment') +
    scale_color_gradient(low = 'red', high = 'blue')

# Save
pdf('go_dotplot.pdf', width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()
```


## Bar Plot


Shows enrichment count or gene ratio.

```r
barplot(ego, showCategory = 20)

# Customize
barplot(ego, showCategory = 15, x = 'GeneRatio', color = 'p.adjust')
```


## Gene-Concept Network (cnetplot)


**Goal:** Visualize which genes contribute to multiple enriched terms, revealing shared biology.

**Approach:** Build a bipartite network connecting enriched terms to their member genes, optionally colored by fold change.

Shows relationships between genes and enriched terms.

```r
# Basic cnetplot
cnetplot(ego)

# With fold change colors
cnetplot(ego, foldChange = gene_list)

# Circular layout
cnetplot(ego, circular = TRUE, colorEdge = TRUE)

# Customize node size
cnetplot(ego, node_label = 'gene', cex_label_gene = 0.8)
```


## Enrichment Map (emapplot)


**Goal:** Identify clusters of related enriched terms by visualizing shared gene overlap.

**Approach:** Compute pairwise term similarity, then plot as a network where edges connect terms sharing genes.

Shows term-term relationships based on shared genes.

```r
# Requires pairwise_termsim first
ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt)

# Customize
emapplot(ego_pt, showCategory = 30, cex_label_category = 0.6)

# Cluster by similarity
emapplot(ego_pt, group_category = TRUE, group_legend = TRUE)
```

### pairwise_termsim() Method Selection

```r
# Default: Jaccard Coefficient (works with any gene set type)
ego_pt <- pairwise_termsim(ego)

# For GO terms: Wang semantic similarity (more biologically meaningful)
ego_pt <- pairwise_termsim(ego, method = 'Wang', semData = godata('org.Hs.eg.db', ont = 'BP'))
```

| Method | Type | When to Use |
|--------|------|-------------|
| JC (Jaccard) | Gene overlap | Default; works with KEGG, Reactome, any gene set |
| Wang | Graph-based | Best for GO; captures biological relationships independent of annotation version |
| Resnik/Lin/Jiang | IC-based | GO only; depends on annotation corpus (results change between database releases) |


## Tree Plot


Hierarchical clustering of enriched terms.

```r
ego_pt <- pairwise_termsim(ego)
treeplot(ego_pt)

# Show more categories
treeplot(ego_pt, showCategory = 30)
```


## Upset Plot


Show overlapping genes between terms.

```r
upsetplot(ego)

# Limit to specific number of terms
upsetplot(ego, n = 10)
```


## GSEA-Specific Plots


### Running Score Plot (gseaplot2)

```r
# Single gene set
gseaplot2(gse, geneSetID = 1, title = gse$Description[1])

# Multiple gene sets
gseaplot2(gse, geneSetID = 1:3)

# With subplots
gseaplot2(gse, geneSetID = 1, subplots = 1:3)

# By term ID
gseaplot2(gse, geneSetID = 'GO:0006955')
```

### Ridge Plot

Distribution of fold changes in gene sets.

```r
ridgeplot(gse)

# Top n gene sets
ridgeplot(gse, showCategory = 15)

# Order by NES
ridgeplot(gse, showCategory = 20) + theme(axis.text.y = element_text(size = 8))
```

**Reading ridge plots:**
- **Shifted right (positive values):** Gene set enriched among upregulated genes
- **Shifted left (negative values):** Gene set enriched among downregulated genes
- **Bimodal distribution:** Pathway contains both strongly up- and down-regulated genes; may indicate heterogeneous pathway with opposing components
- **Narrow peak:** Enrichment driven by a small cluster of similarly ranked genes
- **Broad distribution:** Many genes with varied rankings (more diffuse, less concentrated signal)


## GO-Specific Plot (goplot)


DAG structure of GO terms.

```r
# Only for GO enrichment results
goplot(ego)

# Specific ontology
goplot(ego_bp)  # where ego_bp is enrichGO with ont='BP'
```


## Heatplot


Gene-concept heatmap.

```r
heatplot(ego, foldChange = gene_list)

# Customize
heatplot(ego, showCategory = 15, foldChange = gene_list)
```


## Compare Multiple Analyses


**Goal:** Visualize enrichment results side by side across multiple gene lists or conditions.

**Approach:** Use dotplot on compareCluster output, optionally faceting by cluster.

```r
# Compare clusters (from compareCluster)
dotplot(ck, showCategory = 10)

# Facet by cluster
dotplot(ck) + facet_grid(~Cluster)
```


## Customize ggplot2 Elements


**Goal:** Fine-tune enrichment plots with custom titles, themes, colors, and text sizes.

**Approach:** Chain ggplot2 modifiers onto enrichplot output since all functions return ggplot2 objects.

All enrichplot functions return ggplot2 objects.

```r
p <- dotplot(ego, showCategory = 20)

# Add title
p + ggtitle('GO Biological Process Enrichment')

# Change theme
p + theme_minimal()

# Adjust text
p + theme(axis.text.y = element_text(size = 10))

# Change colors
p + scale_color_viridis_c()
```


## Save Plots


**Goal:** Export enrichment plots as publication-quality PDF or PNG files.

**Approach:** Use base R pdf/png device functions or ggplot2 ggsave to write plots to files.

```r
# PDF (vector, publication quality)
pdf('enrichment_plots.pdf', width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()

# PNG (raster)
png('dotplot.png', width = 800, height = 600, res = 100)
dotplot(ego, showCategory = 20)
dev.off()

# Using ggsave
p <- dotplot(ego)
ggsave('dotplot.pdf', p, width = 10, height = 8)
```


## Visualization Summary


| Function | Best For | Input Type |
|----------|----------|------------|
| dotplot | Overview of enrichment | ORA, GSEA |
| barplot | Simple counts/ratios | ORA |
| cnetplot | Gene-term relationships | ORA |
| emapplot | Term clustering | ORA |
| treeplot | Hierarchical grouping | ORA |
| upsetplot | Term overlap | ORA |
| gseaplot2 | Running enrichment score | GSEA |
| ridgeplot | Fold change distribution | GSEA |
| goplot | GO DAG structure | GO only |
| heatplot | Gene-concept matrix | ORA |


## Choosing the Right Visualization


| Goal | Plot | Key Tip |
|------|------|---------|
| First overview of top enriched terms | dotplot | Best starting point; shows 3 dimensions (ratio, count, p-value) |
| Which genes drive multiple enriched terms | cnetplot | Limit to 5-10 terms; use `circular = TRUE` for crowded networks |
| Identify functional modules among terms | emapplot | Run `pairwise_termsim()` first; if everything connects to everything, results are redundant |
| GSEA: detailed single-pathway view | gseaplot2 | Check where genes cluster in the ranked list |
| GSEA: overview of all enriched sets | ridgeplot | Read direction (left/right shift) and shape (narrow vs broad) |
| Compare enrichment across conditions | dotplot on compareCluster | Use `facet_grid(~Cluster)` for side-by-side panels |


## Common Visualization Mistakes


- **Too many terms**: plots with > 30 terms are unreadable. Use `showCategory = 15-20`.
- **Not simplifying GO first**: showing 15 redundant GO terms (cell cycle, cell cycle process, mitotic cell cycle...) wastes space and misleads. Run `simplify()` before plotting.
- **Missing gene set size**: always show both the overlap count and the total pathway size. A 3/5 overlap (60%) is very different from 30/500 (6%).
- **Bar plots for GSEA**: bar plots show count or enrichment. For GSEA, use NES on the x-axis, not p-value. Use dotplot or ridgeplot instead.
- **Skipping pairwise_termsim()**: emapplot and treeplot will fail or produce meaningless results without it.

