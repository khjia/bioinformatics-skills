---
name: bio-differential-expression
description: "Complete differential expression analysis with DESeq2 and edgeR. Covers count matrix creation, normalization, statistical testing, result extraction, filtering, and visualization (MA plots, volcano plots, heatmaps). Use when comparing gene expression between conditions in RNA-seq experiments."
---

# Differential Expression Analysis

Complete workflow for RNA-seq differential expression using DESeq2 and edgeR, from count matrices to publication-ready figures.

## When to use

- Compare gene expression between conditions (treatment vs control, time points, tissues)
- Input: gene count matrix + sample metadata
- Output: significant DEGs, volcano plots, heatmaps, MA plots

## Tool selection

| Tool | Method | Best for |
|------|--------|----------|
| DESeq2 | Negative binomial, Wald/LRT | General purpose, recommended default |
| edgeR | Negative binomial, QL F-test | Large datasets, quick iteration |

Both are industry standards. DESeq2 is more commonly cited; edgeR is faster for exploration.


## DESeq2: Required Libraries


```r
library(DESeq2)
library(apeglm)  # For lfcShrink with type='apeglm'
```


## DESeq2: Creating DESeqDataSet


**Goal:** Construct a DESeqDataSet object from various input formats for DE analysis.

**Approach:** Wrap count data and sample metadata into the DESeq2 container, specifying the experimental design formula.

**"Load my RNA-seq counts into DESeq2"** → Create a DESeqDataSet from a count matrix, SummarizedExperiment, or tximport object with sample metadata and a design formula.

### From Count Matrix

```r
# counts: matrix with genes as rows, samples as columns
# coldata: data frame with sample metadata (rownames must match colnames of counts)
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)
```

### From SummarizedExperiment

```r
library(SummarizedExperiment)
dds <- DESeqDataSet(se, design = ~ condition)
```

### From tximport (Salmon/Kallisto)

```r
library(tximport)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
```


## DESeq2: Standard DESeq2 Workflow


**Goal:** Run the complete DESeq2 pipeline from raw counts to shrunken log fold change estimates.

**Approach:** Create dataset, pre-filter low-count genes, set reference level, run size factor estimation + dispersion estimation + Wald test, then apply LFC shrinkage.

**"Find differentially expressed genes between treated and control"** → Test for significant expression changes between conditions using negative binomial models with empirical Bayes shrinkage.

```r
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)

# Pre-filter low count genes (recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level for condition
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2 pipeline (estimateSizeFactors, estimateDispersions, nbinomWaldTest)
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Apply log fold change shrinkage (recommended for visualization/ranking)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```


## DESeq2: Design Formulas


**Goal:** Specify the experimental design to model biological and nuisance variables.

**Approach:** Build R formula objects that encode condition, batch, and interaction terms for the GLM.

```r
# Simple two-group comparison
design = ~ condition

# Controlling for batch effects
design = ~ batch + condition

# Interaction model
design = ~ genotype + treatment + genotype:treatment

# Multi-factor without interaction
design = ~ genotype + treatment
```


## DESeq2: Specifying Contrasts


**Goal:** Extract results for specific pairwise or complex comparisons from a fitted DESeq2 model.

**Approach:** Use coefficient names or contrast vectors to define which groups to compare.

```r
# See available coefficients
resultsNames(dds)

# Results by coefficient name
res <- results(dds, name = 'condition_treated_vs_control')

# Results by contrast (compare specific levels)
res <- results(dds, contrast = c('condition', 'treated', 'control'))

# Contrast with list format (for complex designs)
res <- results(dds, contrast = list('conditionB', 'conditionA'))
```


## DESeq2: Log Fold Change Shrinkage


**Goal:** Reduce noisy fold change estimates for low-count genes to improve ranking and visualization.

**Approach:** Apply empirical Bayes shrinkage (apeglm, ashr, or normal) to moderate log fold changes toward zero.

```r
# apeglm method (default, recommended)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# ashr method (alternative)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'ashr')

# normal method (original, less recommended)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'normal')
```


## DESeq2: Setting Significance Thresholds


**Goal:** Control the stringency of differential expression calls using adjusted p-value and fold change cutoffs.

**Approach:** Set alpha for multiple testing correction and optionally apply a minimum log fold change threshold.

```r
# Default: padj < 0.1
res <- results(dds)

# Custom alpha threshold
res <- results(dds, alpha = 0.05)

# With log fold change threshold
res <- results(dds, lfcThreshold = 1)  # |log2FC| > 1
```


## DESeq2: Accessing DESeq2 Results


**Goal:** Retrieve, filter, and sort DE results for downstream use.

**Approach:** Extract results as a data frame, subset by significance, and order by p-value or fold change.

```r
# Summary of results
summary(res)

# Get significant genes
sig <- subset(res, padj < 0.05)

# Order by adjusted p-value
resOrdered <- res[order(res$padj),]

# Order by log fold change
resOrdered <- res[order(abs(res$log2FoldChange), decreasing = TRUE),]

# Convert to data frame
res_df <- as.data.frame(res)
```


## DESeq2: Result Columns


| Column | Description |
|--------|-------------|
| `baseMean` | Mean of normalized counts across all samples |
| `log2FoldChange` | Log2 fold change (treatment vs control) |
| `lfcSE` | Standard error of log2 fold change |
| `stat` | Wald statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (Benjamini-Hochberg) |


## DESeq2: Normalization and Counts


**Goal:** Obtain normalized expression values suitable for visualization and cross-sample comparison.

**Approach:** Extract size-factor-normalized counts or apply variance-stabilizing / rlog transformations.

```r
# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Get size factors
sizeFactors(dds)

# Variance stabilizing transformation (for visualization)
vsd <- vst(dds, blind = FALSE)

# Regularized log transformation (alternative, slower)
rld <- rlog(dds, blind = FALSE)
```


## DESeq2: Multi-Factor Designs


**Goal:** Account for batch or other nuisance variables while testing the effect of interest.

**Approach:** Include batch as a covariate in the design formula so DESeq2 adjusts for it during testing.

```r
# Design with batch correction
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ batch + condition)
dds <- DESeq(dds)

# Extract condition effect (controlling for batch)
res <- results(dds, name = 'condition_treated_vs_control')
```


## DESeq2: Interaction Models


**Goal:** Identify genes whose response to treatment differs between genotypes (or other factor combinations).

**Approach:** Fit a model with interaction terms and test the interaction coefficient for significance.

```r
# Interaction between genotype and treatment
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ genotype + treatment + genotype:treatment)
dds <- DESeq(dds)

# Test interaction term
res_interaction <- results(dds, name = 'genotypeKO.treatmentdrug')

# Or use contrast for difference of differences
res_interaction <- results(dds, contrast = list(
    c('genotypeKO.treatmentdrug'),
    c()
))
```


## DESeq2: Likelihood Ratio Test


**Goal:** Test whether a factor (e.g., condition) explains significant variance compared to a reduced model.

**Approach:** Compare full and reduced GLMs using a likelihood ratio test instead of Wald tests.

```r
# Compare full vs reduced model
dds <- DESeq(dds, test = 'LRT', reduced = ~ batch)

# Results from LRT
res <- results(dds)
```


## DESeq2: Pre-Filtering Strategies


**Goal:** Remove uninformative genes to reduce multiple testing burden and improve statistical power.

**Approach:** Apply count-based filters requiring minimum expression across a threshold number of samples.

```r
# Remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Keep genes with at least n counts in at least k samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Filter by expression level
keep <- rowMeans(counts(dds, normalized = TRUE)) >= 10
dds <- dds[keep,]
```


## DESeq2: Working with Existing Objects


```r
# Update design formula
design(dds) <- ~ batch + condition
dds <- DESeq(dds)

# Subset samples
dds_subset <- dds[, dds$group == 'A']

# Subset genes
dds_genes <- dds[rownames(dds) %in% gene_list,]
```


## DESeq2: Exporting Results


**Goal:** Save DE results and normalized counts to files for sharing or downstream tools.

**Approach:** Convert results to data frames and write as CSV files.

```r
# Write to CSV
write.csv(as.data.frame(resOrdered), file = 'deseq2_results.csv')

# Write normalized counts
write.csv(as.data.frame(normalized_counts), file = 'normalized_counts.csv')
```


## DESeq2: Common Errors


| Error | Cause | Solution |
|-------|-------|----------|
| "design matrix not full rank" | Confounded variables or missing levels | Check coldata for confounding |
| "counts matrix should be integers" | Non-integer counts (e.g., from tximport) | Use DESeqDataSetFromTximport() |
| "all samples have 0 counts" | Gene filtering issue | Check count matrix format |
| "factor levels not in colData" | Typo in design formula | Verify column names in coldata |


## DESeq2: Deprecated Features


| Feature | Status | Alternative |
|---------|--------|-------------|
| No-replicate designs | Removed (v1.22) | Require biological replicates |
| `betaPrior = TRUE` | Deprecated | Use `lfcShrink()` instead |
| `rlog()` for large datasets | Not recommended | Use `vst()` for >100 samples |


## DESeq2: Quick Reference: Workflow Steps


```r
# 1. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)

# 2. Pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 3. Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# 4. Run DESeq2
dds <- DESeq(dds)

# 5. Get results with shrinkage
res <- lfcShrink(dds, coef = resultsNames(dds)[2], type = 'apeglm')

# 6. Filter significant results
sig <- subset(res, padj < 0.05)
```


## DESeq2: Decision Guidance


### Shrinkage Method Selection

| Method | Use When | Limitation |
|--------|----------|------------|
| apeglm (default) | Coefficient-based comparisons (`coef=`) | Cannot use `contrast=` |
| ashr | Arbitrary contrasts needed; many coefficients | Slightly less aggressive shrinkage |
| normal | **Avoid** — over-shrinks large precise effects | Kept for backward compatibility only |

Shrinkage changes LFC estimates only, NOT p-values. Use shrunken LFCs for ranking (GSEA input, heatmap ordering) and visualization (volcano x-axis). Use un-shrunken p-values for significance calls.

### LRT vs Wald Test

| Scenario | Test |
|----------|------|
| Pairwise comparison (A vs B) | Wald (default) |
| Factor with >= 3 levels (any gene changing across conditions) | LRT with `reduced = ~ 1` |
| Time series (any temporal change) | LRT |
| Testing a specific coefficient direction | Wald |

LRT is omnibus (ANOVA-like). The LFC in LRT output is last-level-vs-reference, NOT the omnibus effect. Filter LRT results on padj only, not LFC.

### Why padj = NA

| Cause | baseMean | pvalue | padj |
|-------|----------|--------|------|
| Zero counts across all samples | 0 | NA | NA |
| Cook's distance outlier (automatic when any group >= 7 samples) | > 0 | NA | NA |
| Below independent filtering threshold | > 0 | numeric | NA |

Independent filtering optimizes a mean-expression cutoff at the `results()` step to maximize BH-adjusted rejections. This is separate from manual pre-filtering and uses the fitted model's information.

### Size Factor Alternatives

Default median-of-ratios assumes most genes are NOT differentially expressed.

| Scenario | Solution |
|----------|----------|
| High zero-inflation (single-cell) | `type = 'poscounts'` |
| Very small libraries | `type = 'iterate'` |
| Known stable reference genes available | `controlGenes` parameter |
| Prokaryotic stress (majority DE) | Spike-in normalization or `controlGenes` |

### Pre-filtering

```r
# Minimal (speed only; independent filtering handles statistical optimization)
keep <- rowSums(counts(dds)) >= 10

# Group-aware (recommended): require counts in at least the smallest group
keep <- rowSums(counts(dds) >= 10) >= min(table(dds$condition))
```

### Prokaryotic RNA-seq

Bacterial/archaeal experiments differ from eukaryotic in ways that affect DESeq2:
- Use non-spliced aligners (BWA-MEM, Bowtie2) — no introns in prokaryotes
- Polycistronic operons cause read-through between adjacent genes
- rRNA depletion is essential (80-95% rRNA without poly-A selection)
- Under stress conditions, a majority of genes may be DE, violating normalization assumptions — use spike-in normalization or `controlGenes` with known stable housekeeping genes
- KEGG organism codes are strain-specific (e.g., `pae` for P. aeruginosa PAO1); find codes with `clusterProfiler::search_kegg_organism()`
- Annotation: use Prokka/Bakta GFF files rather than Ensembl/biomaRt

### Choosing DESeq2 vs edgeR

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Standard bulk RNA-seq (n >= 3/group) | Either — expect ~90% overlap | Both perform well; concordance is high |
| Very small sample size (n = 2-3/group) | edgeR QL F-test | QL framework provides tighter FPR control with few replicates |
| Salmon/Kallisto quantification | DESeq2 via tximport | DESeqDataSetFromTximport handles offset matrix natively |
| Need LFC shrinkage for ranking/GSEA | DESeq2 | apeglm/ashr shrinkage built in; edgeR has no equivalent |
| Formal fold-change threshold testing | edgeR `glmTreat` | More flexible than DESeq2 `lfcThreshold` |
| Large datasets (>100 samples) | edgeR | Faster with C++ backend in v4 |
| Python-only environment | DESeq2 (via PyDESeq2) | No edgeR Python equivalent |

If overlap between DESeq2 and edgeR is < 60%, investigate filtering, normalization, or dispersion — discordance usually indicates a modeling issue.

### Python Alternative: PyDESeq2

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_df, metadata=metadata, design='~condition')
dds.deseq2()
stat_res = DeseqStats(dds, contrast=('condition', 'treated', 'control'))
stat_res.summary()
results_df = stat_res.results_df
```

PyDESeq2 (scverse, v0.5+) supports Wald test, multi-factor designs, apeglm shrinkage. No LRT yet. Results closely match R DESeq2.


## DESeq2: Related Skills


- edger-basics - Alternative DE analysis with edgeR
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes


## edgeR: Required Libraries


```r
library(edgeR)
library(limma)  # For design matrices and voom
```


## edgeR: Creating DGEList Object


**Goal:** Construct an edgeR container from a count matrix with sample group information.

**Approach:** Wrap raw counts and group labels into a DGEList object for normalization and testing.

**"Load my RNA-seq counts into edgeR"** → Create a DGEList from a count matrix with sample group assignments and optional gene annotations.

```r
# From count matrix
# counts: matrix with genes as rows, samples as columns
# group: factor indicating sample groups
y <- DGEList(counts = counts, group = group)

# With gene annotation
y <- DGEList(counts = counts, group = group, genes = gene_info)

# Check structure
y
```


## edgeR: Standard edgeR Workflow (Quasi-Likelihood)


**Goal:** Run the complete edgeR QL pipeline from raw counts to differentially expressed gene lists.

**Approach:** Filter, normalize (TMM), estimate dispersions, fit quasi-likelihood GLM, and test coefficients with the QL F-test.

**"Find differentially expressed genes between my groups"** → Test for significant expression differences using negative binomial models with quasi-likelihood F-tests.

```r
# Create DGEList
y <- DGEList(counts = counts, group = group)

# Filter low-expression genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM by default)
y <- calcNormFactors(y)

# Create design matrix
design <- model.matrix(~ group)

# Estimate dispersion (optional in edgeR v4+ but improves BCV plots)
y <- estimateDisp(y, design)

# Fit quasi-likelihood model
fit <- glmQLFit(y, design)

# Perform quasi-likelihood F-test
qlf <- glmQLFTest(fit, coef = 2)

# View top genes
topTags(qlf)
```


## edgeR: Filtering Low-Expression Genes


**Goal:** Remove genes with insufficient expression to reduce noise and multiple testing burden.

**Approach:** Apply automatic or manual CPM/count thresholds requiring expression in a minimum number of samples.

```r
# Automatic filtering (recommended)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Manual filtering: CPM threshold
keep <- rowSums(cpm(y) > 1) >= 2  # At least 2 samples with CPM > 1
y <- y[keep, , keep.lib.sizes = FALSE]

# Filter by minimum counts
keep <- rowSums(y$counts >= 10) >= 3  # At least 3 samples with 10+ counts
y <- y[keep, , keep.lib.sizes = FALSE]
```


## edgeR: Normalization Methods


**Goal:** Correct for differences in library composition between samples.

**Approach:** Compute TMM (or alternative) normalization factors that adjust effective library sizes.

```r
# TMM normalization (default, recommended)
y <- calcNormFactors(y, method = 'TMM')

# Alternative methods
y <- calcNormFactors(y, method = 'RLE')      # Relative Log Expression
y <- calcNormFactors(y, method = 'upperquartile')
y <- calcNormFactors(y, method = 'none')     # No normalization

# View normalization factors
y$samples$norm.factors
```


## edgeR: Design Matrices


**Goal:** Define the linear model structure for the experimental design.

**Approach:** Build model matrices encoding group, batch, and interaction terms for the GLM.

```r
# Simple two-group comparison
design <- model.matrix(~ group)

# With batch correction
design <- model.matrix(~ batch + group)

# Interaction model
design <- model.matrix(~ genotype + treatment + genotype:treatment)

# No intercept (for direct group comparisons)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
```


## edgeR: Dispersion Estimation


**Goal:** Estimate biological variability (dispersion) to parameterize the negative binomial model.

**Approach:** Compute common, trended, and gene-wise dispersions using empirical Bayes moderation.

```r
# Estimate all dispersions
y <- estimateDisp(y, design)

# Or estimate separately
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# View dispersions
y$common.dispersion
y$trended.dispersion
y$tagwise.dispersion

# Plot BCV (biological coefficient of variation)
plotBCV(y)
```


## edgeR: Quasi-Likelihood Testing


**Goal:** Test for differential expression using the quasi-likelihood framework for robust inference.

**Approach:** Fit a QL GLM and test individual coefficients, contrasts, or multiple coefficients simultaneously.

```r
# Fit QL model
fit <- glmQLFit(y, design)

# Test specific coefficient
qlf <- glmQLFTest(fit, coef = 2)

# Test with contrast
contrast <- makeContrasts(groupB - groupA, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Test multiple coefficients (ANOVA-like)
qlf <- glmQLFTest(fit, coef = 2:3)
```


## edgeR: Making Contrasts


**Goal:** Define specific pairwise or complex group comparisons for testing.

**Approach:** Use makeContrasts with a no-intercept design to specify arbitrary between-group differences.

```r
# Design without intercept
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Pairwise comparisons
contrast <- makeContrasts(
    TreatedVsControl = treated - control,
    DrugAVsControl = drugA - control,
    DrugBVsControl = drugB - control,
    DrugAVsDrugB = drugA - drugB,
    levels = design
)

# Test each contrast
qlf_treated <- glmQLFTest(fit, contrast = contrast[, 'TreatedVsControl'])
qlf_drugA <- glmQLFTest(fit, contrast = contrast[, 'DrugAVsControl'])
```


## edgeR: Accessing Results


**Goal:** Retrieve and filter DE results from the fitted model.

**Approach:** Use topTags to extract ranked gene lists with FDR-corrected p-values.

```r
# Top differentially expressed genes
topTags(qlf, n = 20)

# All results as data frame
results <- topTags(qlf, n = Inf)$table

# Summary of DE genes at different thresholds
summary(decideTests(qlf))

# Get DE genes with specific cutoffs
de_genes <- topTags(qlf, n = Inf, p.value = 0.05)$table
```


## edgeR: Result Columns


| Column | Description |
|--------|-------------|
| `logFC` | Log2 fold change |
| `logCPM` | Average log2 counts per million |
| `F` | Quasi-likelihood F-statistic |
| `PValue` | Raw p-value |
| `FDR` | False discovery rate (adjusted p-value) |


## edgeR: Alternative: Exact Test (Classic edgeR)


**Goal:** Perform a simple two-group DE test without a design matrix.

**Approach:** Use the classic edgeR exact test based on the negative binomial distribution.

```r
# For simple two-group comparison only
y <- DGEList(counts = counts, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y)

# Exact test
et <- exactTest(y)
topTags(et)
```


## edgeR: Alternative: glmLRT (Likelihood Ratio Test)


**Goal:** Test for DE using likelihood ratio tests as an alternative to the QL F-test.

**Approach:** Fit a standard GLM and compare nested models via the likelihood ratio statistic.

```r
# Fit GLM
fit <- glmFit(y, design)

# Likelihood ratio test
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
```


## edgeR: Treat Test (Log Fold Change Threshold)


**Goal:** Test whether genes exceed a minimum fold change threshold, not just differ from zero.

**Approach:** Use glmTreat to apply a fold change threshold directly in the statistical test.

```r
# Test for |logFC| > threshold
tr <- glmTreat(fit, coef = 2, lfc = log2(1.5))  # |FC| > 1.5
topTags(tr)
```


## edgeR: Multi-Factor Designs


**Goal:** Test for condition effects while adjusting for batch or other covariates.

**Approach:** Include nuisance variables in the design matrix so the QL test controls for them.

```r
# Design with batch correction
design <- model.matrix(~ batch + condition, data = sample_info)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Test condition effect (controlling for batch)
# Condition coefficient is typically the last
qlf <- glmQLFTest(fit, coef = ncol(design))
```


## edgeR: Getting Normalized Counts


**Goal:** Obtain normalized expression values for visualization and downstream analysis.

**Approach:** Compute CPM or log-CPM values using TMM-adjusted library sizes.

```r
# Counts per million (CPM)
cpm_values <- cpm(y)

# Log2 CPM
log_cpm <- cpm(y, log = TRUE)

# RPM (reads per million, same as CPM)
rpm_values <- cpm(y)

# With prior count for log transformation
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
```


## edgeR: Exporting Results


**Goal:** Save DE results and normalized counts to files for sharing or downstream tools.

**Approach:** Extract all results via topTags and write as CSV alongside CPM values.

```r
# Get all results
all_results <- topTags(qlf, n = Inf)$table

# Add gene IDs as column
all_results$gene_id <- rownames(all_results)

# Write to file
write.csv(all_results, file = 'edger_results.csv', row.names = FALSE)

# Export normalized counts
write.csv(cpm(y), file = 'cpm_values.csv')
```


## edgeR: Common Errors


| Error | Cause | Solution |
|-------|-------|----------|
| "design matrix not full rank" | Confounded variables | Check sample metadata |
| "No residual df" | Too few samples | Need more replicates |
| "NA/NaN/Inf" | Zero counts in all samples | Filter more stringently |


## edgeR: Deprecated/Changed Functions


| Old | Status | New |
|-----|--------|-----|
| `decidetestsDGE()` | Removed (v4.4) | `decideTests()` |
| `glmFit()` + `glmLRT()` | Still works | Prefer `glmQLFit()` + `glmQLFTest()` |
| `estimateDisp()` | Optional (v4+) | `glmQLFit()` estimates internally |
| `mglmLS()`, `mglmSimple()` | Retired | `mglmLevenberg()` or `mglmOneWay()` |

**Note:** `calcNormFactors()` and `normLibSizes()` are synonyms - both work.


## edgeR: Quick Reference: Workflow Steps


```r
# 1. Create DGEList
y <- DGEList(counts = counts, group = group)

# 2. Filter low-expression genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# 3. Normalize
y <- calcNormFactors(y)

# 4. Create design matrix
design <- model.matrix(~ group)

# 5. Estimate dispersion (optional in v4+)
y <- estimateDisp(y, design)

# 6. Fit quasi-likelihood model
fit <- glmQLFit(y, design)

# 7. Test for DE
qlf <- glmQLFTest(fit, coef = 2)

# 8. Get results
topTags(qlf, n = 20)
```


## edgeR: Choosing edgeR vs DESeq2


| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Standard bulk RNA-seq (n >= 3/group) | Either — expect ~90% overlap | Both perform well; concordance is high |
| Very small sample size (n = 2-3/group) | edgeR QL F-test | QL framework provides tighter FPR control with few replicates |
| Large datasets (>100 samples, many conditions) | edgeR | Faster; C++ backend in v4 |
| Salmon/Kallisto quantification | DESeq2 via tximport | DESeqDataSetFromTximport handles offset matrix natively |
| Need LFC shrinkage for ranking/GSEA | DESeq2 | apeglm/ashr shrinkage built in; edgeR has no equivalent |
| Formal fold-change threshold testing | edgeR `glmTreat` | More flexible than DESeq2 `lfcThreshold` |
| Python-only environment | DESeq2 (via PyDESeq2) | No edgeR Python equivalent |
| Omnibus test (>= 3 groups) | Either (LRT) | Both support likelihood ratio tests |

If overlap between DESeq2 and edgeR is < 60%, investigate differences in filtering, normalization, or dispersion estimation — discordance usually indicates a modeling issue, not a tool difference.


## edgeR: Decision Guidance


### Test Selection

| Test | When to Use | Key Property |
|------|-------------|--------------|
| QL F-test (`glmQLFit` + `glmQLFTest`) | **Default for most experiments** | Best FPR control with small n; accounts for dispersion uncertainty |
| LRT (`glmFit` + `glmLRT`) | Large samples (n >= 6 per group); complex designs where QL is too conservative | More powerful but can be anti-conservative with few replicates |
| Exact test (`exactTest`) | Simple two-group comparison only | No design matrix; cannot adjust for covariates |
| TREAT (`glmTreat`) | Testing against a minimum fold-change threshold | Tests H0: \|LFC\| <= threshold, not H0: LFC = 0 |

QL F-test p-values are always >= LRT p-values. In null comparisons (replicates vs replicates), QL consistently returns ~0 false positives while LRT can return many.

### edgeR v4 Changes

- Constant NB dispersion estimated from the most highly expressed genes (v3 used trended dispersions)
- `estimateDisp()` is now optional before `glmQLFit()` (but still needed for BCV plots)
- Fractional count support for transcript quantification uncertainty (Gibbs sampling)
- C++ backend for model fitting (faster)
- `decidetestsDGE()` removed; use `decideTests()` instead
- `mglmLS()`, `mglmSimple()` retired; use `mglmLevenberg()` or `mglmOneWay()`

### filterByExpr Internals

Default parameters: `min.count = 10`, `min.total.count = 15`, `large.n = 10`, `min.prop = 0.7`.

Algorithm:
1. Convert `min.count` to CPM cutoff: `min.count / median(lib.size) * 1e6`
2. Determine minimum number of samples `n` from design (smallest group size)
3. If group size > `large.n`: `n = large.n + (group_size - large.n) * min.prop`
4. Keep gene if CPM >= cutoff in >= n samples AND total count >= `min.total.count`

Example: median library 51M reads -> CPM cutoff ~0.2; 3 replicates per group -> need CPM >= 0.2 in >= 3 samples.

### Normalization Caveats

TMM/RLE both assume most genes are NOT DE. This assumption breaks in:
- Prokaryotic stress responses (majority of genes may be DE)
- Extreme perturbations or cross-species comparisons
- Single-cell with many zeros

When violated, consider spike-in normalization or supplying known stable reference genes. For prokaryotic RNA-seq, also note: use non-spliced aligners (BWA-MEM, Bowtie2), rRNA depletion is essential, and KEGG organism codes are strain-specific.


## edgeR: Related Skills


- deseq2-basics - Alternative DE analysis with DESeq2
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes


## Required Libraries


```r
library(DESeq2)  # or library(edgeR)
library(dplyr)   # For data manipulation
```


## Extracting DESeq2 Results


**Goal:** Retrieve DE statistics from a fitted DESeq2 model as a usable data frame.

**Approach:** Call results() with optional shrinkage, then convert to a data frame with gene identifiers.

```r
# Basic results
res <- results(dds)

# With specific alpha (adjusted p-value threshold)
res <- results(dds, alpha = 0.05)

# With log fold change shrinkage
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
```


## Extracting edgeR Results


**Goal:** Retrieve DE statistics from a fitted edgeR model as a data frame.

**Approach:** Use topTags with n=Inf to extract all gene-level results.

```r
# Get all results
results <- topTags(qlf, n = Inf)$table

# Add gene column
results$gene <- rownames(results)
```


## Filtering Significant Genes


**Goal:** Identify genes meeting statistical significance and biological effect size criteria.

**Approach:** Subset results by adjusted p-value, fold change magnitude, and expression level thresholds.

**"Get the significant differentially expressed genes"** → Filter DE results by adjusted p-value and fold change cutoffs to produce up- and down-regulated gene lists.

### By Adjusted P-value

```r
# DESeq2
sig_genes <- subset(res, padj < 0.05)

# edgeR
sig_genes <- subset(results, FDR < 0.05)

# Using dplyr
sig_genes <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(padj)
```

### By Fold Change

```r
# Absolute log2 fold change > 1 (2-fold change)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Up-regulated only
up_genes <- subset(res, padj < 0.05 & log2FoldChange > 1)

# Down-regulated only
down_genes <- subset(res, padj < 0.05 & log2FoldChange < -1)
```

### Combined Filters

```r
# Stringent filtering
sig_genes <- res_df %>%
    filter(padj < 0.01,
           abs(log2FoldChange) > 1,
           baseMean > 10) %>%
    arrange(padj)
```


## Ordering Results


**Goal:** Rank DE genes by statistical significance or biological effect size.

**Approach:** Sort results by adjusted p-value, absolute fold change, or mean expression.

```r
# By adjusted p-value (most significant first)
res_ordered <- res[order(res$padj), ]

# By absolute fold change (largest changes first)
res_ordered <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

# By base mean expression
res_ordered <- res[order(res$baseMean, decreasing = TRUE), ]

# Combined: significant genes ordered by fold change
sig_ordered <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
```


## Summary Statistics


**Goal:** Quantify the number of up- and down-regulated genes at chosen thresholds.

**Approach:** Count genes passing significance filters and report directional breakdown.

```r
# DESeq2 summary
summary(res)

# Manual counts
n_tested <- sum(!is.na(res$padj))
n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf('Tested: %d genes\n', n_tested))
cat(sprintf('Significant (padj < 0.05): %d genes\n', n_sig))
cat(sprintf('Up-regulated: %d genes\n', n_up))
cat(sprintf('Down-regulated: %d genes\n', n_down))

# edgeR summary
summary(decideTests(qlf))
```


## Adding Gene Annotations


**Goal:** Enrich DE results with gene symbols, descriptions, and cross-database identifiers.

**Approach:** Map Ensembl or Entrez IDs to human-readable annotations using org.db, biomaRt, or custom files.

**"Add gene names to my DE results"** → Map gene identifiers to symbols and descriptions using annotation databases, then merge with the results table.

### From Bioconductor Annotation Package

```r
library(org.Hs.eg.db)  # Human; use org.Mm.eg.db for mouse

# If gene IDs are Ensembl
res_df$symbol <- mapIds(org.Hs.eg.db,
                         keys = rownames(res_df),
                         column = 'SYMBOL',
                         keytype = 'ENSEMBL',
                         multiVals = 'first')

res_df$entrez <- mapIds(org.Hs.eg.db,
                         keys = rownames(res_df),
                         column = 'ENTREZID',
                         keytype = 'ENSEMBL',
                         multiVals = 'first')

res_df$description <- mapIds(org.Hs.eg.db,
                              keys = rownames(res_df),
                              column = 'GENENAME',
                              keytype = 'ENSEMBL',
                              multiVals = 'first')
```

### From BioMart

```r
library(biomaRt)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

annotations <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
    filters = 'ensembl_gene_id',
    values = rownames(res_df),
    mart = mart
)

# Merge with results
res_annotated <- merge(res_df, annotations,
                        by.x = 'row.names', by.y = 'ensembl_gene_id',
                        all.x = TRUE)
```

### From Custom File

```r
# Load annotation file
gene_info <- read.csv('gene_annotations.csv')

# Merge with results
res_annotated <- merge(res_df, gene_info, by = 'gene', all.x = TRUE)
```


## Exporting Results


**Goal:** Save DE results in formats suitable for sharing, publication, or downstream tools.

**Approach:** Write filtered and annotated results to CSV, Excel workbooks, or ranked gene lists for pathway analysis.

### To CSV

```r
# All results
write.csv(res_df, file = 'deseq2_all_results.csv', row.names = FALSE)

# Significant only
sig_genes <- res_df %>% filter(padj < 0.05)
write.csv(sig_genes, file = 'deseq2_significant.csv', row.names = FALSE)
```

### To Excel

```r
library(openxlsx)

# Create workbook with multiple sheets
wb <- createWorkbook()

addWorksheet(wb, 'All Results')
writeData(wb, 'All Results', res_df)

addWorksheet(wb, 'Significant')
writeData(wb, 'Significant', sig_genes)

addWorksheet(wb, 'Up-regulated')
writeData(wb, 'Up-regulated', up_genes)

addWorksheet(wb, 'Down-regulated')
writeData(wb, 'Down-regulated', down_genes)

saveWorkbook(wb, 'de_results.xlsx', overwrite = TRUE)
```

### Gene Lists for Pathway Analysis

```r
# Just gene IDs for GO/KEGG analysis
sig_gene_list <- rownames(subset(res, padj < 0.05))
write.table(sig_gene_list, file = 'significant_genes.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# With fold changes for GSEA
gsea_input <- res_df %>%
    filter(!is.na(log2FoldChange)) %>%
    select(gene, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
write.table(gsea_input, file = 'gsea_input.rnk',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
```


## Comparing Results Between Methods


**Goal:** Assess concordance between DESeq2 and edgeR results to identify robust DE genes.

**Approach:** Compute set overlaps and visualize with a Venn diagram.

```r
# Get significant genes from both methods
deseq2_sig <- rownames(subset(deseq2_res, padj < 0.05))
edger_sig <- rownames(subset(edger_results, FDR < 0.05))

# Overlap
common <- intersect(deseq2_sig, edger_sig)
deseq2_only <- setdiff(deseq2_sig, edger_sig)
edger_only <- setdiff(edger_sig, deseq2_sig)

cat(sprintf('DESeq2 significant: %d\n', length(deseq2_sig)))
cat(sprintf('edgeR significant: %d\n', length(edger_sig)))
cat(sprintf('Common: %d\n', length(common)))
cat(sprintf('DESeq2 only: %d\n', length(deseq2_only)))
cat(sprintf('edgeR only: %d\n', length(edger_only)))

# Venn diagram
library(VennDiagram)
venn.diagram(
    x = list(DESeq2 = deseq2_sig, edgeR = edger_sig),
    filename = 'de_overlap.png',
    fill = c('steelblue', 'coral')
)
```


## Multiple Testing Correction


**Goal:** Apply or compare multiple testing correction methods for DE p-values.

**Approach:** Use Benjamini-Hochberg (default), Bonferroni, or IHW for adjusted p-values.

```r
# DESeq2 uses Benjamini-Hochberg by default
# To use different methods:

# Independent Hypothesis Weighting (more powerful)
library(IHW)
res_ihw <- results(dds, filterFun = ihw)

# Manual p-value adjustment
res_df$padj_bonferroni <- p.adjust(res_df$pvalue, method = 'bonferroni')
res_df$padj_bh <- p.adjust(res_df$pvalue, method = 'BH')
res_df$padj_fdr <- p.adjust(res_df$pvalue, method = 'fdr')
```


## Handling NA Values


**Goal:** Understand and handle missing values in DE results caused by filtering or outlier detection.

**Approach:** Identify the source of NAs (zero counts, independent filtering, outliers) and remove or investigate them.

```r
# Count NAs
sum(is.na(res$padj))

# Remove genes with NA padj
res_complete <- res[!is.na(res$padj), ]

# Understand why NAs occur
# - baseMean = 0: No counts
# - NA only in padj: Outlier or low count filtered by independent filtering

# Check outliers
res[which(is.na(res$pvalue) & res$baseMean > 0), ]
```


## Quick Reference: Result Columns


### DESeq2

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized counts |
| `log2FoldChange` | Log2 fold change |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (BH) |

### edgeR

| Column | Description |
|--------|-------------|
| `logFC` | Log2 fold change |
| `logCPM` | Average log2 CPM |
| `F` | Quasi-likelihood F-statistic |
| `PValue` | Raw p-value |
| `FDR` | False discovery rate |


## Interpretation Guidance


### Typical DE Gene Proportions

| Experiment Type | Expected % DE (padj < 0.05, \|LFC\| > 1) |
|----------------|-------------------------------------------|
| Subtle perturbation (low-dose drug, mild stress) | 0.5-3% |
| Standard treatment vs control | 3-10% |
| Different tissues or cell types | 15-40% |
| Cancer vs normal | 10-30% |
| Prokaryotic stress response | 10-50%+ |

If >50% of genes are DE in a standard comparison, suspect a technical issue (batch effect, normalization failure, sample swap). Prokaryotic stress experiments are the exception — bacteria can rewire large portions of their transcriptome.

### LFC Cutoff Selection

| Cutoff | When to Use | Rationale |
|--------|------------|-----------|
| \|LFC\| > 0 (padj only) | Exploratory; generating ranked lists for GSEA | Captures all statistically significant changes |
| \|LFC\| > 0.5 (~1.4-fold) | Default for most experiments | Filters trivially small but statistically significant changes |
| \|LFC\| > 1 (~2-fold) | Standard stringent cutoff | Conventional in literature; good for large-effect studies |
| \|LFC\| > 2 (~4-fold) | Drug screens, very high-signal comparisons | May miss biologically important small changes (e.g., transcription factors) |

Prefer formal threshold testing (`lfcThreshold` in DESeq2, `glmTreat` in edgeR) over post-hoc filtering. Formal tests control the false positive rate at the threshold boundary; post-hoc filtering does not.

### P-value Histogram Diagnostics

Check the raw p-value distribution before trusting DE results:

| Shape | Interpretation | Action |
|-------|---------------|--------|
| Uniform + spike near 0 | Correct: null genes uniform, true DE near 0 | Proceed normally |
| Anti-conservative (U-shape or spike at both ends) | Inflated significance; unmodeled batch or violated assumptions | Check for batch effects, verify model |
| Conservative (spike near 1, depleted near 0) | Over-correction; too many covariates or wrong dispersion | Simplify model, check dispersion plot |
| Spike at p = 1 only | Discrete artifact from low-count genes | Pre-filter more aggressively |

### Shrunken vs Un-shrunken LFCs

| Task | Use |
|------|-----|
| Significance calls (which genes are DE) | Un-shrunken p-values (padj/FDR) |
| Ranking genes by effect size | Shrunken LFCs (apeglm/ashr) |
| GSEA input (ranked gene list) | Shrunken LFCs or Wald statistic |
| Volcano plot x-axis | Shrunken LFCs |
| Post-hoc LFC filtering | Apply to shrunken LFCs for more stable gene lists |

### Preparing Gene Lists for Pathway Analysis

| Method | Input Required | How to Prepare |
|--------|---------------|----------------|
| ORA (enrichGO, enrichKEGG) | Significant gene list + background | `sig_genes <- subset(res, padj < 0.05)`; background = all tested genes |
| GSEA (fgsea, clusterProfiler::GSEA) | ALL genes ranked, no cutoff | Rank by `stat` (DESeq2 Wald) or `sign(logFC) * -log10(PValue)` (edgeR) |

Never use ORA on a ranked list or GSEA on a filtered list. For ORA, always supply the background (all genes that were tested), not just the genome — pre-filtering and independent filtering reduce the tested set.

```r
# GSEA ranking from DESeq2
gsea_ranks <- res_df$stat
names(gsea_ranks) <- res_df$gene
gsea_ranks <- sort(gsea_ranks[!is.na(gsea_ranks)], decreasing = TRUE)

# GSEA ranking from edgeR
gsea_ranks <- sign(results$logFC) * -log10(results$PValue)
names(gsea_ranks) <- rownames(results)
gsea_ranks <- sort(gsea_ranks[is.finite(gsea_ranks)], decreasing = TRUE)
```

### Prokaryotic Gene Annotation

For bacterial/archaeal organisms, Ensembl and org.db packages are unavailable. Use:

```r
# Load annotation from Prokka/Bakta GFF
library(rtracklayer)
gff <- import('annotation.gff3')
gene_info <- as.data.frame(gff[gff$type == 'gene', c('locus_tag', 'Name', 'product')])

# Merge with DE results
res_annotated <- merge(res_df, gene_info, by.x = 'gene', by.y = 'locus_tag', all.x = TRUE)

# KEGG enrichment with bacterial organism code
library(clusterProfiler)
# Find strain-specific KEGG code
search_kegg_organism('Pseudomonas aeruginosa', by = 'scientific_name')
# Use the code (e.g., 'pae' for PAO1)
kegg_res <- enrichKEGG(gene = sig_gene_ids, organism = 'pae', keyType = 'kegg')
```


## Related Skills


- deseq2-basics - Run DESeq2 analysis
- edger-basics - Run edgeR analysis
- de-visualization - Visualize results
- pathway-analysis/go-enrichment - GO over-representation analysis
- pathway-analysis/kegg-pathways - KEGG pathway enrichment
- pathway-analysis/gsea - Gene set enrichment analysis


## Required Libraries


```r
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)  # For labeled points
```


## MA Plot


**Goal:** Visualize the relationship between mean expression and log fold change to assess DE results.

**Approach:** Plot log fold change against mean normalized counts, highlighting significant genes.

**"Make an MA plot of my DE results"** → Plot mean expression vs. fold change with significant genes colored, using plotMA or ggplot2.

### DESeq2 MA Plot

```r
# Built-in MA plot
plotMA(res, ylim = c(-5, 5), main = 'MA Plot')

# With custom alpha
plotMA(res, alpha = 0.05, ylim = c(-5, 5))

# Highlight specific genes
plotMA(res, ylim = c(-5, 5))
with(subset(res, padj < 0.01 & abs(log2FoldChange) > 2),
     points(baseMean, log2FoldChange, col = 'red', pch = 20))
```

### Custom ggplot2 MA Plot

```r
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & !is.na(res_df$padj)

ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c('grey60', 'red')) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(x = 'log10(Mean Expression)', y = 'log2 Fold Change', title = 'MA Plot') +
    theme_bw() +
    theme(legend.position = 'bottom')
```

### edgeR MA Plot

```r
# Using plotMD (mean-difference plot)
plotMD(qlf, main = 'MD Plot')
abline(h = c(-1, 1), col = 'blue', lty = 2)
```


## Volcano Plot


**Goal:** Display statistical significance against fold change magnitude to identify the most important DE genes.

**Approach:** Plot -log10(p-value) vs. log2 fold change with threshold lines and optional gene labels.

**"Create a volcano plot of differentially expressed genes"** → Scatter plot of fold change vs. significance with colored significance regions and labeled top hits.

### Basic Volcano Plot

```r
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c('grey60', 'red')) +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'blue') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'blue') +
    labs(x = 'log2 Fold Change', y = '-log10(p-value)', title = 'Volcano Plot') +
    theme_bw()
```

### Volcano with Gene Labels

```r
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

# Label top genes
top_genes <- head(res_df[order(res_df$padj), ], 10)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significant), alpha = 0.5, size = 1) +
    scale_color_manual(values = c('grey60', 'red')) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 20) +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    labs(x = 'log2 Fold Change', y = '-log10(p-value)') +
    theme_bw()
```

### EnhancedVolcano

```r
library(EnhancedVolcano)

EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = 0.05,
    FCcutoff = 1,
    title = 'Differential Expression',
    subtitle = 'Treatment vs Control')
```


## PCA Plot


**Goal:** Assess sample clustering and identify batch effects or outliers via dimensionality reduction.

**Approach:** Apply variance-stabilizing transformation then project samples onto principal components, coloring by experimental variables.

**"Show me a PCA plot of my samples"** → Perform PCA on transformed expression data and visualize sample separation by condition and batch.

### DESeq2 PCA

```r
# Variance stabilizing transformation first
vsd <- vst(dds, blind = FALSE)

# Basic PCA
plotPCA(vsd, intgroup = 'condition')

# With more options
plotPCA(vsd, intgroup = c('condition', 'batch'), ntop = 500)
```

### Custom PCA with ggplot2

```r
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, 'percentVar'))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
    geom_point(size = 4) +
    xlab(paste0('PC1: ', percentVar[1], '% variance')) +
    ylab(paste0('PC2: ', percentVar[2], '% variance')) +
    ggtitle('PCA Plot') +
    theme_bw() +
    theme(legend.position = 'right')
```

### edgeR PCA (via limma)

```r
library(limma)
log_cpm <- cpm(y, log = TRUE)
plotMDS(log_cpm, col = as.numeric(group), pch = 16)
legend('topright', legend = levels(group), col = 1:nlevels(group), pch = 16)
```


## Heatmaps


**Goal:** Visualize expression patterns of significant genes across samples to reveal clusters and condition effects.

**Approach:** Z-score normalize VST-transformed counts for significant genes and cluster with pheatmap, annotating by condition.

**"Make a heatmap of the top differentially expressed genes"** → Extract significant genes, z-score normalize, and create a clustered heatmap with sample annotations.

### Top DE Genes Heatmap

```r
library(pheatmap)

# Get top significant genes
sig_genes <- rownames(subset(res, padj < 0.01))

# Get normalized counts
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes, ]

# Scale by row (z-score)
mat_scaled <- t(scale(t(mat)))

# Create annotation
annotation_col <- data.frame(
    condition = colData(dds)$condition,
    row.names = colnames(mat)
)

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         color = colorRampPalette(c('blue', 'white', 'red'))(100),
         main = 'Top DE Genes')
```

### Sample Distance Heatmap

```r
vsd <- vst(dds, blind = FALSE)

# Calculate sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Annotation
annotation <- data.frame(
    condition = colData(dds)$condition,
    row.names = colnames(dds)
)

pheatmap(sampleDistMatrix,
         annotation_col = annotation,
         annotation_row = annotation,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(c('white', 'steelblue'))(100),
         main = 'Sample Distance Matrix')
```

### Gene Expression Heatmap

```r
# Select genes of interest
genes_of_interest <- c('gene1', 'gene2', 'gene3', 'gene4', 'gene5')
mat <- assay(vsd)[genes_of_interest, ]

pheatmap(mat,
         scale = 'row',
         annotation_col = annotation_col,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         main = 'Genes of Interest')
```


## Dispersion Plot


**Goal:** Assess the fit of the dispersion model to verify DE analysis assumptions.

**Approach:** Plot gene-wise, fitted, and final dispersion estimates against mean expression.

### DESeq2

```r
plotDispEsts(dds, main = 'Dispersion Estimates')
```

### edgeR

```r
plotBCV(y, main = 'Biological Coefficient of Variation')
```


## Counts Plot for Individual Genes


**Goal:** Visualize expression of a specific gene across samples and conditions.

**Approach:** Extract per-sample counts for a gene and plot by condition using plotCounts or ggplot2.

### DESeq2

```r
# Plot counts for a specific gene
plotCounts(dds, gene = 'GENE_NAME', intgroup = 'condition')

# With ggplot2
d <- plotCounts(dds, gene = 'GENE_NAME', intgroup = 'condition', returnData = TRUE)
ggplot(d, aes(x = condition, y = count, color = condition)) +
    geom_point(position = position_jitter(width = 0.1), size = 3) +
    scale_y_log10() +
    ggtitle('GENE_NAME Expression') +
    theme_bw()
```

### edgeR

```r
# Get CPM for a gene
gene_idx <- which(rownames(y) == 'GENE_NAME')
cpm_gene <- cpm(y)[gene_idx, ]

# Plot
df <- data.frame(cpm = cpm_gene, group = group)
ggplot(df, aes(x = group, y = cpm, color = group)) +
    geom_point(position = position_jitter(width = 0.1), size = 3) +
    scale_y_log10() +
    labs(y = 'CPM', title = 'GENE_NAME Expression') +
    theme_bw()
```


## P-value Histogram


**Goal:** Diagnose the quality of the DE analysis by examining the raw p-value distribution.

**Approach:** Histogram of raw p-values; a uniform distribution with a peak near zero indicates a well-calibrated test.

```r
# Check p-value distribution (should be uniform under null with peak near 0)
res_df <- as.data.frame(res)
ggplot(res_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = 'steelblue', color = 'white') +
    labs(x = 'P-value', y = 'Frequency', title = 'P-value Distribution') +
    theme_bw()
```


## Saving Plots


**Goal:** Export publication-quality plots in vector or raster formats.

**Approach:** Use pdf/png devices or ggsave with appropriate resolution and dimensions.

```r
# Save as PDF (vector)
pdf('volcano_plot.pdf', width = 8, height = 6)
# ... plot code ...
dev.off()

# Save as PNG (raster)
png('volcano_plot.png', width = 800, height = 600, res = 150)
# ... plot code ...
dev.off()

# Using ggsave for ggplot objects
p <- ggplot(...) + ...
ggsave('plot.pdf', p, width = 8, height = 6)
ggsave('plot.png', p, width = 8, height = 6, dpi = 300)
```


## Color Palettes


**Goal:** Select appropriate color schemes for heatmaps and categorical data.

**Approach:** Use RColorBrewer palettes -- diverging for expression, sequential for distances, qualitative for groups.

```r
# For heatmaps
library(RColorBrewer)

# Diverging (for expression: blue-white-red)
colorRampPalette(rev(brewer.pal(n = 7, name = 'RdBu')))(100)

# Sequential (for distances)
colorRampPalette(brewer.pal(n = 9, name = 'Blues'))(100)

# For categorical groups
brewer.pal(n = 8, name = 'Set1')
```


## Quick Reference: Common Plots


| Plot | Purpose | Function |
|------|---------|----------|
| MA plot | LFC vs mean expression | `plotMA()`, `plotMD()` |
| Volcano | LFC vs significance | ggplot2, EnhancedVolcano |
| PCA | Sample clustering | `plotPCA()`, `plotMDS()` |
| Heatmap | Gene patterns | `pheatmap()` |
| Dispersion | Model fit | `plotDispEsts()`, `plotBCV()` |
| Counts | Individual genes | `plotCounts()` |


## Diagnostic Interpretation


### P-value Histogram

Check raw p-value distribution before trusting DE results:

| Shape | Meaning | Action |
|-------|---------|--------|
| Uniform + spike near 0 | Correct: null genes uniform, true DE near 0 | Proceed normally |
| Anti-conservative (U-shape, spikes at 0 and 1) | Inflated significance; unmodeled batch or violated assumptions | Check for batch effects, verify model specification |
| Conservative (depleted near 0, spike near 1) | Over-correction; too many covariates or wrong dispersion | Simplify model, check dispersion plot |
| Spike at p = 1 only | Discrete artifact from low-count genes | Pre-filter more aggressively |

### MA Plot Diagnostics

| Pattern | Meaning |
|---------|---------|
| Symmetric cloud centered at LFC = 0 | Correct normalization |
| Cloud shifted up or down | Normalization failure; majority-DE experiment may violate assumptions |
| Funnel shape widening at low expression | Expected — low-count genes have noisier fold changes |
| Discrete horizontal bands | Low-count artifacts; consider stronger pre-filtering |

### Volcano Plot: Shrunken LFCs

Use shrunken LFCs (apeglm/ashr) on the x-axis and un-shrunken p-values on the y-axis. This combination gives stable fold change estimates while preserving the original significance assessment. Without shrinkage, low-count genes with extreme but unreliable fold changes dominate the plot edges.

### Dispersion Plot Diagnostics

| Pattern | Meaning |
|---------|---------|
| Gene-wise points scattered around fitted line | Good model fit |
| Gene-wise points far above fitted line | Possible outlier genes or unmodeled batch effects |
| Fitted line flat (no trend) | Unusual — check if data is over-filtered or has unusual structure |
| Final estimates much lower than gene-wise | Expected — shrinkage toward the fitted trend |

### PCA Diagnostics

| Pattern | Meaning | Action |
|---------|---------|--------|
| Clear separation by condition on PC1/PC2 | Strong biological signal | Proceed |
| Separation by batch, not condition | Batch effect dominates | Include batch in model; do NOT use corrected counts for DE |
| One sample far from its group | Potential outlier or sample swap | Check library QC metrics; consider removing |
| No separation on PC1/PC2 but present on PC3+ | Subtle effects | May still find DE genes; check dispersion estimates |


## Related Skills


- deseq2-basics - Generate DESeq2 results for visualization
- edger-basics - Generate edgeR results for visualization
- de-results - Filter genes before visualization
- data-visualization/specialized-omics-plots - Custom ggplot2 volcano/MA/PCA functions
- data-visualization/heatmaps-clustering - Advanced heatmap customization

