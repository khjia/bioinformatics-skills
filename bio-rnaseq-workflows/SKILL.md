---
name: bio-rnaseq-workflows
description: "End-to-end RNA-seq analysis workflows: from raw FASTQ to DEGs (Salmon+DESeq2 or STAR+featureCounts), time-course analysis, splicing analysis, small RNA-seq, and expression-to-pathways integration."
---


## Workflow Overview


```
FASTQ files
    |
    v
[1. QC & Trimming] -----> fastp
    |
    v
[2. Quantification] ----> Salmon (recommended) or STAR + featureCounts
    |
    v
[3. Import to R] -------> tximport (for Salmon) or direct counts
    |
    v
[4. DE Analysis] -------> DESeq2
    |
    v
[5. Visualization] -----> Volcano, MA, heatmaps
    |
    v
Significant gene list
```


## Experimental Design & Reproducibility


A defensible RNA-seq analysis starts with **design**, not tools. These decisions determine whether results are trustworthy:

### Replication
- **≥3 biological replicates per group** is the minimum for dispersion estimation
- More replicates beat deeper sequencing for statistical power
- Technical replicates (same library resequenced) don't count — DESeq2 treats them as biological
- No-replicate designs are **not supported** by DESeq2 (removed since v1.22)

### Batch & Confounding
- **Confounded batch and condition** is unrecoverable — if every treated sample was processed on a different day/lane than controls, the effect can't be separated from batch
- Randomize sample processing across conditions
- Model known batches in the design formula: `~ batch + condition`

### Strandedness
- Choosing the wrong strandedness silently discards ~half the reads
- Use Salmon `-l A` (auto-detect) or infer strandedness with `infer_experiment.py` (RSeQC)
- Verify the assigned-reads fraction after quantification

### Reproducibility
- **Pin tool versions**: record exact versions of fastp, STAR/Salmon, DESeq2, reference genome + annotation release
- These belong in the Methods section
- Unpinned "latest" pipelines make results irreproducible

### Gene-ID Consistency
- DESeq2 output is often Ensembl IDs; enrichment tools (Enrichr, MSigDB) want symbols
- Map IDs **before** enrichment, or "nothing is significant"
- Record the ID mapping source and version


## Primary Path: Salmon + DESeq2


### Step 1: Quality Control with fastp

```bash
# Single sample
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o sample_R1.trimmed.fq.gz -O sample_R2.trimmed.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 35 \
    --html sample_fastp.html

# Batch processing
for sample in sample1 sample2 sample3; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o trimmed/${sample}_R1.fq.gz -O trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --html qc/${sample}_fastp.html
done
```

**QC Checkpoint 1:** Check fastp reports
- Q30 bases >80%
- Adapter content <5%
- Duplication rate reasonable for library type

### Step 2: Salmon Quantification

```bash
# Build index (once per transcriptome)
salmon index -t transcriptome.fa -i salmon_index -k 31

# Quantify each sample
for sample in sample1 sample2 sample3; do
    salmon quant -i salmon_index \
        -l A \
        -1 trimmed/${sample}_R1.fq.gz \
        -2 trimmed/${sample}_R2.fq.gz \
        -o quants/${sample} \
        --validateMappings \
        --gcBias \
        --seqBias \
        -p 8
done
```

**QC Checkpoint 2:** Check Salmon logs
- Mapping rate >70%
- >10 million reads mapped

### Step 3: Import with tximport

```r
library(tximport)
library(DESeq2)

# Create tx2gene mapping (Ensembl example)
tx2gene <- read.csv('tx2gene.csv')  # columns: TXNAME, GENEID

# List quantification files
samples <- c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6')
files <- file.path('quants', samples, 'quant.sf')
names(files) <- samples

# Import transcript-level estimates
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

# Create sample metadata
coldata <- data.frame(
    condition = factor(c('control', 'control', 'control', 'treated', 'treated', 'treated')),
    row.names = samples
)
```

### Step 4: DESeq2 Analysis

```r
# Create DESeqDataSet from tximport
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2
dds <- DESeq(dds)

# Get results with shrinkage
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Summary
summary(res)
```

**QC Checkpoint 3:** Check DESeq2 diagnostics
- Dispersion plot shows expected trend
- PCA separates conditions
- No severe outliers in sample distances

### Step 5: Visualization and Export

```r
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c('grey', 'red')) +
    theme_minimal() +
    labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 P-value')

# Heatmap of top genes
vsd <- vst(dds, blind = FALSE)
top_genes <- head(order(res$padj), 50)
pheatmap(assay(vsd)[top_genes,], scale = 'row', show_rownames = FALSE)

# Export significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes), 'significant_genes.csv')
```


## Alternative Path: STAR + featureCounts + DESeq2


### Step 2 Alternative: STAR Alignment

```bash
# Build STAR index (once)
STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf \
    --sjdbOverhang 100 \
    --runThreadN 8

# Align each sample
for sample in sample1 sample2 sample3; do
    STAR --genomeDir star_index \
        --readFilesIn trimmed/${sample}_R1.fq.gz trimmed/${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix aligned/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --runThreadN 8
done
```

### Step 3 Alternative: featureCounts

```bash
# Count reads per gene
featureCounts -T 8 -p --countReadPairs \
    -a genes.gtf \
    -o counts.txt \
    aligned/*_Aligned.sortedByCoord.out.bam
```

### Step 4 Alternative: Load Counts Directly

```r
# Load featureCounts output
counts <- read.table('counts.txt', header = TRUE, row.names = 1, skip = 1)
counts <- counts[, 6:ncol(counts)]  # Remove annotation columns
colnames(counts) <- gsub('_Aligned.sortedByCoord.out.bam', '', colnames(counts))

# Create DESeqDataSet directly
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
```


## Parameter Recommendations


| Step | Parameter | Recommendation |
|------|-----------|----------------|
| fastp | --qualified_quality_phred | 20 (standard) |
| fastp | --length_required | 35 for 2x100, 50 for 2x150 |
| Salmon | -l | A (auto-detect library type) |
| Salmon | --gcBias | Enable for better accuracy |
| STAR | --sjdbOverhang | read_length - 1 |
| featureCounts | -s | 0=unstranded, 1=stranded, 2=reversely stranded |
| DESeq2 | lfcShrink type | apeglm (recommended) |
| DESeq2 | alpha | 0.05 (standard significance) |


## Common Pitfalls


These cause most wrong or irreproducible bulk RNA-seq results:

| # | Pitfall | Why It Matters | Prevention |
|---|---------|---------------|------------|
| 1 | **Too few replicates** | <3 biological replicates gives almost no power and unstable dispersion estimates | ≥3 biological replicates per group; more replicates > deeper sequencing |
| 2 | **Confounded batch and condition** | If every treated sample was processed on a different day/lane than controls, the effect is unrecoverable | Randomize sample processing; model known batches (`~batch + condition`) |
| 3 | **Wrong strandedness** | Choosing the wrong STAR column or featureCounts `-s`/Salmon library type silently discards ~half the reads | Use Salmon `-l A` or infer strandedness; verify assigned-reads fraction |
| 4 | **Feeding TPM/FPKM to DESeq2** | DESeq2 needs raw (or length-scaled) **counts**, never TPM/FPKM/normalized values | Use tximport with `countsFromAbundance="lengthScaledTPM"` or raw featureCounts |
| 5 | **Gene-ID mismatch into enrichment** | DESeq2 output is often Ensembl IDs; Enrichr/MSigDB want symbols → "nothing is significant" | Map IDs before enrichment; verify match rate >85% |
| 6 | **Skipping post-quant QC** | Swapped labels, outliers, and hidden batches are invisible until you look at PCA and sample-distance heatmaps | Always run PCA + sample-distance heatmap before trusting DE results |
| 7 | **Mixing tools across samples** | Different aligners/versions/references across samples introduces systematic bias | Quantify every sample with the same tool, version, reference, and parameters |
| 8 | **Low mapping rate (<50%)** | Wrong reference genome, contamination, or severe degradation | Check species, run FastQ Screen; investigate before proceeding |
| 9 | **All genes DE / no genes DE** | Usually indicates normalization failure, sample swap, or confounded design | Check sample metadata, PCA, dispersion plot; verify design matrix is full rank |

### Post-Quantification QC (Critical)
Always inspect before trusting DE results:
- **PCA plot**: should separate by condition, not by batch/date/lane
- **Sample-distance heatmap**: replicates should cluster together
- **Dispersion plot**: should show expected decreasing trend with mean
- **P-value histogram**: should be mostly flat with a peak near 0 (if real signal exists)

```r
# Post-quant QC
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists))
plotDispEsts(dds)
hist(res$pvalue, breaks=50, main="P-value distribution")
```


## Complete Bash Pipeline Script


```bash
#!/bin/bash
set -e

THREADS=8
SAMPLES="sample1 sample2 sample3 sample4 sample5 sample6"
SALMON_INDEX="salmon_index"
OUTDIR="results"

mkdir -p ${OUTDIR}/{trimmed,quants,qc}

# Step 1: QC and trim
for sample in $SAMPLES; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        -w ${THREADS}
done

# Step 2: Quantify
for sample in $SAMPLES; do
    salmon quant -i ${SALMON_INDEX} -l A \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        -o ${OUTDIR}/quants/${sample} \
        --validateMappings --gcBias -p ${THREADS}
done

echo "Quantification complete. Run R script for DE analysis."
```


## Complete R Analysis Script


```r
library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)

# Configuration
samples <- c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6')
conditions <- c('control', 'control', 'control', 'treated', 'treated', 'treated')
quant_dir <- 'results/quants'

# Import
tx2gene <- read.csv('tx2gene.csv')
files <- file.path(quant_dir, samples, 'quant.sf')
names(files) <- samples
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

# DESeq2
coldata <- data.frame(condition = factor(conditions), row.names = samples)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)

# Results
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

cat('Significant genes:', nrow(sig), '\n')
write.csv(as.data.frame(sig), 'significant_genes.csv')
```


## Related Skills


- read-qc/fastp-workflow - Detailed QC options and parameters
- rna-quantification/alignment-free-quant - Salmon and kallisto details
- rna-quantification/tximport-workflow - tximport options and tx2gene creation
- differential-expression/deseq2-basics - Complete DESeq2 reference
- differential-expression/de-visualization - Advanced visualization options
- pathway-analysis/go-enrichment - Next step: functional enrichment


## Pipeline Overview


```
Expression matrix + time metadata
    |
    v
[1. Temporal DE] ---------> limma splines / DESeq2 LRT
    |
    v
[2. Filter] --------------> Significant temporal genes (FDR <0.05)
    |
    v
[3. Mfuzz Clustering] ----> Soft clustering of expression profiles
    |                            |
    |                            +---> QC: membership >0.5, no empty clusters
    |
    +--- Circadian design? ---> [4a. Rhythm Detection] (MetaCycle / CosinorPy)
    |                               |
    |                               v
    |                           Rhythmic genes + period/phase estimates
    |
    v
[4b. GAM Trajectory] -----> mgcv GAM fitting for top clusters
    |
    v
[5. Pathway Enrichment] --> clusterProfiler per-cluster GO/KEGG
    |
    v
Temporal gene modules + enriched pathways + trajectory plots
```


## Step 1: Temporal Differential Expression


### R (limma splines)

```r
library(limma)
library(splines)

expr <- as.matrix(read.csv('counts_normalized.csv', row.names = 1))
meta <- read.csv('metadata.csv')

time_points <- meta$time
design <- model.matrix(~ ns(time_points, df = 3))

fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Test all spline coefficients jointly for temporal significance
temporal_results <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = 'F')
# topTable already returns adj.P.Val (BH-corrected); use it directly
```

### R (DESeq2 LRT)

```r
library(DESeq2)

counts <- as.matrix(read.csv('raw_counts.csv', row.names = 1))
meta <- read.csv('metadata.csv')
meta$time <- factor(meta$time)

dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ 1)
res <- results(dds)
```

### Python (statsmodels)

```python
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from scipy import stats

expr = pd.read_csv('counts_normalized.csv', index_col=0)
meta = pd.read_csv('metadata.csv')

spline_basis = dmatrix('bs(time, df=3)', data=meta, return_type='dataframe')
design_full = np.column_stack([np.ones(len(meta)), spline_basis.values])
design_reduced = np.ones((len(meta), 1))

pvals = []
for gene in expr.index:
    y = expr.loc[gene].values
    ss_full = np.sum((y - design_full @ np.linalg.lstsq(design_full, y, rcond=None)[0]) ** 2)
    ss_red = np.sum((y - design_reduced @ np.linalg.lstsq(design_reduced, y, rcond=None)[0]) ** 2)
    df_diff = design_full.shape[1] - design_reduced.shape[1]
    df_resid = len(y) - design_full.shape[1]
    f_stat = ((ss_red - ss_full) / df_diff) / (ss_full / df_resid)
    pvals.append(1 - stats.f.cdf(f_stat, df_diff, df_resid))

_, fdr, _, _ = multipletests(pvals, method='fdr_bh')
temporal_genes = expr.index[fdr < 0.05].tolist()
```

### QC Checkpoint: Temporal DE

```r
# Gate 1: Sufficient temporal genes detected
sig_genes <- temporal_results[temporal_results$adj.P.Val < 0.05, ]
n_sig <- nrow(sig_genes)
message(sprintf('Significant temporal genes: %d', n_sig))
if (n_sig < 100) message('WARNING: Few temporal genes. Check time point spacing or consider relaxing FDR.')
if (n_sig > 10000) message('WARNING: Many temporal genes. Consider stricter FDR or inspect batch effects.')

# Gate 2: Residual distribution check
residuals <- residuals(fit, expr)
message(sprintf('Residual mean: %.4f, SD: %.4f', mean(residuals), sd(residuals)))
```


## Step 2: Filter Significant Genes


```r
# FDR <0.05: Standard threshold for temporal DE
# More permissive (0.1) acceptable for exploratory clustering
sig_genes <- rownames(temporal_results[temporal_results$adj.P.Val < 0.05, ])
expr_sig <- expr[sig_genes, ]
message(sprintf('Genes passing FDR <0.05: %d', length(sig_genes)))
```


## Step 3: Mfuzz Soft Clustering


```r
library(Mfuzz)

eset <- ExpressionSet(assayData = as.matrix(expr_sig))

# Standardize expression profiles (mean=0, sd=1 per gene)
eset <- standardise(eset)

# Estimate fuzzifier m
# mestimate() calculates optimal m from data geometry; typical range 1.5-2.5
m <- mestimate(eset)
message(sprintf('Estimated fuzzifier m = %.2f', m))

# Cluster count: start with sqrt(n_genes/2), refine with gap statistic
# Typical range 4-20 depending on temporal complexity
n_clusters <- 8
cl <- mfuzz(eset, c = n_clusters, m = m)

# Filter low-membership genes
# Membership >0.5: gene clearly belongs to one cluster
# Lower to 0.3 for exploratory analysis with more overlap
core_genes <- acore(eset, cl, min.acore = 0.5)
```

### Python Alternative (tslearn)

```python
from tslearn.clustering import TimeSeriesKMeans

# Row-wise z-scoring: normalize each gene across its timepoints (not per-timepoint)
expr_scaled = (expr_sig.values - expr_sig.values.mean(axis=1, keepdims=True)) / expr_sig.values.std(axis=1, keepdims=True)

# n_clusters: 4-20 depending on complexity; evaluate with silhouette score
model = TimeSeriesKMeans(n_clusters=8, metric='softdtw', metric_params={'gamma': 0.01},
                         max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled.reshape(expr_scaled.shape[0], expr_scaled.shape[1], 1))
```

### QC Checkpoint: Clustering

```r
# Gate 1: No empty clusters
cluster_sizes <- table(cl$cluster)
message('Cluster sizes:')
print(cluster_sizes)
if (any(cluster_sizes == 0)) message('WARNING: Empty clusters found. Reduce n_clusters.')

# Gate 2: Membership quality
for (i in seq_along(core_genes)) {
    n_core <- nrow(core_genes[[i]])
    message(sprintf('Cluster %d: %d core genes (membership >0.5)', i, n_core))
}

# Gate 3: Silhouette score (optional validation)
library(cluster)
hard_labels <- cl$cluster
sil <- silhouette(hard_labels, dist(exprs(eset)))
message(sprintf('Mean silhouette: %.3f', mean(sil[, 3])))
```


## Step 4a: Rhythm Detection (Optional - Circadian Designs)


Only applicable when sampling covers 24h+ cycles with sufficient resolution (every 2-4h).

### R (MetaCycle)

```r
library(MetaCycle)

# Expects genes as rows, time points as columns
# Column names must be numeric time values (hours)
expr_for_meta <- expr_sig
colnames(expr_for_meta) <- meta$time_hours

write.csv(expr_for_meta, 'expr_for_metacycle.csv')

# Period range 20-28h: standard circadian search window
# Adjust for ultradian (4-12h) or infradian (>28h) rhythms
meta2d('expr_for_metacycle.csv', filestyle = 'csv',
       minper = 20, maxper = 28,
       timepoints = sort(unique(meta$time_hours)),
       outdir = 'metacycle_results')
```

### Python (CosinorPy)

```python
from cosinorpy import file_parser, cosinor

# fit_group expects long-format DataFrame with columns 'x' (time), 'y' (expression), 'test' (gene name)
# Reshape expression matrix to long format before passing
# period=24: standard circadian; adjust for other periodicities
results = cosinor.fit_group(expr_long, period=24, n_components=1)
rhythmic = results[results['p'] < 0.05]
```


## Step 4b: GAM Trajectory Fitting


### R (mgcv)

```r
library(mgcv)

cluster_trajectories <- list()
for (cl_id in 1:n_clusters) {
    cl_genes <- names(cl$cluster[cl$cluster == cl_id])
    mean_profile <- colMeans(expr_sig[cl_genes, ])

    df_gam <- data.frame(time = meta$time, expr = mean_profile)

    # k: basis dimension; k=5 sufficient for most time courses
    # Increase to k=10 for >20 time points; decrease to k=3 for <6 time points
    gam_fit <- gam(expr ~ s(time, k = 5), data = df_gam)

    cluster_trajectories[[cl_id]] <- list(
        fit = gam_fit,
        r_squared = summary(gam_fit)$r.sq,
        edf = summary(gam_fit)$edf
    )
    message(sprintf('Cluster %d: R^2 = %.3f, EDF = %.2f', cl_id,
                    summary(gam_fit)$r.sq, summary(gam_fit)$edf))
}
```

### Python (pygam)

```python
from pygam import LinearGAM, s
import numpy as np

for cl_id in range(n_clusters):
    cl_mask = labels == cl_id
    mean_profile = expr_scaled[cl_mask].mean(axis=0)

    # n_splines=5: sufficient for most time courses
    gam = LinearGAM(s(0, n_splines=5)).fit(meta['time'].values.reshape(-1, 1), mean_profile)
    print(f'Cluster {cl_id}: GCV = {gam.statistics_["GCV"]:.4f}')
```


## Step 5: Per-Cluster Pathway Enrichment


### R (clusterProfiler)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# Background = all temporal genes tested (not the full genome)
all_temporal_entrez <- bitr(rownames(expr_sig), fromType = 'SYMBOL', toType = 'ENTREZID',
                            OrgDb = org.Hs.eg.db)

enrichment_results <- list()
for (i in seq_along(core_genes)) {
    genes <- core_genes[[i]]$NAME

    entrez <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

    # GO Biological Process enrichment with proper background
    ego <- enrichGO(gene = entrez$ENTREZID,
                    universe = all_temporal_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = 'BP', pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                    readable = TRUE)

    # Simplify redundant GO terms (parent-child hierarchy creates redundancy)
    if (nrow(as.data.frame(ego)) > 0) {
        ego <- simplify(ego, cutoff = 0.7, by = 'p.adjust')
    }

    enrichment_results[[i]] <- ego
    message(sprintf('Cluster %d: %d significant GO terms', i, nrow(as.data.frame(ego))))
}
```

### Python (gseapy)

```python
import gseapy as gp

# Use all temporal genes as background for proper enrichment statistics
all_temporal_genes = list(expr_sig.index)

for cl_id in range(n_clusters):
    cl_genes = [g for g, l in zip(expr_sig.index, labels) if l == cl_id]

    enr = gp.enrichr(gene_list=cl_genes, gene_sets='GO_Biological_Process_2023',
                     organism='human', background=all_temporal_genes,
                     outdir=f'enrichr_cluster_{cl_id}')
    sig_terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
    print(f'Cluster {cl_id}: {len(sig_terms)} significant GO terms')
```

### QC Checkpoint: Enrichment

```r
# Gate: At least 3 clusters should have significant GO terms
clusters_with_terms <- sum(sapply(enrichment_results, function(x) nrow(as.data.frame(x)) > 0))
message(sprintf('Clusters with significant GO terms: %d / %d', clusters_with_terms, length(enrichment_results)))
if (clusters_with_terms < 3) {
    message('WARNING: Few clusters enriched. Check gene ID conversion or relax thresholds.')
}
```


## Troubleshooting


| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| < 100 temporal genes | Insufficient replicates or noisy data | Add replicates; use DESeq2 LRT instead of limma |
| Empty Mfuzz clusters | Too many clusters | Reduce n_clusters; check gap statistic |
| All genes in one cluster | Fuzzifier too low or too few clusters | Increase m or n_clusters |
| No rhythmic genes | Non-circadian design or low power | Verify 24h+ sampling; increase resolution |
| GAM overfitting | k too high for time points | Set k = min(n_timepoints - 1, 5) |
| Few enriched clusters | Gene ID conversion failure | Check species; verify Entrez ID mapping |
| Low membership scores | High expression noise | Increase fuzzifier m; apply stricter gene filtering |


## Step 1: Read Quality Control


```bash
# fastp for adapter trimming and quality filtering
fastp \
    -i sample_R1.fastq.gz \
    -I sample_R2.fastq.gz \
    -o sample_clean_R1.fastq.gz \
    -O sample_clean_R2.fastq.gz \
    --detect_adapter_for_pe \
    --thread 8 \
    -h sample_fastp.html
```


## Step 2: STAR 2-Pass Alignment


```bash
# First pass to detect novel junctions
STAR \
    --runThreadN 8 \
    --genomeDir star_index/ \
    --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix sample_pass1_ \
    --outSAMtype BAM Unsorted \
    --outSJfilterOverhangMin 8 8 8 8 \
    --alignSJDBoverhangMin 1

# Generate new index with discovered junctions
# (Combine SJ.out.tab files from all samples)
cat *_SJ.out.tab > combined_SJ.out.tab

# Second pass with combined junctions
STAR \
    --runThreadN 8 \
    --genomeDir star_index/ \
    --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
    --readFilesCommand zcat \
    --sjdbFileChrStartEnd combined_SJ.out.tab \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSJfilterOverhangMin 8 8 8 8 \
    --alignSJDBoverhangMin 1 \
    --quantMode GeneCounts
```


## Step 3: Junction QC Checkpoint


```python
import subprocess

def check_junction_saturation(bam_file, bed_file, output_prefix):
    '''
    QC Checkpoint: Verify junction detection saturation.
    Plateau indicates sufficient depth for splicing analysis.
    '''
    subprocess.run([
        'junction_saturation.py',
        '-i', bam_file,
        '-r', bed_file,
        '-o', output_prefix
    ], check=True)

    # Manual check: curves should plateau
    print(f'Check {output_prefix}.junctionSaturation_plot.pdf')
    print('If curves still rising, consider deeper sequencing')
```


## Step 4: Differential Splicing with rMATS-turbo


```bash
# Create sample list files
# condition1_bams.txt: sample1.bam,sample2.bam,sample3.bam
# condition2_bams.txt: sample4.bam,sample5.bam,sample6.bam

rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp
```


## Step 5: Filter Results


```python
import pandas as pd

def filter_differential_splicing(rmats_dir, event_type='SE',
                                  fdr_cutoff=0.05, dpsi_cutoff=0.1, min_reads=10):
    '''
    Filter rMATS results for significant events.

    Thresholds:
    - |deltaPSI| > 0.1 (lenient) or > 0.2 (stringent)
    - FDR < 0.05
    - Junction reads >= 10
    '''
    jc_file = f'{rmats_dir}/{event_type}.MATS.JC.txt'
    df = pd.read_csv(jc_file, sep='\t')

    significant = df[
        (df['FDR'] < fdr_cutoff) &
        (df['IncLevelDifference'].abs() > dpsi_cutoff)
    ].copy()

    print(f'Significant {event_type} events: {len(significant)}')

    # Sort by significance and effect size
    significant['score'] = -significant['FDR'].apply(lambda x: max(x, 1e-300)).apply(
        lambda x: __import__('numpy').log10(x)
    ) * significant['IncLevelDifference'].abs()

    return significant.sort_values('score', ascending=False)
```


## Step 6: Optional Isoform Switching


```r
library(IsoformSwitchAnalyzeR)

# Import Salmon quantification if available
switchList <- importRdata(
    isoformCountMatrix = counts,
    isoformRepExpression = tpm,
    designMatrix = design,
    isoformExonAnnoation = 'annotation.gtf',
    isoformNtFasta = 'transcripts.fa'
)

# Analyze switches
switchList <- isoformSwitchTestDEXSeq(switchList, reduceToSwitchingGenes = TRUE)
```


## Step 7: Sashimi Visualization


```python
import subprocess

def visualize_top_events(rmats_dir, grouping_file, gtf_file, output_dir, n_top=20):
    '''Generate sashimi plots for top differential events.'''
    import pandas as pd
    from pathlib import Path

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for event_type in ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']:
        jc_file = f'{rmats_dir}/{event_type}.MATS.JC.txt'
        df = pd.read_csv(jc_file, sep='\t')
        sig = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.1)]

        for idx, event in sig.head(n_top).iterrows():
            chrom = event['chr']
            start = event.get('upstreamES', event.get('1stExonStart_0base', 0)) - 500
            end = event.get('downstreamEE', event.get('2ndExonEnd', 0)) + 500
            gene = event['geneSymbol']

            subprocess.run([
                'ggsashimi.py',
                '-b', grouping_file,
                '-c', f'{chrom}:{start}-{end}',
                '-o', f'{output_dir}/{event_type}_{gene}',
                '-g', gtf_file,
                '--shrink',
                '--fix-y-scale',
                '-M', '5'
            ], check=True)
```


## Complete Pipeline Script


```bash
#!/bin/bash
set -e

# Configuration
SAMPLES="sample1 sample2 sample3 sample4 sample5 sample6"
CONDITIONS="control control control treatment treatment treatment"
GTF="annotation.gtf"
STAR_INDEX="star_index/"
THREADS=8

# Step 1: QC and trimming
for sample in $SAMPLES; do
    fastp -i ${sample}_R1.fq.gz -I ${sample}_R2.fq.gz \
          -o ${sample}_clean_R1.fq.gz -O ${sample}_clean_R2.fq.gz \
          --thread $THREADS
done

# Step 2: STAR 2-pass alignment
# ... (as above)

# Step 3: Junction QC
for sample in $SAMPLES; do
    junction_saturation.py -i ${sample}.bam -r annotation.bed -o ${sample}_junc
done

# Step 4: rMATS differential splicing
rmats.py --b1 control_bams.txt --b2 treatment_bams.txt \
         --gtf $GTF -t paired --readLength 150 --nthread $THREADS \
         --od rmats_output --tmp rmats_tmp

echo "Pipeline complete. Check rmats_output/ for results."
```


## When NOT to Use This Pipeline (Pipeline Variants)


This pipeline targets **bulk short-read RNA-seq differential splicing between two groups**. For other regimes, use the dedicated skill:

| Question | Use instead |
|----------|-------------|
| "Does this DNA variant alter splicing?" | alternative-splicing/splice-variant-prediction (SpliceAI, Pangolin, MMSplice, ClinGen SVI 2023) |
| "What is aberrant in this single rare-disease patient?" | alternative-splicing/outlier-splicing-detection (FRASER 2.0, OUTRIDER, DROP) |
| "Full-isoform analysis from PacBio Iso-Seq / ONT" | alternative-splicing/long-read-splicing (FLAIR, IsoQuant, Bambu, SQANTI3, rMATS-long) |
| "Single-cell splicing analysis" | alternative-splicing/single-cell-splicing (chemistry-first decision; MARVEL, BRIE2 plate; long-read SC) |
| "Heterogeneous cohort, n>=10 vs n>=10" | This pipeline + MAJIQ V3 HET module (see alternative-splicing/differential-splicing) |
| "Microexon-focused (3-27 nt)" | This pipeline with VAST-TOOLS or MicroExonator; see alternative-splicing/splicing-quantification |


## Step 1: Preprocessing


```bash
# Adapter trimming and size selection
cutadapt -a TGGAATTCTCGGGTGCCAAGG \
    --minimum-length 18 --maximum-length 30 \
    -o trimmed.fastq.gz reads.fastq.gz
```


## Step 2: miRDeep2 Analysis


```bash
# Align to genome
mapper.pl trimmed.fastq.gz -e -h -i -j -l 18 \
    -m -p genome_index -s reads_collapsed.fa \
    -t reads_collapsed_vs_genome.arf

# miRNA quantification and novel prediction
miRDeep2.pl reads_collapsed.fa genome.fa \
    reads_collapsed_vs_genome.arf \
    mature_ref.fa none hairpin_ref.fa
```


## Step 3: Differential Expression


```r
library(DESeq2)
counts <- read.csv('mirna_counts.csv', row.names = 1)
dds <- DESeqDataSetFromMatrix(counts, colData, ~condition)
dds <- DESeq(dds)
results <- results(dds)
```


## Step 4: Target Prediction


```bash
# miRanda for target prediction
miranda mature_mirnas.fa target_3utrs.fa -out targets.txt
```


## QC Checkpoints


1. **After trimming**: Size distribution should peak at 21-23nt
2. **After alignment**: >70% mapping rate expected
3. **After DE**: Check volcano plot and PCA


## Method Selection


| Scenario | Method | Why |
|----------|--------|-----|
| Have DE results with Wald stat / t-stat for all genes | GSEA (Step 4) | Uses full ranking; no arbitrary cutoff; ~35% higher F1 than ORA |
| Clear gene list from non-DE source (co-expression, GWAS) | ORA (Steps 1-3) | No ranking available |
| RNA-seq with known gene length bias | GOseq (goseq package) | Standard ORA ignores length bias |
| Bacterial / prokaryotic data | KEGG with locus tags | No org.*.eg.db; use keyType='kegg' |
| Multiple conditions to compare | compareCluster or mitch | Never compare p-values across separate enrichments |

When in doubt, run both ORA and GSEA and compare. Concordant results are more trustworthy.


## Input Preparation


### From DESeq2 Results

```r
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load DE results
res <- read.csv('deseq2_results.csv', row.names = 1)

# Significant genes for ORA
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))

# Background = all tested genes (NOT the full genome)
# Pre-filtering and independent filtering reduce the tested set; use only genes that were tested
background_genes <- rownames(res[!is.na(res$pvalue), ])

# Ranked list for GSEA — prefer Wald statistic (combines magnitude + precision)
# Alternatives: shrunken LFC, or sign(logFC) * -log10(PValue) for edgeR
ranked_genes <- res$stat
names(ranked_genes) <- rownames(res)
ranked_genes <- sort(ranked_genes[!is.na(ranked_genes)], decreasing = TRUE)
```

### Gene ID Conversion

```r
# Convert gene symbols to Entrez IDs
sig_entrez <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID',
                   OrgDb = org.Hs.eg.db)

# For ranked list
ranked_entrez <- bitr(names(ranked_genes), fromType = 'SYMBOL', toType = 'ENTREZID',
                      OrgDb = org.Hs.eg.db)
ranked_list <- ranked_genes[ranked_entrez$SYMBOL]
names(ranked_list) <- ranked_entrez$ENTREZID
```


## Step 1: GO Over-representation Analysis


```r
# Convert background genes too
bg_entrez <- bitr(background_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Biological Process — always specify universe (background)
go_bp <- enrichGO(gene = sig_entrez$ENTREZID,
                  universe = bg_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = 'BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1,
                  readable = TRUE)

# Molecular Function
go_mf <- enrichGO(gene = sig_entrez$ENTREZID,
                  universe = bg_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = 'MF',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# Cellular Component
go_cc <- enrichGO(gene = sig_entrez$ENTREZID,
                  universe = bg_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = 'CC',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# Simplify redundant terms
go_bp_simple <- simplify(go_bp, cutoff = 0.7, by = 'p.adjust')
```


## Step 2: KEGG Pathway Enrichment


```r
kegg <- enrichKEGG(gene = sig_entrez$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)

# Convert KEGG IDs to readable names
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```


## Step 3: Reactome Pathway Enrichment


```r
library(ReactomePA)

reactome <- enrichPathway(gene = sig_entrez$ENTREZID,
                          organism = 'human',
                          pvalueCutoff = 0.05,
                          readable = TRUE)
```


## Step 4: Gene Set Enrichment Analysis (GSEA)


```r
# GO GSEA
gsea_go <- gseGO(geneList = ranked_list,
                 OrgDb = org.Hs.eg.db,
                 ont = 'BP',
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

# KEGG GSEA
gsea_kegg <- gseKEGG(geneList = ranked_list,
                     organism = 'hsa',
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)
```


## Step 5: Visualization


```r
library(enrichplot)
library(ggplot2)

# Dot plot
dotplot(go_bp_simple, showCategory = 20) +
    ggtitle('GO Biological Process Enrichment')
ggsave('go_bp_dotplot.pdf', width = 10, height = 8)

# Bar plot
barplot(kegg, showCategory = 15) +
    ggtitle('KEGG Pathway Enrichment')
ggsave('kegg_barplot.pdf', width = 9, height = 6)

# Enrichment map (network of related terms)
go_bp_simple <- pairwise_termsim(go_bp_simple)
emapplot(go_bp_simple, showCategory = 30) +
    ggtitle('GO Term Similarity Network')
ggsave('go_network.pdf', width = 10, height = 10)

# Concept network (gene-term connections)
cnetplot(go_bp, showCategory = 5, categorySize = 'pvalue') +
    ggtitle('Gene-Concept Network')
ggsave('cnet_plot.pdf', width = 12, height = 10)

# GSEA plot for specific pathway
gseaplot2(gsea_kegg, geneSetID = 1:3, pvalue_table = TRUE)
ggsave('gsea_plot.pdf', width = 10, height = 8)

# Ridge plot for GSEA
ridgeplot(gsea_go, showCategory = 15)
ggsave('gsea_ridge.pdf', width = 8, height = 10)
```


## Step 6: Export Results


```r
# Export enrichment results
write.csv(as.data.frame(go_bp), 'go_bp_enrichment.csv', row.names = FALSE)
write.csv(as.data.frame(kegg), 'kegg_enrichment.csv', row.names = FALSE)
write.csv(as.data.frame(reactome), 'reactome_enrichment.csv', row.names = FALSE)
write.csv(as.data.frame(gsea_go), 'gsea_go_results.csv', row.names = FALSE)

# Combine key results
combined <- rbind(
    data.frame(Database = 'GO_BP', as.data.frame(go_bp_simple)[1:10,]),
    data.frame(Database = 'KEGG', as.data.frame(kegg)[1:10,]),
    data.frame(Database = 'Reactome', as.data.frame(reactome)[1:10,])
)
write.csv(combined, 'top_enriched_pathways.csv', row.names = FALSE)
```


## Complete Workflow Script


```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# Configuration
de_file <- 'deseq2_results.csv'
output_dir <- 'pathway_analysis'
dir.create(output_dir, showWarnings = FALSE)

# Load and prepare data
res <- read.csv(de_file, row.names = 1)
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
cat('Significant genes:', length(sig_genes), '\n')

# Convert IDs
sig_entrez <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
cat('Converted to Entrez:', nrow(sig_entrez), '\n')

# Ranked list for GSEA (Wald statistic preferred over LFC)
ranked <- res$stat
names(ranked) <- rownames(res)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)
ranked_entrez <- bitr(names(ranked), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
ranked_list <- ranked[ranked_entrez$SYMBOL]
names(ranked_list) <- ranked_entrez$ENTREZID

# Background genes (all tested, not full genome)
bg_entrez <- bitr(rownames(res[!is.na(res$pvalue), ]), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# GO enrichment with background
go_bp <- enrichGO(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = TRUE)
go_bp_simple <- simplify(go_bp, cutoff = 0.7)

# KEGG
kegg <- enrichKEGG(sig_entrez$ENTREZID, organism = 'hsa')
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# Reactome
reactome <- enrichPathway(sig_entrez$ENTREZID, organism = 'human', readable = TRUE)

# GSEA
gsea_go <- gseGO(ranked_list, OrgDb = org.Hs.eg.db, ont = 'BP', verbose = FALSE)

# Plots
pdf(file.path(output_dir, 'enrichment_plots.pdf'), width = 10, height = 8)
print(dotplot(go_bp_simple, showCategory = 20) + ggtitle('GO Biological Process'))
print(barplot(kegg, showCategory = 15) + ggtitle('KEGG Pathways'))
if (nrow(as.data.frame(reactome)) > 0) {
    print(dotplot(reactome, showCategory = 15) + ggtitle('Reactome Pathways'))
}
dev.off()

# Export
write.csv(as.data.frame(go_bp_simple), file.path(output_dir, 'go_bp.csv'), row.names = FALSE)
write.csv(as.data.frame(kegg), file.path(output_dir, 'kegg.csv'), row.names = FALSE)
write.csv(as.data.frame(reactome), file.path(output_dir, 'reactome.csv'), row.names = FALSE)

cat('\nResults saved to:', output_dir, '\n')
cat('GO BP terms:', nrow(as.data.frame(go_bp_simple)), '\n')
cat('KEGG pathways:', nrow(as.data.frame(kegg)), '\n')
cat('Reactome pathways:', nrow(as.data.frame(reactome)), '\n')
```


## Prokaryotic Organisms


For bacteria/archaea, standard org.db annotation packages are unavailable. Use KEGG directly with strain-specific organism codes:

```r
# Find organism code
search_kegg_organism('Pseudomonas aeruginosa', by = 'scientific_name')

# KEGG ORA with bacterial organism code (e.g., 'pae' for P. aeruginosa PAO1)
kegg_bac <- enrichKEGG(gene = sig_gene_ids, organism = 'pae', keyType = 'kegg',
                       pvalueCutoff = 0.05)

# For organisms without KEGG annotation, use KEGG Orthology
# Map genes to KO IDs via eggNOG-mapper or KOALA, then:
kegg_ko <- enrichKEGG(gene = ko_ids, organism = 'ko', keyType = 'kegg')
```

GO enrichment for prokaryotes: use `enricher()` with custom GO-to-gene mapping from eggNOG-mapper or InterProScan output, rather than org.db packages.


## Multi-Condition Enrichment Comparison


When comparing enrichment across conditions (e.g., treatment A vs B vs C):

```r
# compareCluster: run ORA across multiple gene lists
gene_clusters <- list(
    ConditionA = sig_genes_A,
    ConditionB = sig_genes_B,
    ConditionC = sig_genes_C
)
cc <- compareCluster(gene_clusters, fun = 'enrichKEGG', organism = 'hsa')
dotplot(cc, showCategory = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Do not compare raw -log10(p-values) across conditions — they scale with sample size. Compare NES (normalized enrichment scores) for GSEA, or use compareCluster for ORA.

