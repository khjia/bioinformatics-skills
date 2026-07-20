---
name: bio-rna-quantification
description: "Complete RNA-seq quantification pipeline: Salmon/kallisto pseudo-alignment, featureCounts, tximport import, and count matrix QC. Covers index building, sample processing, strand-specific libraries, batch processing, and count matrix validation."
---


## Salmon Workflow


### Build Index

```bash
# Download transcriptome FASTA
# Ensembl: Homo_sapiens.GRCh38.cdna.all.fa.gz

# Basic index (fast, less accurate)
salmon index -t transcripts.fa -i salmon_index

# Decoy-aware index (recommended for accuracy)
# First, create decoys from genome
grep "^>" genome.fa | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
cat transcripts.fa genome.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i salmon_index -p 8
```

### Quantify Samples

```bash
# Paired-end reads
salmon quant -i salmon_index -l A \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -o sample_quant -p 8

# Single-end reads
salmon quant -i salmon_index -l A \
    -r sample.fastq.gz \
    -o sample_quant -p 8
```

**Key flags:**
- `-l A` - Automatically detect library type
- `-p` - Number of threads
- `--validateMappings` - More accurate (default in recent versions)
- `--gcBias` - Correct for GC bias
- `--seqBias` - Correct for sequence-specific bias

### Library Types

| Code | Description |
|------|-------------|
| `A` | Automatic detection (recommended) |
| `ISR` | Inward, stranded, read 1 from reverse |
| `ISF` | Inward, stranded, read 1 from forward |
| `IU` | Inward, unstranded |

### Batch Processing

```bash
for sample in sample1 sample2 sample3; do
    salmon quant -i salmon_index -l A \
        -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz \
        -o ${sample}_quant -p 8
done
```

### Output Files

```
sample_quant/
├── quant.sf           # Main quantification file
├── aux_info/          # Auxiliary information
├── cmd_info.json      # Command used
├── lib_format_counts.json  # Library format detection
└── logs/              # Log files
```

**quant.sf format:**
```
Name                    Length  EffectiveLength TPM         NumReads
ENST00000456328.2       1657    1477.000        0.000000    0.000
ENST00000450305.2       632     452.000         12.345678   156.789
```


## kallisto Workflow


### Build Index

```bash
kallisto index -i kallisto_index transcripts.fa
```

### Quantify Samples

```bash
# Paired-end
kallisto quant -i kallisto_index -o sample_quant \
    sample_R1.fastq.gz sample_R2.fastq.gz

# Single-end (must specify fragment length)
kallisto quant -i kallisto_index -o sample_quant \
    --single -l 200 -s 20 sample.fastq.gz

# With bootstraps (for sleuth)
kallisto quant -i kallisto_index -o sample_quant -b 100 \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

**Key flags:**
- `-b` - Number of bootstrap samples
- `-t` - Number of threads
- `--single` - Single-end mode
- `-l` - Estimated fragment length (single-end)
- `-s` - Fragment length standard deviation

### Output Files

```
sample_quant/
├── abundance.tsv      # Main quantification (text)
├── abundance.h5       # HDF5 format (for sleuth)
└── run_info.json      # Run information
```

**abundance.tsv format:**
```
target_id               length  eff_length  est_counts  tpm
ENST00000456328.2       1657    1477.00     0.00        0.000000
ENST00000450305.2       632     452.00      156.79      12.345678
```


## Salmon vs kallisto


| Feature | Salmon | kallisto |
|---------|--------|----------|
| Speed | Fast | Fastest |
| Accuracy | Higher | Good |
| GC bias correction | Yes | No |
| Decoy sequences | Yes | No |
| Memory usage | Moderate | Low |

**Recommendation:** Use Salmon for production, kallisto for quick exploratory analysis.


## Combining Results


```bash
# Salmon: use tximport in R
# kallisto: use tximport or sleuth

# Quick Python combination
python << 'EOF'
import pandas as pd
from pathlib import Path

samples = ['sample1', 'sample2', 'sample3']
tpm_data = {}
counts_data = {}

for sample in samples:
    quant_file = Path(f'{sample}_quant/quant.sf')  # Salmon
    # quant_file = Path(f'{sample}_quant/abundance.tsv')  # kallisto
    df = pd.read_csv(quant_file, sep='\t', index_col=0)
    tpm_data[sample] = df['TPM']
    counts_data[sample] = df['NumReads']  # or est_counts for kallisto

tpm_matrix = pd.DataFrame(tpm_data)
counts_matrix = pd.DataFrame(counts_data)
tpm_matrix.to_csv('tpm_matrix.csv')
counts_matrix.to_csv('counts_matrix.csv')
EOF
```


## Quality Checks


```bash
# Check mapping rate from Salmon logs
grep "Mapping rate" sample_quant/logs/salmon_quant.log

# Check library type detection
cat sample_quant/lib_format_counts.json
```

**Good metrics:**
- Mapping rate > 70%
- Consistent library type across samples


## Common Issues


**Low mapping rate:**
- Wrong transcriptome version
- Contamination in samples
- Wrong library type

**Inconsistent library types:**
- Mixed library preparations
- Sample swap


## Related Skills


- read-qc/fastp-workflow - Upstream preprocessing
- rna-quantification/tximport-workflow - Import results to R
- rna-quantification/count-matrix-qc - QC of quantification
- differential-expression/deseq2-basics - Downstream analysis


## Basic Usage


```bash
# Single sample
featureCounts -a annotation.gtf -o counts.txt aligned.bam

# Multiple samples (recommended - single matrix output)
featureCounts -a annotation.gtf -o counts.txt sample1.bam sample2.bam sample3.bam

# All BAMs in directory
featureCounts -a annotation.gtf -o counts.txt *.bam
```


## Paired-End Data


```bash
# Count fragments, not reads (required for paired-end)
featureCounts -p --countReadPairs -a annotation.gtf -o counts.txt *.bam

# Check proper pairs only
featureCounts -p --countReadPairs -B -C -a annotation.gtf -o counts.txt *.bam
```

**Flags:**
- `-p` - Input is paired-end
- `--countReadPairs` - Count fragments instead of reads
- `-B` - Only count properly paired reads
- `-C` - Don't count chimeric fragments


## Strand-Specific Libraries


```bash
# Unstranded (default)
featureCounts -s 0 -a annotation.gtf -o counts.txt *.bam

# Forward stranded (e.g., dUTP, NSR)
featureCounts -s 1 -a annotation.gtf -o counts.txt *.bam

# Reverse stranded (e.g., Illumina TruSeq, most common)
featureCounts -s 2 -a annotation.gtf -o counts.txt *.bam
```

**Determining strandedness:** Use `infer_experiment.py` from RSeQC or check library prep protocol.


## Feature Types


```bash
# Count at gene level (default)
featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt *.bam

# Count at transcript level
featureCounts -t exon -g transcript_id -a annotation.gtf -o counts.txt *.bam

# Count CDS only
featureCounts -t CDS -g gene_id -a annotation.gtf -o counts.txt *.bam
```

**Flags:**
- `-t` - Feature type in GTF (default: exon)
- `-g` - Meta-feature attribute (default: gene_id)


## Multi-Mapping Reads


```bash
# Discard multi-mappers (default, recommended for DE)
featureCounts -a annotation.gtf -o counts.txt *.bam

# Count multi-mappers (fractional)
featureCounts -M --fraction -a annotation.gtf -o counts.txt *.bam

# Count multi-mappers (full count to each location)
featureCounts -M -a annotation.gtf -o counts.txt *.bam
```


## Overlapping Features


```bash
# Discard reads overlapping multiple features (default)
featureCounts -a annotation.gtf -o counts.txt *.bam

# Count reads overlapping multiple features
featureCounts -O -a annotation.gtf -o counts.txt *.bam

# Fractional count for overlaps
featureCounts -O --fraction -a annotation.gtf -o counts.txt *.bam
```


## Performance Options


```bash
# Use multiple threads
featureCounts -T 8 -a annotation.gtf -o counts.txt *.bam

# Use less memory (slower)
featureCounts --largeBAM -a annotation.gtf -o counts.txt *.bam
```


## Output Files


featureCounts produces two files:

1. **counts.txt** - Main count matrix
```
Geneid  Chr     Start   End     Strand  Length  sample1.bam  sample2.bam
GENE1   chr1    100     500     +       400     1523         1891
GENE2   chr1    1000    2000    -       1000    892          756
```

2. **counts.txt.summary** - Assignment statistics
```
Status                  sample1.bam  sample2.bam
Assigned                1523456      1678234
Unassigned_Unmapped     12345        11234
Unassigned_NoFeatures   234567       245678
```


## Extract Count Matrix


```bash
# Remove first 6 columns (metadata) to get just counts
cut -f1,7- counts.txt | tail -n +2 > count_matrix.txt
```


## Python Processing


```python
import pandas as pd

counts = pd.read_csv('counts.txt', sep='\t', comment='#')
count_matrix = counts.set_index('Geneid').iloc[:, 5:]  # Skip metadata columns
count_matrix.columns = [c.replace('.bam', '') for c in count_matrix.columns]
count_matrix.to_csv('count_matrix.csv')
```


## R Processing


```r
counts <- read.table('counts.txt', header=TRUE, row.names=1, skip=1)
count_matrix <- counts[, 6:ncol(counts)]  # Skip metadata columns
colnames(count_matrix) <- gsub('.bam', '', colnames(count_matrix))
```


## Basic tximport


**Goal:** Import transcript-level quantifications from Salmon or kallisto into R as gene-level counts with proper length-offset correction for DESeq2 or edgeR.

**Approach:** Create a transcript-to-gene mapping from a GTF or biomaRt, then run tximport on the quantification files to produce a gene-level count matrix with length-scaled TPM offsets.

```r
library(tximport)

# Define sample files
files <- c(
    sample1 = 'sample1_quant/quant.sf',
    sample2 = 'sample2_quant/quant.sf',
    sample3 = 'sample3_quant/quant.sf'
)

# Load transcript-to-gene mapping
tx2gene <- read.csv('tx2gene.csv')  # columns: TXNAME, GENEID

# Import at gene level
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
```


## Creating tx2gene Mapping


### From GTF (using GenomicFeatures)

```r
library(GenomicFeatures)

txdb <- makeTxDbFromGFF('annotation.gtf')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, k, 'GENEID', 'TXNAME')
```

### From Ensembl (using biomaRt)

```r
library(biomaRt)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
tx2gene <- getBM(
    attributes = c('ensembl_transcript_id_version', 'ensembl_gene_id_version'),
    mart = mart
)
colnames(tx2gene) <- c('TXNAME', 'GENEID')
```

### From Salmon quant.sf

```r
quant <- read.table('sample1_quant/quant.sf', header = TRUE)
tx2gene <- data.frame(
    TXNAME = quant$Name,
    GENEID = gsub('\\..*', '', quant$Name)  # Remove version
)
```


## Import Types


### Gene-Level Summarization (Default)

```r
# Summarize transcripts to gene level
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
# Returns: counts, abundance (TPM), length at gene level
```

### Transcript-Level (No Summarization)

```r
# Keep transcript-level estimates
txi <- tximport(files, type = 'salmon', txOut = TRUE)
# Returns: counts, abundance, length at transcript level
```

### Scaled TPM (for visualization)

```r
# Gene-level TPM
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene,
                countsFromAbundance = 'scaledTPM')
```


## Source-Specific Import


### Salmon

```r
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
```

### kallisto

```r
txi <- tximport(files, type = 'kallisto', tx2gene = tx2gene)
```

### RSEM

```r
txi <- tximport(files, type = 'rsem', tx2gene = tx2gene)
```

### StringTie

```r
txi <- tximport(files, type = 'stringtie', tx2gene = tx2gene)
```


## Using with DESeq2


```r
library(DESeq2)

# Create sample metadata
coldata <- data.frame(
    condition = factor(c('control', 'control', 'treated', 'treated')),
    row.names = names(files)
)

# Create DESeqDataSet from tximport
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
```


## Using with edgeR


```r
library(edgeR)

# Create DGEList with offset
cts <- txi$counts
normMat <- txi$length
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts / normMat)) + log(colSums(cts / normMat))

y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)

# Continue with edgeR analysis
y <- estimateDisp(y, design)
```


## tximeta: Metadata-Aware Import


tximeta automatically attaches transcript and gene information from the original annotation.

```r
library(tximeta)

# First time: link transcriptome to annotation
makeLinkedTxome(
    indexDir = 'salmon_index',
    source = 'Ensembl',
    organism = 'Homo sapiens',
    release = '110',
    genome = 'GRCh38',
    fasta = 'transcripts.fa',
    gtf = 'annotation.gtf'
)

# Import with full metadata
coldata <- data.frame(
    files = files,
    names = names(files),
    condition = c('control', 'control', 'treated', 'treated')
)

se <- tximeta(coldata)

# Summarize to gene level
gse <- summarizeToGene(se)

# Convert to DESeqDataSet
dds <- DESeqDataSet(gse, design = ~ condition)
```


## tximport Output Structure


```r
names(txi)
# [1] "abundance"           "counts"              "length"
# [4] "countsFromAbundance"

# abundance: TPM values (genes x samples)
# counts: estimated counts (genes x samples)
# length: effective gene lengths (genes x samples)
```


## Handling Version Numbers


```r
# Remove version from transcript IDs
tx2gene$TXNAME <- gsub('\\.\\d+$', '', tx2gene$TXNAME)

# Or ignore version during import
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene,
                ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
```


## Load and Inspect Counts


**Goal:** Assess count matrix quality before differential expression by detecting outliers, batch effects, and sample relationship problems.

**Approach:** Load counts into DESeq2 or pandas, compute per-sample library size statistics, apply variance-stabilizing transformation, then run PCA and sample-sample correlation to identify outliers and batch structure.

### R

```r
library(DESeq2)

# From tximport
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# From count matrix
counts <- read.csv('count_matrix.csv', row.names = 1)
coldata <- data.frame(condition = factor(c('ctrl', 'ctrl', 'treat', 'treat')),
                      row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~ condition)
```

### Python

```python
import pandas as pd
import numpy as np

counts = pd.read_csv('count_matrix.csv', index_col=0)
metadata = pd.read_csv('sample_info.csv', index_col=0)
```


## Basic Statistics


### R

```r
# Total counts per sample
colSums(counts(dds))

# Genes detected per sample
colSums(counts(dds) > 0)

# Counts summary
summary(colSums(counts(dds)))
```

### Python

```python
total_counts = counts.sum()
genes_detected = (counts > 0).sum()

print('Total counts per sample:')
print(total_counts)
print('\nGenes detected:')
print(genes_detected)
```


## Filter Low-Count Genes


### R

```r
# Remove genes with low counts across samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# More stringent: at least N samples with count >= M
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
```

### Python

```python
min_counts = 10
min_samples = 3

gene_filter = (counts >= min_counts).sum(axis=1) >= min_samples
counts_filtered = counts[gene_filter]
```


## Normalize for Visualization


### R (DESeq2 VST)

```r
# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Or regularized log (slower, better for small n)
rld <- rlog(dds, blind = TRUE)

# Get transformed values
vst_matrix <- assay(vsd)
```

### Python (log2 CPM)

```python
from sklearn.preprocessing import StandardScaler

cpm = counts * 1e6 / counts.sum()
log_cpm = np.log2(cpm + 1)
```


## Sample Correlation


### R

```r
library(pheatmap)

# Sample correlation heatmap
sample_cor <- cor(assay(vsd))
pheatmap(sample_cor, annotation_col = coldata)

# Sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dist), annotation_col = coldata)
```

### Python

```python
import seaborn as sns
import matplotlib.pyplot as plt

sample_cor = log_cpm.corr()
sns.clustermap(sample_cor, annot=True, cmap='RdBu_r', center=0.9,
               vmin=0.8, vmax=1.0)
plt.savefig('sample_correlation.png')
```


## PCA Analysis


### R

```r
# PCA plot
plotPCA(vsd, intgroup = 'condition')

# Custom PCA
pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     condition = coldata$condition)

library(ggplot2)
ggplot(pca_df, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text(aes(label = rownames(pca_df)), vjust = -0.5)
```

### Python

```python
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_cpm.T)

plt.figure(figsize=(8, 6))
for condition in metadata['condition'].unique():
    mask = metadata['condition'] == condition
    plt.scatter(pca_result[mask, 0], pca_result[mask, 1], label=condition)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
plt.legend()
plt.savefig('pca_plot.png')
```


## Detect Outliers


### R

```r
# Cook's distance (after DESeq)
dds <- DESeq(dds)
W <- results(dds)$cooksd
boxplot(W, main = "Cook's Distance")

# Identify outlier samples from PCA
pca <- prcomp(t(assay(vsd)))
outliers <- abs(scale(pca$x[,1])) > 3 | abs(scale(pca$x[,2])) > 3
```

### Python

```python
from scipy import stats

z_scores = stats.zscore(pca_result, axis=0)
outliers = (np.abs(z_scores) > 3).any(axis=1)
print('Potential outliers:', counts.columns[outliers].tolist())
```


## Check for Batch Effects


### R

```r
# Color PCA by batch
plotPCA(vsd, intgroup = c('condition', 'batch'))

# Test for batch effect
design(dds) <- ~ batch + condition
dds <- DESeq(dds)
```

### Python

```python
# Color by batch in PCA
for batch in metadata['batch'].unique():
    mask = metadata['batch'] == batch
    plt.scatter(pca_result[mask, 0], pca_result[mask, 1],
                marker=['o', 's', '^'][list(metadata['batch'].unique()).index(batch)],
                label=f'Batch {batch}')
```


## Library Complexity


### R

```r
# Genes detected vs library size
plot(colSums(counts(dds)), colSums(counts(dds) > 0),
     xlab = 'Library Size', ylab = 'Genes Detected')

# Saturation check
```

### Python

```python
plt.scatter(counts.sum(), (counts > 0).sum())
plt.xlabel('Total Counts')
plt.ylabel('Genes Detected')
plt.savefig('library_complexity.png')
```


## Gene-Level QC


### R

```r
# Most variable genes
rv <- rowVars(assay(vsd))
top_var <- order(rv, decreasing = TRUE)[1:500]

# Expression distribution
boxplot(log2(counts(dds) + 1), las = 2)
```

### Python

```python
gene_var = log_cpm.var(axis=1).sort_values(ascending=False)
top_var_genes = gene_var.head(500).index

counts[top_var_genes].boxplot(figsize=(12, 6))
plt.xticks(rotation=45)
plt.savefig('gene_expression_dist.png')
```


## Summary Report


```r
# Quick summary
cat('Samples:', ncol(dds), '\n')
cat('Genes before filter:', nrow(counts), '\n')
cat('Genes after filter:', nrow(dds), '\n')
cat('Median library size:', median(colSums(counts(dds))), '\n')
cat('Median genes detected:', median(colSums(counts(dds) > 0)), '\n')
```

