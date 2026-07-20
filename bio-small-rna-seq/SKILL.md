---
name: bio-small-rna-seq
description: "Complete small RNA-seq pipeline: adapter trimming, miRNA discovery (miRDeep2), quantification (miRge3), differential miRNA expression, and target prediction. Covers 18-30nt sRNA workflows."
---


## Adapter Trimming with Cutadapt


**Goal:** Remove 3' adapter sequences and size-select reads in the small RNA range.

**Approach:** Run cutadapt with the kit-specific adapter, minimum/maximum length filters, and discard reads without adapter.

Small RNA libraries have specific 3' adapters that must be removed:

```bash
# Standard Illumina TruSeq small RNA adapter
cutadapt \
    -a TGGAATTCTCGGGTGCCAAGG \
    -m 18 \
    -M 30 \
    --discard-untrimmed \
    -o trimmed.fastq.gz \
    input.fastq.gz

# -a: 3' adapter sequence
# -m 18: Minimum length (miRNAs are 18-25 nt)
# -M 30: Maximum length (exclude longer fragments)
# --discard-untrimmed: Remove reads without adapter (likely not small RNA)
```


## Common Small RNA Adapters


| Kit | 3' Adapter Sequence |
|-----|---------------------|
| Illumina TruSeq | TGGAATTCTCGGGTGCCAAGG |
| NEBNext | AGATCGGAAGAGCACACGTCT |
| QIAseq | AACTGTAGGCACCATCAAT |
| Lexogen | TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC |


## Size Selection


```bash
# Filter by length after trimming
cutadapt \
    -a TGGAATTCTCGGGTGCCAAGG \
    -m 18 -M 26 \
    -o mirna_length.fastq.gz \
    input.fastq.gz

# miRNA: 18-26 nt (typically 21-23 nt)
# piRNA: 26-32 nt
# snoRNA: variable, typically longer
```


## Quality Trimming


```bash
# Trim low-quality bases from 3' end before adapter removal
cutadapt \
    -q 20 \
    -a TGGAATTCTCGGGTGCCAAGG \
    -m 18 \
    -o trimmed.fastq.gz \
    input.fastq.gz
```


## Using fastp for Small RNA


```bash
# fastp with small RNA settings
fastp \
    --in1 input.fastq.gz \
    --out1 trimmed.fastq.gz \
    --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
    --length_required 18 \
    --length_limit 30 \
    --html report.html

# Note: fastp auto-detects adapters but specifying is more reliable
```


## Collapse Identical Reads


For small RNAs, collapsing identical sequences reduces computation:

```bash
# Using seqkit
seqkit rmdup -s trimmed.fastq.gz -o collapsed.fasta

# Using fastx_toolkit (legacy)
fastx_collapser -i trimmed.fastq -o collapsed.fasta
```


## Python Preprocessing


```python
import gzip
from collections import Counter

def collapse_reads(fastq_path):
    '''Collapse identical sequences and count occurrences'''
    counts = Counter()

    with gzip.open(fastq_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual

            # Only keep reads in miRNA size range
            if 18 <= len(seq) <= 26:
                counts[seq] += 1

    return counts

# Write collapsed FASTA
def write_collapsed_fasta(counts, output_path):
    with open(output_path, 'w') as f:
        for i, (seq, count) in enumerate(counts.most_common()):
            f.write(f'>seq_{i}_x{count}\n{seq}\n')
```


## QC Metrics for Small RNA


Key metrics to check:
- Read length distribution (should peak at 21-23 nt for miRNA)
- Adapter content (high if library is good)
- Percentage of reads in target size range

```python
import matplotlib.pyplot as plt
from collections import Counter

def plot_length_distribution(fastq_path):
    lengths = Counter()
    with gzip.open(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                lengths[len(line.strip())] += 1

    plt.bar(lengths.keys(), lengths.values())
    plt.xlabel('Read Length')
    plt.ylabel('Count')
    plt.title('Small RNA Length Distribution')
    plt.savefig('length_dist.png')
```


## Related Skills


- mirdeep2-analysis - Novel miRNA discovery
- mirge3-analysis - Fast miRNA quantification
- read-qc/adapter-trimming - General adapter trimming


## Workflow Overview


```
Collapsed reads (FASTA)
    |
    v
mapper.pl ---------> Align to genome, create ARF file
    |
    v
miRDeep2.pl -------> Predict novel miRNAs, quantify known
    |
    v
quantifier.pl -----> Quantify known miRNAs only (optional)
```


## Step 1: Prepare Genome Index


**Goal:** Build a bowtie index from the reference genome for miRDeep2 read mapping.

**Approach:** Run bowtie-build on the genome FASTA to create the index files required by mapper.pl.

```bash
# Build bowtie index for miRDeep2 mapper
bowtie-build genome.fa genome_index
```


## Step 2: Map Reads with mapper.pl


**Goal:** Collapse identical reads and align them to the reference genome.

**Approach:** Use mapper.pl to clip adapters, filter by length, collapse duplicates, and map with bowtie to produce ARF alignment files.

```bash
# Collapse reads and map to genome
mapper.pl reads.fastq \
    -e \
    -h \
    -i \
    -j \
    -k TGGAATTCTCGGGTGCCAAGG \
    -l 18 \
    -m \
    -p genome_index \
    -s reads_collapsed.fa \
    -t reads_vs_genome.arf \
    -v

# Key options:
# -e: Input is FASTQ
# -h: Parse Illumina headers
# -k: Clip 3' adapter
# -l 18: Discard reads < 18 nt
# -m: Collapse reads
# -p: Bowtie index prefix
# -s: Output collapsed FASTA
# -t: Output ARF alignment file
```


## Step 3: Run miRDeep2 Prediction


**Goal:** Predict novel miRNAs and quantify known miRNAs from aligned small RNA reads.

**Approach:** Run miRDeep2.pl with collapsed reads, genome, alignments, and miRBase references to score candidate miRNA loci.

```bash
# Predict novel miRNAs
miRDeep2.pl \
    reads_collapsed.fa \
    genome.fa \
    reads_vs_genome.arf \
    mature_ref.fa \
    mature_other.fa \
    hairpin_ref.fa \
    -t Human \
    2> report.log

# Arguments:
# 1. Collapsed reads FASTA
# 2. Genome FASTA
# 3. Alignment ARF file
# 4. Known mature miRNAs (same species)
# 5. Known mature miRNAs (other species, for conservation)
# 6. Known hairpin precursors
# -t: Species for miRBase lookup
```


## Prepare miRBase References


**Goal:** Download and extract species-specific miRNA references from miRBase.

**Approach:** Fetch mature and hairpin FASTA files from miRBase, then grep species-specific entries by prefix.

```bash
# Download from miRBase
wget https://www.mirbase.org/download/mature.fa
wget https://www.mirbase.org/download/hairpin.fa

# Extract species-specific sequences
grep -A1 ">hsa-" mature.fa > mature_human.fa
grep -A1 ">hsa-" hairpin.fa > hairpin_human.fa
```


## Step 4: Quantify Known miRNAs Only


**Goal:** Quantify expression of known miRNAs without running novel discovery.

**Approach:** Run quantifier.pl with hairpin and mature references against collapsed reads for fast quantification.

```bash
# If not doing novel discovery
quantifier.pl \
    -p hairpin_human.fa \
    -m mature_human.fa \
    -r reads_collapsed.fa \
    -t hsa

# Output: miRNAs_expressed_all_samples.csv
```


## Output Files


| File | Description |
|------|-------------|
| result_*.html | Interactive results report |
| result_*.csv | Predicted novel miRNAs with scores |
| miRNAs_expressed_all_samples*.csv | Expression quantification |
| pdfs_*.pdf | Secondary structure plots |


## Interpret miRDeep2 Scores


```
Score interpretation:
>10: High confidence novel miRNA
5-10: Medium confidence
1-5: Low confidence, needs validation
<1: Likely false positive

Key metrics:
- miRDeep2 score: Overall confidence
- Total read count: Expression level
- Mature/star ratio: Strand bias (expect asymmetry)
- Randfold p-value: Structural stability
```


## Parse Results in Python


**Goal:** Load miRDeep2 prediction and quantification results into pandas DataFrames.

**Approach:** Parse tab-delimited output files and filter novel miRNA predictions by confidence score threshold.

```python
import pandas as pd

def parse_mirdeep2_results(csv_path):
    '''Parse miRDeep2 novel miRNA predictions'''
    df = pd.read_csv(csv_path, sep='\t', skiprows=1)

    # Filter high-confidence predictions
    # Score > 10 indicates high confidence novel miRNA
    high_conf = df[df['miRDeep2 score'] > 10]

    return high_conf

# Parse quantification results
def parse_quantifier_output(csv_path):
    '''Parse quantifier.pl expression matrix'''
    df = pd.read_csv(csv_path, sep='\t')
    return df
```


## Basic Quantification


**Goal:** Quantify known miRNA expression from small RNA-seq FASTQ files.

**Approach:** Run miRge3 annotation pipeline with adapter trimming, organism-specific libraries, and multi-sample input.

```bash
# Run miRge3 on FASTQ files
miRge3.0 annotate \
    -s sample1.fastq.gz,sample2.fastq.gz \
    -lib miRge3_libs \
    -on human \
    -db mirbase \
    -o output_dir \
    -a TGGAATTCTCGGGTGCCAAGG \
    --threads 8

# Key options:
# -s: Input FASTQ files (comma-separated)
# -lib: Path to miRge3 library
# -on: Organism name
# -db: Database (mirbase or mirgenedb)
# -a: 3' adapter sequence
```


## Install miRge3 Libraries


**Goal:** Download organism-specific reference libraries required for miRge3 annotation.

**Approach:** Use miRge3 built-in download command to fetch pre-built bowtie indices and annotations.

```bash
# Download pre-built libraries
miRge3.0 --download-library human mirbase

# Libraries include:
# - Bowtie indices for miRNAs, tRNAs, rRNAs
# - miRBase or MirGeneDB annotations
# - A-to-I editing sites
```


## IsomiR Detection


**Goal:** Identify and quantify isomiR variants including 5'/3' additions, deletions, and internal modifications.

**Approach:** Enable miRge3 isomiR mode to classify reads by their deviation from canonical miRNA sequences.

```bash
# Enable isomiR analysis
miRge3.0 annotate \
    -s sample.fastq.gz \
    -lib miRge3_libs \
    -on human \
    -db mirbase \
    --isomir \
    -o output_dir

# IsomiRs include:
# - 5' variants (templated and non-templated)
# - 3' variants (templated and non-templated)
# - Internal modifications
```


## A-to-I RNA Editing


**Goal:** Detect adenosine-to-inosine RNA editing events in miRNA sequences.

**Approach:** Enable miRge3 A-to-I detection mode which identifies editing sites and calculates editing frequencies.

```bash
# Detect A-to-I editing
miRge3.0 annotate \
    -s sample.fastq.gz \
    -lib miRge3_libs \
    -on human \
    -db mirbase \
    --AtoI \
    -o output_dir

# Outputs editing sites and frequencies
```


## Python API


**Goal:** Run miRge3 quantification programmatically from Python.

**Approach:** Call the miRge3 annotate function directly with configuration parameters instead of CLI invocation.

```python
from mirge3.annotate import annotate

# Run programmatically
annotate(
    samples=['sample1.fastq.gz', 'sample2.fastq.gz'],
    lib_path='miRge3_libs',
    organism='human',
    database='mirbase',
    adapter='TGGAATTCTCGGGTGCCAAGG',
    output_dir='results',
    threads=8
)
```


## Parse miRge3 Output


**Goal:** Load miRge3 count matrices and isomiR tables into pandas for downstream analysis.

**Approach:** Read CSV output files and apply minimum count filtering to remove lowly-expressed miRNAs.

```python
import pandas as pd

def load_mirge3_counts(output_dir):
    '''Load miRge3 count matrix'''
    counts = pd.read_csv(f'{output_dir}/miR.Counts.csv', index_col=0)
    return counts

def load_isomirs(output_dir):
    '''Load isomiR-level counts'''
    isomirs = pd.read_csv(f'{output_dir}/isomiR.Counts.csv', index_col=0)
    return isomirs

# Filter low-expressed miRNAs
def filter_low_counts(counts, min_total=10):
    '''Keep miRNAs with total count >= threshold'''
    return counts[counts.sum(axis=1) >= min_total]
```


## Compare Multiple Samples


**Goal:** Normalize and transform miRNA counts for cross-sample comparison.

**Approach:** Apply RPM normalization to account for library size, then log2-transform for variance stabilization.

```python
def normalize_rpm(counts):
    '''Normalize to reads per million'''
    total_per_sample = counts.sum(axis=0)
    rpm = counts / total_per_sample * 1e6
    return rpm

def log_transform(rpm, pseudocount=1):
    '''Log2 transform with pseudocount'''
    import numpy as np
    return np.log2(rpm + pseudocount)
```


## IsomiR Analysis


**Goal:** Summarize isomiR diversity metrics per canonical miRNA.

**Approach:** Group isomiR-level counts by parent miRNA and compute total reads, variant count, and dominant isoform.

```python
def summarize_isomirs(isomir_counts):
    '''Summarize isomiR diversity per miRNA'''
    # Group by canonical miRNA
    isomir_counts['miRNA'] = isomir_counts.index.str.extract(r'(hsa-\w+-\d+[a-z]*)')[0]

    summary = isomir_counts.groupby('miRNA').agg({
        'count': ['sum', 'count', lambda x: x.idxmax()]
    })
    summary.columns = ['total_reads', 'n_isomirs', 'dominant_isomir']
    return summary
```


## Load miRNA Count Data


```r
library(DESeq2)

# Load miRge3 or miRDeep2 counts
counts <- read.csv('miR.Counts.csv', row.names = 1)

# Create sample metadata
coldata <- data.frame(
    sample = colnames(counts),
    condition = factor(c('control', 'control', 'treated', 'treated')),
    row.names = colnames(counts)
)
```


## DESeq2 Analysis


**Goal:** Identify miRNAs with significant expression changes between experimental conditions, accounting for small RNA-specific normalization.

**Approach:** Create a DESeqDataSet from miRNA counts, filter low-expressed miRNAs using a lower threshold than mRNA (10 reads total), run the DESeq2 pipeline, and extract results with BH-corrected p-values.

```r
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),  # DESeq2 requires integers
    colData = coldata,
    design = ~ condition
)

# Filter low-expressed miRNAs
# miRNAs typically have fewer total counts than mRNAs
# Keep miRNAs with at least 10 reads across samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c('condition', 'treated', 'control'))
res <- res[order(res$padj), ]
```


## Apply Shrinkage for Effect Sizes


```r
# apeglm shrinkage for more accurate log2 fold changes
# Particularly important for low-count miRNAs
library(apeglm)

res_shrunk <- lfcShrink(
    dds,
    coef = 'condition_treated_vs_control',
    type = 'apeglm'
)
```


## Filter Significant miRNAs


```r
# Standard thresholds for miRNA DE
# padj < 0.05: FDR-corrected significance
# |log2FC| > 1: 2-fold change minimum
sig <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
sig <- sig[order(sig$padj), ]

# Separate up and down-regulated
up <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)

cat('Upregulated:', nrow(up), '\n')
cat('Downregulated:', nrow(down), '\n')
```


## edgeR Alternative


```r
library(edgeR)

# Create DGEList
dge <- DGEList(counts = counts, group = coldata$condition)

# Filter low expression
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ condition, data = coldata)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit model and test
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

# Get results
res_edger <- topTags(qlf, n = Inf)$table
```


## Visualization


```r
library(ggplot2)
library(EnhancedVolcano)

# Volcano plot
EnhancedVolcano(
    res_shrunk,
    lab = rownames(res_shrunk),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    title = 'Differential miRNA Expression'
)

# MA plot
plotMA(res_shrunk, ylim = c(-4, 4))
```


## Heatmap of DE miRNAs


```r
library(pheatmap)

# Get normalized counts
vsd <- vst(dds, blind = FALSE)

# Select significant miRNAs
sig_mirnas <- rownames(sig)
mat <- assay(vsd)[sig_mirnas, ]

# Z-score scale rows
mat_scaled <- t(scale(t(mat)))

pheatmap(
    mat_scaled,
    annotation_col = coldata['condition'],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = nrow(mat) < 50
)
```


## Export Results


```r
# Full results with normalized counts
res_df <- as.data.frame(res_shrunk)
res_df$miRNA <- rownames(res_df)
res_df$baseMean_norm <- rowMeans(counts(dds, normalized = TRUE)[rownames(res_df), ])

write.csv(res_df, 'DE_miRNAs_full.csv', row.names = FALSE)

# Significant only
write.csv(as.data.frame(sig), 'DE_miRNAs_significant.csv')
```


## miRanda Algorithm


**Goal:** Predict miRNA-mRNA target interactions using thermodynamic alignment scoring.

**Approach:** Run miRanda to align miRNA sequences against 3' UTR sequences with minimum score and energy thresholds.

```bash
# Run miRanda for target prediction
miranda miRNA.fa UTRs.fa \
    -sc 140 \
    -en -20 \
    -out predictions.txt

# Options:
# -sc 140: Minimum alignment score (default 140)
# -en -20: Maximum free energy threshold (kcal/mol)
# Higher score and lower energy = stronger prediction
```


## Parse miRanda Output


**Goal:** Extract miRNA-target interaction records from miRanda output into a structured DataFrame.

**Approach:** Parse the tab-delimited output lines starting with '>' to extract miRNA, target, score, energy, and position fields.

```python
import pandas as pd

def parse_miranda(output_file):
    '''Parse miRanda output file'''
    results = []
    with open(output_file) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    results.append({
                        'mirna': parts[0].lstrip('>'),
                        'target': parts[1],
                        'score': float(parts[2]),
                        'energy': float(parts[3]),
                        'position': parts[4]
                    })
    return pd.DataFrame(results)
```


## TargetScan Database Lookup


**Goal:** Retrieve conserved miRNA target predictions from the TargetScan database.

**Approach:** Query the downloadable TargetScan context++ score file by miRNA family name and rank by prediction score.

```python
import requests
import pandas as pd

def query_targetscan(mirna_family):
    '''Query TargetScan for predicted targets

    Note: TargetScan uses miRNA family names (e.g., miR-21-5p)
    '''
    # TargetScan provides downloadable files
    # For human: https://www.targetscan.org/vert_80/vert_80_data_download/
    targetscan_file = 'Predicted_Targets_Context_Scores.txt'

    df = pd.read_csv(targetscan_file, sep='\t')
    targets = df[df['miRNA family'] == mirna_family]
    return targets.sort_values('context++ score')
```


## miRDB Database Lookup


**Goal:** Retrieve machine-learning-based miRNA target predictions from miRDB.

**Approach:** Query the miRDB prediction file by miRNA ID and filter for high-confidence targets (score >= 80).

```python
def query_mirdb(mirna_id):
    '''Query miRDB for target predictions

    miRDB uses machine learning for target prediction
    Score > 80 indicates high confidence
    '''
    # Download from http://mirdb.org/download.html
    mirdb_file = 'miRDB_v6.0_prediction_result.txt'

    df = pd.read_csv(mirdb_file, sep='\t', header=None,
                     names=['mirna', 'target', 'score'])
    targets = df[df['mirna'] == mirna_id]
    return targets[targets['score'] >= 80].sort_values('score', ascending=False)
```


## Combine Multiple Databases


**Goal:** Identify high-confidence miRNA targets predicted by multiple independent algorithms.

**Approach:** Compute the intersection of predictions across miRanda, TargetScan, and miRDB, keeping targets found in at least N databases.

```python
def consensus_targets(mirna, min_databases=2):
    '''Find targets predicted by multiple databases

    More reliable targets are predicted by multiple algorithms
    '''
    miranda_targets = set(query_miranda_targets(mirna))
    targetscan_targets = set(query_targetscan_targets(mirna))
    mirdb_targets = set(query_mirdb_targets(mirna))

    # Count predictions per target
    all_targets = miranda_targets | targetscan_targets | mirdb_targets
    consensus = []

    for target in all_targets:
        count = sum([
            target in miranda_targets,
            target in targetscan_targets,
            target in mirdb_targets
        ])
        if count >= min_databases:
            consensus.append({
                'target': target,
                'n_databases': count,
                'miranda': target in miranda_targets,
                'targetscan': target in targetscan_targets,
                'mirdb': target in mirdb_targets
            })

    return pd.DataFrame(consensus).sort_values('n_databases', ascending=False)
```


## Python miRNA Target Prediction


**Goal:** Retrieve experimentally validated miRNA-target interactions from miRTarBase.

**Approach:** Load the miRTarBase Excel download and filter by miRNA name to get validated targets with experimental evidence types.

```python
# Using mirtarbase package for validated targets
def get_validated_targets(mirna):
    '''Get experimentally validated targets from miRTarBase'''
    # Download from https://mirtarbase.cuhk.edu.cn/
    mirtarbase_file = 'miRTarBase_MTI.xlsx'

    df = pd.read_excel(mirtarbase_file)
    validated = df[df['miRNA'] == mirna]
    return validated[['Target Gene', 'Experiments', 'Support Type']]
```


## Seed Match Analysis


**Goal:** Find miRNA seed region complementary matches within a 3' UTR sequence.

**Approach:** Extract the 7-mer seed (positions 2-8), compute its reverse complement, and scan the UTR for all occurrences.

```python
from Bio.Seq import Seq

def find_seed_matches(mirna_seq, utr_seq):
    '''Find seed matches in UTR sequence

    Seed region: positions 2-8 of miRNA (7-mer)
    '''
    mirna = Seq(mirna_seq)
    utr = Seq(utr_seq)

    # Get seed (positions 2-8, 0-indexed: 1-7)
    seed = str(mirna[1:8])
    seed_rc = str(Seq(seed).reverse_complement())

    matches = []
    start = 0
    while True:
        pos = str(utr).find(seed_rc, start)
        if pos == -1:
            break
        matches.append(pos)
        start = pos + 1

    return matches
```


## Functional Enrichment of Targets


**Goal:** Identify biological functions enriched among predicted miRNA target genes.

**Approach:** Run GO and KEGG enrichment analysis on the target gene list using Enrichr via gseapy.

```python
def enrich_target_genes(targets, background=None):
    '''Run GO enrichment on predicted target genes'''
    import gseapy as gp

    enr = gp.enrichr(
        gene_list=targets,
        gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human'],
        organism='Human'
    )
    return enr.results
```

