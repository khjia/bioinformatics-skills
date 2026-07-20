---
name: bio-epitranscriptomics-m6a
description: "m6A epitranscriptomics analysis: MeRIP-seq preprocessing, peak calling (exomePeak2, MACS2), differential m6A analysis, m6Anet ONT-based detection, and visualization (metagene plots, browser tracks)."
---


## Alignment with STAR


**Goal:** Align MeRIP-seq IP and input samples to the genome with splice-aware mapping for downstream peak calling.

**Approach:** Build a STAR genome index with gene annotations, then loop through all IP and input samples to produce coordinate-sorted BAM files.

```bash
# Build index (once)
STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf

# Align IP and input samples
for sample in IP_rep1 IP_rep2 Input_rep1 Input_rep2; do
    STAR --genomeDir star_index \
        --readFilesIn ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_
done
```


## QC Metrics


```bash
# Index BAMs
for bam in *Aligned.sortedByCoord.out.bam; do
    samtools index $bam
done

# Check IP enrichment
# Good MeRIP: IP should have peaks, input should be uniform
samtools flagstat IP_rep1_Aligned.sortedByCoord.out.bam
```


## IP/Input Correlation


```python
import deeptools.plotCorrelation as pc

# Check replicate correlation
multiBamSummary bins \
    -b IP_rep1.bam IP_rep2.bam Input_rep1.bam Input_rep2.bam \
    -o results.npz

plotCorrelation -in results.npz \
    --corMethod spearman \
    -o correlation.png
```


## Related Skills


- read-qc/quality-reports - Raw read quality assessment
- read-alignment/star-alignment - General alignment concepts
- m6a-peak-calling - Next step after preprocessing


## exomePeak2 (Recommended)


**Goal:** Identify m6A-enriched regions by comparing IP and input samples with GC-bias correction and replicate-aware statistical testing.

**Approach:** Provide IP and input BAM files along with a gene annotation to exomePeak2, which models read counts in sliding windows across the transcriptome and calls significant enrichment peaks.

```r
library(exomePeak2)

# Peak calling with biological replicates
result <- exomePeak2(
    bam_ip = c('IP_rep1.bam', 'IP_rep2.bam'),
    bam_input = c('Input_rep1.bam', 'Input_rep2.bam'),
    gff = 'genes.gtf',
    genome = 'hg38',
    paired_end = TRUE
)

# Export peaks
exportResults(result, format = 'BED')
```


## MACS3 Alternative


```bash
# Call peaks treating input as control
macs3 callpeak \
    -t IP_rep1.bam IP_rep2.bam \
    -c Input_rep1.bam Input_rep2.bam \
    -f BAMPE \
    -g hs \
    -n m6a_peaks \
    --nomodel \
    --extsize 150 \
    -q 0.05
```


## MeTPeak


```r
library(MeTPeak)

# GTF-aware peak calling
metpeak(
    IP_BAM = c('IP_rep1.bam', 'IP_rep2.bam'),
    INPUT_BAM = c('Input_rep1.bam', 'Input_rep2.bam'),
    GENE_ANNO_GTF = 'genes.gtf',
    OUTPUT_DIR = 'metpeak_output'
)
```


## Peak Filtering


```bash
# Filter by fold enrichment and q-value
# FC > 2, q < 0.05 typical thresholds
awk '$7 > 2 && $9 < 0.05' peaks.xls > filtered_peaks.bed
```


## exomePeak2 Differential Analysis


**Goal:** Identify m6A sites that differ in methylation level between experimental conditions from MeRIP-seq data.

**Approach:** Run exomePeak2 with a contrast design matrix comparing IP/input ratios across conditions, which accounts for GC bias and biological replicates.

```r
library(exomePeak2)

# Define sample design
# condition: factor for comparison
design <- data.frame(
    condition = factor(c('ctrl', 'ctrl', 'treat', 'treat'))
)

# Differential peak calling
result <- exomePeak2(
    bam_ip = c('ctrl_IP1.bam', 'ctrl_IP2.bam', 'treat_IP1.bam', 'treat_IP2.bam'),
    bam_input = c('ctrl_Input1.bam', 'ctrl_Input2.bam', 'treat_Input1.bam', 'treat_Input2.bam'),
    gff = 'genes.gtf',
    genome = 'hg38',
    experiment_design = design
)

# Get differential sites
diff_sites <- results(result, contrast = c('condition', 'treat', 'ctrl'))
```


## QNB for Differential Methylation


```r
library(QNB)

# Requires count matrices from peak regions
# IP and input counts per sample
qnb_result <- qnbtest(
    IP_count_matrix,
    Input_count_matrix,
    group = c(1, 1, 2, 2)  # 1=ctrl, 2=treat
)

# Filter significant
# padj < 0.05, |log2FC| > 1
sig <- qnb_result[qnb_result$padj < 0.05 & abs(qnb_result$log2FC) > 1, ]
```


## Visualization


```r
library(ggplot2)

# Volcano plot
ggplot(diff_sites, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed')
```


## Data Preparation


```bash
# Basecall with Guppy (requires FAST5 files)
guppy_basecaller \
    -i fast5_dir \
    -s basecalled \
    --flowcell FLO-MIN106 \
    --kit SQK-RNA002

# Align to transcriptome
minimap2 -ax map-ont -uf transcriptome.fa reads.fastq > aligned.sam
```


## Run m6Anet


```python
from m6anet.utils import preprocess
from m6anet import run_inference

# Preprocess: extract features from FAST5
preprocess.run(
    fast5_dir='fast5_pass',
    out_dir='m6anet_data',
    reference='transcriptome.fa',
    n_processes=8
)

# Run m6A inference
run_inference.run(
    input_dir='m6anet_data',
    out_dir='m6anet_results',
    n_processes=4
)
```


## CLI Workflow


**Goal:** Run the complete m6Anet pipeline from FAST5 signal data to per-site m6A modification probabilities.

**Approach:** First extract features from FAST5 files with dataprep (signal-to-feature extraction), then run neural network inference to classify each DRACH motif site as modified or unmodified.

```bash
# Preprocess
m6anet dataprep \
    --input_dir fast5_pass \
    --output_dir m6anet_data \
    --reference transcriptome.fa \
    --n_processes 8

# Inference
m6anet inference \
    --input_dir m6anet_data \
    --output_dir m6anet_results \
    --n_processes 4
```


## Interpret Results


```python
import pandas as pd

results = pd.read_csv('m6anet_results/data.site_proba.csv')

# Filter high-confidence m6A sites
# probability > 0.9: High confidence threshold
m6a_sites = results[results['probability_modified'] > 0.9]
```


## Metagene Plots with Guitar


```r
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load m6A peaks
peaks <- import('m6a_peaks.bed')

# Create metagene plot
# Shows distribution relative to transcript features
GuitarPlot(
    peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    saveToPDFprefix = 'm6a_metagene'
)
```


## Custom Metagene with deepTools


**Goal:** Create a metagene profile showing m6A enrichment distribution relative to gene body landmarks (TSS, TES).

**Approach:** Compute the log2 IP/input ratio as a bigWig track with bamCompare, then build a signal matrix over scaled gene regions with computeMatrix and render as a profile plot.

```bash
# Create bigWig from IP/Input ratio
bamCompare -b1 IP.bam -b2 Input.bam \
    --scaleFactors 1:1 \
    --ratio log2 \
    -o IP_over_Input.bw

# Metagene around stop codons
computeMatrix scale-regions \
    -S IP_over_Input.bw \
    -R genes.bed \
    --regionBodyLength 2000 \
    -a 500 -b 500 \
    -o matrix.gz

plotProfile -m matrix.gz -o metagene.pdf
```


## Browser Tracks


```bash
# Create normalized bigWig for genome browser
bamCoverage -b IP.bam \
    --normalizeUsing CPM \
    -o IP_normalized.bw

# Peak BED to bigBed
bedToBigBed m6a_peaks.bed chrom.sizes m6a_peaks.bb
```


## Heatmaps


```r
library(ComplexHeatmap)

# m6A signal around peaks
Heatmap(
    signal_matrix,
    name = 'm6A signal',
    cluster_rows = TRUE,
    show_row_names = FALSE
)
```

