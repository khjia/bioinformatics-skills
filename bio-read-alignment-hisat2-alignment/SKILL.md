---
name: bio-read-alignment-hisat2-alignment
description: Align RNA-seq reads with HISAT2, a memory-efficient splice-aware aligner. Use when STAR's memory requirements are too high or for general RNA-seq alignment.
tool_type: cli
primary_tool: HISAT2
---

## Version Compatibility

Reference examples tested with: samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# HISAT2 RNA-seq Alignment

**"Align RNA-seq reads with HISAT2"** → Map RNA-seq reads to a reference genome with splice-aware alignment. Suitable for gene expression quantification workflows.
- CLI: `hisat2 -x index -1 R1.fq -2 R2.fq | samtools sort -o aligned.bam`

## Build Index

### GFF3 vs GTF (critical pitfall)

**`hisat2_extract_splice_sites.py` and `hisat2_extract_exons.py` require GTF format with `gene_id` and `transcript_id` attributes.** They do NOT understand GFF3 format (which uses `ID` and `Parent` instead). If you feed them a GFF3 file, they silently produce 0-byte output with no error — and the index builds without splice site info.

**Always convert GFF3 to GTF first:**

```bash
# Install gffread (if not already in env)
mamba install -n hisat2 -c bioconda gffread -y

# Convert
gffread genome.gff3 -T -o genome.gtf

# Verify conversion has gene_id/transcript_id
head -3 genome.gtf
# Should show: ... gene_id "XXX"; transcript_id "XXX";
```

After conversion, the extraction scripts will work correctly and produce non-empty splice_sites.txt and exons.txt.

### With annotation (recommended)

```bash
hisat2_extract_splice_sites.py annotation.gtf > splice_sites.txt
hisat2_extract_exons.py annotation.gtf > exons.txt

# Verify output is non-empty BEFORE building index
wc -l splice_sites.txt exons.txt
# Both should have thousands of lines

hisat2-build -p 16 \
    --ss splice_sites.txt \
    --exon exons.txt \
    reference.fa hisat2_index
```

### Without annotation

```bash
hisat2-build -p 16 reference.fa hisat2_index
```

## Basic Alignment

```bash
# Paired-end reads
hisat2 -p 8 -x hisat2_index \
    -1 reads_1.fq.gz -2 reads_2.fq.gz \
    -S aligned.sam

# Single-end reads
hisat2 -p 8 -x hisat2_index \
    -U reads.fq.gz \
    -S aligned.sam
```

## Direct to Sorted BAM

```bash
# Pipe to samtools
hisat2 -p 8 -x hisat2_index \
    -1 r1.fq.gz -2 r2.fq.gz | \
    samtools sort -@ 4 -o aligned.sorted.bam -

samtools index aligned.sorted.bam
```

## Stranded Libraries

```bash
# Forward stranded (e.g., Ligation)
hisat2 -p 8 -x hisat2_index \
    --rna-strandness FR \
    -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam

# Reverse stranded (e.g., dUTP, TruSeq - most common)
hisat2 -p 8 -x hisat2_index \
    --rna-strandness RF \
    -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam

# Single-end stranded
hisat2 -p 8 -x hisat2_index \
    --rna-strandness F \    # or R for reverse
    -U reads.fq.gz -S aligned.sam
```

## Novel Splice Junction Discovery

```bash
# Output novel splice junctions
hisat2 -p 8 -x hisat2_index \
    --novel-splicesite-outfile novel_splices.txt \
    -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam

# Use known + novel junctions for subsequent alignments
hisat2 -p 8 -x hisat2_index \
    --novel-splicesite-infile novel_splices.txt \
    -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam
```

## Two-Pass Alignment (Manual)

**Goal:** Improve splice junction sensitivity by discovering novel junctions across all samples in a first pass, then realigning with the combined junction set.

**Approach:** Run HISAT2 on each sample to extract novel splice sites, merge and deduplicate junctions across samples, then realign all samples using the combined junction catalog.

```bash
# Pass 1: Discover junctions from all samples
for r1 in *_R1.fq.gz; do
    r2=${r1/_R1/_R2}
    base=$(basename $r1 _R1.fq.gz)
    hisat2 -p 8 -x hisat2_index \
        --novel-splicesite-outfile ${base}_splices.txt \
        -1 $r1 -2 $r2 -S /dev/null
done

# Combine and filter junctions
cat *_splices.txt | sort -u > combined_splices.txt

# Pass 2: Realign with all junctions
for r1 in *_R1.fq.gz; do
    r2=${r1/_R1/_R2}
    base=$(basename $r1 _R1.fq.gz)
    hisat2 -p 8 -x hisat2_index \
        --novel-splicesite-infile combined_splices.txt \
        -1 $r1 -2 $r2 | \
        samtools sort -@ 4 -o ${base}.sorted.bam -
done
```

## Read Group Information

```bash
hisat2 -p 8 -x hisat2_index \
    --rg-id sample1 \
    --rg SM:sample1 \
    --rg PL:ILLUMINA \
    --rg LB:lib1 \
    -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam
```

## Downstream Quantification

```bash
# Output name-sorted BAM for htseq-count
hisat2 -p 8 -x hisat2_index -1 r1.fq.gz -2 r2.fq.gz | \
    samtools sort -n -@ 4 -o aligned.namesorted.bam -

# Or coordinate-sorted for featureCounts
hisat2 -p 8 -x hisat2_index -1 r1.fq.gz -2 r2.fq.gz | \
    samtools sort -@ 4 -o aligned.sorted.bam -
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| -p | 1 | Number of threads |
| -x | - | Index basename |
| --rna-strandness | unstranded | FR/RF/F/R |
| --dta | off | Downstream transcriptome assembly |
| --dta-cufflinks | off | For Cufflinks |
| --min-intronlen | 20 | Minimum intron length |
| --max-intronlen | 500000 | Maximum intron length |
| -k | 5 | Max alignments to report |

## For StringTie/Cufflinks

```bash
# Use --dta for StringTie
hisat2 -p 8 -x hisat2_index \
    --dta \
    -1 r1.fq.gz -2 r2.fq.gz | \
    samtools sort -@ 4 -o aligned.sorted.bam -
```

## Alignment Summary

```bash
# HISAT2 prints summary to stderr
hisat2 -p 8 -x hisat2_index -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam 2> summary.txt
```

Example:
```
50000000 reads; of these:
  50000000 (100.00%) were paired; of these:
    2500000 (5.00%) aligned concordantly 0 times
    45000000 (90.00%) aligned concordantly exactly 1 time
    2500000 (5.00%) aligned concordantly >1 times
95.00% overall alignment rate
```

## Memory Comparison

| Aligner | Human Genome Memory |
|---------|-------------------|
| STAR | ~30GB |
| HISAT2 | ~8GB |

## Related Skills

- read-alignment/star-alignment - Alternative with more features
- rna-quantification/featurecounts-counting - Count aligned reads
- rna-quantification/alignment-free-quant - Skip alignment entirely
- differential-expression/deseq2-basics - Downstream DE analysis

## Pitfalls (proven on CNS cluster, 2026)

### 1. NFS corrupts large HISAT2 index files

Building the index directly on NFS (e.g., /media/nfs1 or /media/nfs2) can produce corrupted .ht2 files — the final .2.ht2 or .5/.6 files end up as 0 bytes or wrong size. Error: `Index is corrupt: File size for ... should have been X but is actually 0`.

**Fix:** Build on local `/tmp` on CNS2, then copy to NFS:

```bash
TMPDIR=/tmp/hisat2_build_$$ && mkdir -p $TMPDIR
cp splice_sites.txt exons.txt $TMPDIR/
hisat2-build -p 12 --ss $TMPDIR/splice_sites.txt --exon $TMPDIR/exons.txt genome.fa $TMPDIR/index_name
cp $TMPDIR/index_name.*.ht2 /path/to/nfs/index_dir/
rm -rf $TMPDIR
```

NFS1 at >95% capacity is especially prone to this corruption. Even at lower utilization, parallel writes from multiple processes exacerbate it.

### 2. BAI index fails for chromosomes > 512 Mb

BAI format has a hard limit of 2^29 = 536,870,912 bases. Chromosomes exceeding this (e.g., pea chr03 = 553Mb) cause `samtools index` to fail with `Numerical result out of range`.

**Fix:** Always use CSI index (`samtools index -c`) for large plant genomes:

```bash
samtools sort -@ 4 -o aligned.bam -
samtools index -c aligned.bam    # produces .csi, not .bai
```

For HTCondor wrappers, check for `.csi` not `.bai` in resume logic:

```bash
if [ -s "$bam" ] && [ -s "${bam}.csi" ]; then
  echo "[skip] $sample already done"; exit 0
fi
```

### 3. `--dta` memory impact on large genomes

For genomes > 3Gb, the `--dta` flag adds ~10+ GB of memory usage. On the 4.5Gb pea genome with 4 threads, HISAT2 needed >24GB with `--dta`. Without `--dta`, 24GB is sufficient.

**Only use `--dta` if downstream needs StringTie/Cufflinks transcript assembly. For featureCounts-based quantification, omit it.**

### 4. Condor defaults memory to 128 MB

If `request_memory` is not specified in a `.condor` file, HTCondor on this cluster defaults to 128 MB per job — instantly holding any bioinformatics tool.

**Always specify `request_memory`.** Start with 24 GB for HISAT2 on large plant genomes without `--dta`. With `--dta`, start at 32-40 GB.

### 5. Clean up failed condor logs before resubmitting

When a cluster fails and is removed, old log files in `logs/` persist. The next submission appends to them if filenames collide (e.g., `hisat2.$(idx).out`). Always `rm -f logs/*` before re-submitting the same cluster with fixes.
