---
name: rna-seq-pipeline
description: End-to-end bulk RNA-seq pipeline from ENA download to DEG interaction analysis. Covers download (parallel curl), audit (metadata check), fastp QC summary, PCA-based sample selection, and DEG interaction pattern classification.
tags: [rna-seq, pipeline, download, audit, fastp, pca, deg, interaction]
---

# Bulk RNA-seq Pipeline

Covers the complete RNA-seq analysis flow from **public data download** to **differential expression interpretation**, organized as 5 independent stages you can run sequentially or standalone.

## Pipeline Overview

```
Stage 0: Download     → ENA parallel curl + MD5 verify
Stage 1: Audit        → library_strategy check, pairing, completeness  
Stage 2: QC Summary   → fastp JSON + HISAT2 log → unified TSV
Stage 3: PCA Selection→ representative sample picking for annotation
Stage 4: DEG Patterns → interaction-term classification (A/B/C/D)
```

Each stage writes to a numbered project directory (e.g. `00.rawdata`, `02.cleandata`, `03.Mapping`, `06.DEG`).

---

## Stage 0 — ENA Parallel Download

Download FASTQ from ENA/EBI over HTTPS, parallel with MD5 verification.

### Why ENA over NCBI SRA

ENA mirrors the same runs as `.fastq.gz` over HTTPS with MD5 checksums. Skip the `.sra → fasterq-dump` step entirely.

**Bandwidth (CNS cluster, sing-box VLESS proxy, May 2026):**

| Mode | MB/s |
|------|------|
| 1× curl direct (no proxy) | 1.09 |
| 1× FTP protocol | 0.13 |
| 1× HTTPS via proxy | 1.0–5.2 |
| **4× parallel HTTPS via proxy** | **4.7–5.9 ✅** |
| 8× parallel HTTPS via proxy | 5.4 (no benefit) |

### Step 0.1 — Get URLs + MD5 from ENA filereport

```bash
# Single accession
curl -sS -x http://127.0.0.1:10809 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR21846941&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,instrument_platform&format=tsv"

# Batch via build_worklist.py
python3 scripts/build_worklist.py \
  --projects PRJNA267198 PRJNA887921 \
  --dst /path/to/01.rawdata \
  --out metadata/worklist.txt \
  --layout PAIRED
```

### Step 0.2 — Parallel download

```bash
META=/path/to/metadata \
WORKLIST=/path/to/worklist.txt \
LOG=/path/to/download.log \
nohup bash scripts/download_parallel.sh > /dev/null 2>&1 &
```

### Step 0.3 — Monitor

```bash
# Concurrency check
pgrep -c "^curl$"

# Speed sampling (60s)
bash scripts/speed_check.sh /path/to/01.rawdata 60

# Progress
tail -5 /path/to/download.log
```

### Step 0.4 — Verify

```bash
expected=$(wc -l < metadata/worklist.txt)
got=$(ls 01.rawdata/*.fastq.gz | wc -l)
echo "expected=$expected got=$got"

# MD5 sweep
awk -F'|' '{print $3"  "$2"/"gensub(/.*\//,"","g",$1)}' metadata/worklist.txt \
  | md5sum -c - 2>&1 | grep -v ': OK$'
```

### Pitfalls

1. **Orphan xargs**: `pkill -f download_parallel.sh` does NOT kill xargs children. Use `kill -9 <xargs_pid>`.
2. **xargs field splitting**: Use `|` delimiter, NOT TAB. Pass `{}` as `$1` via `_ {}`.
3. **--max-time**: Set to 7200s for large files. Each retry restarts from scratch.
4. **HTTPS not FTP**: Even with `ftp.sra.ebi.ac.uk`, use `https://` prefix (40× faster).
5. **FAKE PAIRED in ENA**: Some runs tagged PAIRED but have only one `.fastq.gz` — check `fastq_ftp` for `;` separator or NCBI Trace `nreads`.
6. **New Bioproject (<24h)**: ENA may be blank. Use NCBI Trace per-run API: `trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=SRR...`.
7. **Path buckets**: Never construct ENA URLs by hand. Always read `fastq_ftp` from filereport.

---

## Stage 1 — Download Audit

Verify downloaded data is genuine mRNA-seq, properly paired, and complete.

### 1.1 File completeness

```bash
RAWDIR=/path/to/01.rawdata

# R1 vs R2 parity
ls $RAWDIR/*_1.fastq.gz 2>/dev/null | sed 's/_1.fastq.gz//' | sort > /tmp/r1.txt
ls $RAWDIR/*_2.fastq.gz 2>/dev/null | sed 's/_2.fastq.gz//' | sort > /tmp/r2.txt
echo "R1 without R2:"; comm -23 /tmp/r1.txt /tmp/r2.txt
echo "R2 without R1:"; comm -13 /tmp/r1.txt /tmp/r2.txt

# Empty files
find $RAWDIR -name "*.fastq.gz" -size 0

# Truncated files (<10MB for PE RNA-seq)
find $RAWDIR -name "*.fastq.gz" -size -10M
```

### 1.2 Metadata verification (public data only)

Batch-query ENA for library_strategy, library_source, library_layout:

```bash
cut -f2 samples.tsv | grep -P '^(SRR|ERR|DRR)' | sort > /tmp/all_runs.txt
bash scripts/check_ena_strategy.sh /tmp/all_runs.txt /tmp/ena_check.tsv
```

Analyze:

```bash
# Strategy distribution
tail -n +2 /tmp/ena_check.tsv | cut -f2 | sort | uniq -c | sort -rn

# Non-RNA-Seq runs (exclude candidates)
tail -n +2 /tmp/ena_check.tsv | awk -F'\t' '$2!="RNA-Seq" && $2!="ssRNA-seq"'
```

### 1.3 Metadata interpretation

| library_strategy | Verdict |
|------------------|---------|
| RNA-Seq | ✅ mRNA |
| ssRNA-seq | ✅ mRNA (strand-specific) |
| WGS / WES / OTHER | ❌ Exclude |

**Key rule**: `library_strategy` is ground truth, NOT `library_source`. Some submitters tag RNA-Seq as GENOMIC by mistake.

### 1.4 Cleanup orphan files

```bash
comm -13 /tmp/target_runs.txt /tmp/actual_runs.txt | while read run; do
  echo "Removing orphan: $run"
  rm -f rawdata/${run}_*.fastq.gz
done
```

---

## Stage 2 — QC Summary Table

After running fastp + HISAT2, generate a unified QC summary table.

### Prerequisites

- fastp JSON reports in `02.cleandata/reports/{sample}.fastp.json`
- HISAT2 logs in `03.Mapping/stats/{sample}.hisat2.log`
- Sample annotation TSV with tissue/treatment columns

### Workflow

```python
import json, os, pandas as pd
from glob import glob

reports_dir = "02.cleandata/reports"
stats_dir = "03.Mapping/stats"
ann_file = "download/metadata/sample_annotation.tsv"

ann = pd.read_csv(ann_file, sep='\t', index_col=0) if os.path.exists(ann_file) else {}

rows = []
for f in sorted(glob(f"{reports_dir}/*.fastp.json")):
    sample = os.path.basename(f).replace(".fastp.json","")
    with open(f) as fh:
        j = json.load(fh)
    raw = j['summary']['before_filtering']['total_reads']
    clean = j['summary']['after_filtering']['total_reads']
    retention = clean / raw * 100 if raw > 0 else 0
    
    mapping = None
    logfile = f"{stats_dir}/{sample}.hisat2.log"
    if os.path.exists(logfile):
        with open(logfile) as fh:
            for line in fh:
                if 'overall alignment' in line:
                    mapping = float(line.split()[0].replace('%',''))
                    break
    
    meta = ann.loc[sample] if sample in ann.index else {}
    rows.append({
        'Sample': sample, 'Raw_reads': raw, 'Clean_reads': clean,
        'Retention%': round(retention,1), 'Mapping%': mapping or '',
        'Tissue': meta.get('tissue',''), 'Treatment': meta.get('treatment',''),
    })

pd.DataFrame(rows).to_csv('qc_summary.tsv', sep='\t', index=False)
```

### Interpretation

| Mapping% | Assessment |
|----------|------------|
| >95% | Excellent |
| 90-95% | Good |
| 80-90% | Acceptable (check rRNA/viral) |
| <80% | Investigate |

---

## Stage 3 — PCA Sample Selection

Pick representative samples from a large count matrix for functional annotation (eggNOG, GO, KEGG).

### Criteria

- Clean data 4-10 GB
- Mapping rate ≥90%
- Diverse in PCA space (cover expression landscape)
- Cover multiple tissues and treatments
- Local samples prioritized

### 3.1 Run PCA on count matrix

```python
import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

df = pd.read_csv('raw_counts_matrix.tsv', sep='\t', index_col=0)
df = df[df.sum(axis=1) > 0]

lib = df.sum(axis=0)
cpm = df.div(lib, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

pca = PCA(n_components=5)
X = StandardScaler().fit_transform(log_cpm.T)
pcs = pca.fit_transform(X)

pca_df = pd.DataFrame(pcs[:,:2], index=log_cpm.columns, columns=['PC1','PC2'])
pca_df = pca_df.join(ann[['tissue','treatment_category']], how='left')
pca_df.to_csv('pca_with_metadata.tsv', sep='\t')
```

### 3.2 Filter and select

```python
# Load QC data
sizes = {}  # Sample → clean GB
mrates = {}  # Sample → mapping%

candidates = {s for s in pca_df.index 
              if 4 <= sizes.get(s, 0)/1e9 <= 10 
              and mrates.get(s, 0) >= 90}

# Prioritize local samples, then pick diverse public samples
selected = {s for s in candidates if not s.startswith(('SRR','ERR','DRR'))}

# Pick 1-2 per tissue×treatment, maximizing PCA distance
grouped = pca_df.loc[list(candidates - selected)].groupby(['tissue','treatment_category'])
for (tissue, trt), grp in grouped:
    grp = grp.copy()
    grp['dist'] = np.sqrt(grp['PC1']**2 + grp['PC2']**2)
    for s in grp.nlargest(2, 'dist').index:
        if len(selected) < 20:
            selected.add(s)

# Output
out = '\n'.join([f'{s}\tlocal' if not s.startswith(('SRR','ERR','DRR')) 
                 else f'{s}\tpublic' for s in sorted(selected)])
with open('selected_samples.tsv', 'w') as f:
    f.write('Sample\tType\n' + out)
```

### Pitfalls

1. **Merging tissue names**: Normalize variant tissue names (leaf/leaves/叶片) with a mapping dict.
2. **Low PC variance**: For diverse datasets PC1+PC2 may be 15-20%. Normal — within-group variance captured by higher PCs.
3. **Manual curation**: Auto-selection is a starting point. Review and adjust manually.

---

## Stage 4 — DEG Interaction Pattern Classification

After running DESeq2 with a `genotype:treatment` interaction term, classify significant interaction genes by biological response pattern.

### The Problem

DESeq2 interaction log2FoldChange sign alone only says "up" or "down" — insufficient to interpret how genotypes differ in their treatment response.

### 4 Patterns

| Pattern | Name | Criteria |
|---------|------|----------|
| **A** | Genotype1-specific | Genotype1 `padj<0.05 & |LFC|≥1`, Genotype2 not sig |
| **B** | Genotype2-specific | Genotype2 `padj<0.05 & |LFC|≥1`, Genotype1 not sig |
| **C** | Opposite direction | Both sig, opposite sign |
| **D** | Magnitude diff | Both sig, same sign, `|LFC_diff| ≥ 1` |
| Other | Weak | Interaction sig but neither genotype individually sig |

### Input files

```
DESeq2 results:
  interaction_DEGs.tsv           (interaction-term sig genes)
  genotype1_salt_vs_ctrl_all.tsv (per-genotype LFC)
  genotype2_salt_vs_ctrl_all.tsv

→ LEFT JOIN on gene ID → classify by pattern → output
```

### Parameters

```r
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
LFC_DIFF_CUTOFF <- 1
```

### Output

| File | Content |
|------|---------|
| `interaction_patterns.tsv` | Full table with per-genotype LFC, sig, direction, pattern |
| `interaction_pattern_summary.tsv` | Count per pattern |
| `interaction_A_*.txt` | Gene list per pattern |

### Verified example (Setaria italica, guzi project)

2,236 interaction-significant genes:
- A (genotype1-specific): 555 (24.8%)
- B (genotype2-specific): 653 (29.2%)
- C (opposite direction): 32 (1.4%)
- D (magnitude diff): 725 (32.4%)
- Other (weak): 271 (12.1%)

### Pitfalls

1. Gene ID format must match between interaction and per-genotype files
2. C pattern naturally rare (opposite response to same stress is unusual)
3. Sankey plots need ≥2 ontologies to work

---

## Scripts

The following scripts are bundled with this skill:

| Script | Stage | Purpose |
|--------|-------|---------|
| `scripts/download_parallel.sh` | 0 | Parallel ENA FASTQ download with MD5 verify |
| `scripts/build_worklist.py` | 0 | Build worklist from BioProject accessions |
| `scripts/speed_check.sh` | 0 | Live download speed monitoring |
| `scripts/check_ena_strategy.sh` | 1 | Batch ENA metadata query |

## Citations

### Tools

- fastp: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 2018, 34(17):i884–i890. doi:10.1093/bioinformatics/bty560
- HISAT2: Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 2019, 37:907–915. doi:10.1038/s41587-019-0201-4
- DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 2014, 15:550. doi:10.1186/s13059-014-0550-8

### Data sources

- ENA: Burgin J, Ahamed A, Cummins C, et al. The European Nucleotide Archive in 2022. Nucleic Acids Research, 2023, 51(D1):D121–D125. doi:10.1093/nar/gkac1051

### Python libraries

- scikit-learn: Pedregosa F, et al. Scikit-learn: Machine Learning in Python. JMLR, 2011, 12:2825–2830.
- pandas: McKinney W. Data Structures for Statistical Computing in Python. SciPy, 2010.

### Methods template

> RNA-seq raw FASTQ files were downloaded from the European Nucleotide Archive (Burgin et al., 2023) using parallel curl with MD5 verification. Adapter trimming and quality filtering were performed with fastp (Chen et al., 2018) using default parameters. Reads were aligned to the reference genome with HISAT2 (Kim et al., 2019). Gene-level counts were quantified with featureCounts and differential expression analysis was performed using DESeq2 (Love et al., 2014) with the design ~ condition. PCA-based sample selection was used to select representative samples for functional annotation.

## 中文使用说明

本 skill 将 RNA-seq 分析流程分为 5 个独立阶段，每个阶段可以单独跑也可以串联跑：

| 阶段 | 功能 | 输入 | 输出 |
|------|------|------|------|
| Stage 0 | 从 ENA 并行下载 fastq.gz | BioProject/SRR 号列表 | rawdata/ + 校验 |
| Stage 1 | 审计下载数据 | rawdata/ 目录 | 元数据核查表 |
| Stage 2 | QC 汇总表 | fastp JSON + HISAT2 日志 | qc_summary.tsv |
| Stage 3 | PCA 选代表性样本 | count 矩阵 + QC 指标 | selected_samples.tsv |
| Stage 4 | DEG 交互效应分类 | DESeq2 结果 | interaction_patterns.tsv |

**设计原则**：按实际工作流顺序编号，但每个阶段独立可复用。Stage 0 和 1 用于公共数据下载场景，Stage 2-4 适用于任何 RNA-seq 项目（公共或自测）。
