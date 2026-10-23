---
name: rnaseq-hisat2-featurecounts-deseq2
description: "End-to-end bulk RNA-seq analysis pipeline: HISAT2 alignment → featureCounts quantification + TPM + tissue-level low-expression filtering → DESeq2 differential expression → clusterProfiler enrichment (GO/KEGG). Covers GFF3→GTF conversion, NFS pipe pitfalls, CSI indexing for large genomes, and tissue-aware TPM filtering."
tags: [rna-seq, hisat2, featureCounts, TPM, DESeq2, enrichment, plant, pipeline]
---

# RNA-seq HISAT2 + featureCounts + DESeq2 Pipeline

End-to-end bulk RNA-seq analysis pipeline for plant genomes using HISAT2 alignment, featureCounts quantification with tissue-level TPM filtering, DESeq2 differential expression, and clusterProfiler enrichment.

## Pipeline Overview

```
1.rawdata/    → fastp QC
2.clean/      → HISAT2 alignment
3.Mapping/    → featureCounts quantification
4.Quantification/  → TPM + tissue-level filtering
5.DEG/        → DESeq2 differential expression
6.Enrichment/ → clusterProfiler GO/KEGG enrichment
```

## Stage 1: Read QC (fastp)

### Pitfall: Report sequencing volume in BASES (Gb/Tb), NOT compressed-file bytes (GB)

When reporting "how much data" (raw or clean), do NOT use the `.fastq.gz` file size on disk — gz compression shrinks it to ~1/3 of the true base count. Users expect sequencing volume in **bases** (e.g. "~6 Gb/sample" typical for RNA-seq). Reporting compressed bytes (e.g. "0.59 TB") when the real base count is ~1.02 Tb looks like the data is far smaller than it should be, and the user will catch it.

**Correct source of base counts — read from the fastp JSON, not by re-decompressing:**

```python
import json, os
d = json.load(open(f"{sid}_fastp.json"))
after = d["summary"]["after_filtering"]
clean_bases = after["total_bases"]      # TRUE clean base count (bp)
clean_reads = after["total_reads"]      # single-end reads (R1+R2 both counted)
q30 = after["q30_rate"]
# Raw (before) is in d["summary"]["before_filtering"]["total_bases"]
```

- `total_bases` = true base count (bp). Sum across samples → total clean volume in Gb/Tb.
- The fastp JSON has both `before_filtering` and `after_filtering` blocks → get raw AND clean bases without re-reading files.
- **Do NOT re-run `zcat | awk 'NR%4==2{...}'`** to count bases for many samples — it re-decompresses everything (slow). fastp already recorded it.
- A sample with low clean bases is usually because **it was shallow-sequenced or short-read** (e.g. 75bp reads, 4M pairs → 0.3 Gb), NOT because fastp filtered heavily. Verify via `after.total_reads/before.total_reads` — if pass rate is ~98%, no data was lost; the sample was just small to begin with. Also note compressed raw file "3.4G" ≠ 0.3 Gb bases (compression again).

For the QC summary table, report **clean bases (Gb)** per sample (from `total_bases`), plus read count, Q30, and mapping rate. See the fig project's `3.Mapping/summary/clean_and_mapping_summary.tsv` as a template.

### Pitfall: fastp JSON `total_reads` counts single-end reads, NOT PE pairs — verify data volume via HISAT2

When reporting raw/filtered read counts for **paired-end** data, do NOT trust `fastp` JSON `summary.before_filtering.total_reads` directly as the PE-pair count. In fastp, for PE input `total_reads` = **single-end read count (R1+R2 both counted)** — i.e. it is ~2× the number of PE pairs. Confusion here produces wrong data-volume numbers.

**Correct way to compute PE pairs & data volume:**
- HISAT2 log `total_reads` (e.g. "18236078 reads") = **PE pairs** (each pair counted once). Sum across samples = total PE pairs.
- fastp `total_reads` ≈ HISAT2 total_reads × 2 (single-end reads).
- Data volume = PE pairs × 2 × read_length. Verify read length by peeking at raw FASTQ (`gzip.open(...)`; don't assume 150bp — check `Counter` of read lengths), not by `total_bases/total_reads` (fastp truncates low-quality tails so that division gives a misleading shorter length).
- 珊瑚菜 example: fastp showed total_reads ≈ 556M (single-end); actual PE pairs = 278M; PE150 → 82.6 Gb (NOT 166.9 Gb which is the wrong single-end×2 calculation).

In the QC summary table report the **PE pairs** column (from HISAT2) and label it clearly, so the read-count column is unambiguous.

### Pitfall: Detect "fake-PE" samples before mapping (public data QC)

Public RNA-seq runs are sometimes **labeled PAIRED in ENA/NCBI but actually single-end in SRA** (submitter mislabeled). This is common and must be caught before deciding which samples to keep for paired-end analysis. The tell-tale is the FASTQ read header read-index field:

```bash
# True interleaved-PE (one spot, two reads): headers alternate read-index 1,2,1,2...
# @SRR.1 1 length=150
# @SRR.1 2 length=150
# @SRR.2 1 length=150
# @SRR.2 2 length=150

# Fake-PE / actually single-end: each spot has ONE read, read-index increments 1,2,3,4...
# @SRR.1 1 length=150
# @SRR.2 2 length=150
# @SRR.3 3 length=150
# @SRR.4 4 length=150
```

Detection script — check the read-index field pattern on a sample:

```python
import gzip, re
def read_index_sequence(fname, n=12):
    idxs = []
    with gzip.open(fname, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                m = re.match(r'^\S+\.\d+ (\d+) length', line)
                if m:
                    idxs.append(m.group(1))
                if len(idxs) >= n:
                    break
    return idxs
# idxs == ['1','2','3',...] incrementing (not 1,2,1,2) => single-end despite PAIRED label
```

Consequences & handling:
- `fasterq-dump --split-files` on these produces only ONE file (no `_1`/`_2` pair) — you cannot manufacture a mate pair that doesn't exist.
- If the project requires **strictly paired-end**, these samples must be **excluded** (delete their files + `.done` marker), not kept as broken pairs.
- ENA `fastq_ftp` for these shows a single `.fastq.gz` (no `_1`/`_2`), confirming the mislabel even before download.
- When the user says "keep only paired-end," verify **all** samples this way; a few samples may be genuinely single-end despite being in the PE project list.

### Getting platform info when ENA/NCBI web APIs are blocked

On the CNS cluster (outbound HTTPS to ENA/NCBI often blocked), extract platform from downloaded `.sra` files directly:
```bash
vdb-dump /path/to/SRR123.sra --info   # prints "platf : SRA_PLATFORM_ILLUMINA" / _PACBIO_SMRT etc.
```
NCBI `esummary` needs the **numeric UID** (get via `esearch` first), not the SRR accession directly (`Invalid uid` error otherwise). ENA `filereport` returns empty when the network is blocked — don't rely on it.

```bash
# Standard paired-end RNA-seq filtering — 与重测序流程统一采用 duzhong 标准配置
fastp \
    -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o sample_clean_R1.fastq.gz -O sample_clean_R2.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --n_base_limit 5 \
    --average_qual 20 \
    --length_required 60 \
    --compression 6 \
    -w 16 \
    -h sample_fastp.html -j sample_fastp.json
```

DNBSEQ/MGI 数据接头不是 Illumina（reads header `@ML…`），要显式指定时加：

```bash
    --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
```

### Pitfall: fastp 默认参数几乎不过滤（2026-10 实测，duzhong 面板 SRR37245294 前 40 万 read 对）

**`-q`/`-u` 只做"整条 read 判定"，不裁单个碱基。** 裁碱基必须显式给 `--cut_front`/`--cut_tail`/
`--cut_right`。fastp 的默认值 `-q 15 -u 40 -l 15` 极其宽松——`-q 15` 意为"只有不到 40% 的碱基
低于 Q15 才丢整条 read"，实测下来**整批一条不丢**（lowq=0），Q30 一点不提升。

同一份真实数据（Q30 94.62%）跑 4 套参数：

| 配置 | bases 保留 | 平均读长 | Q30 | 被质量剔除 | 3'端 Q<20 碱基数/20万read |
|------|-----------|---------|-----|-----------|----------------------|
| **fastp 默认**（无参数） | 99.62% | 149.0 | 94.62 | 0 | 2803 |
| `--detect_adapter_for_pe` only | 99.61% | 149.0 | 94.62 | 0 | 2803 |
| **旧转录组参数**（`-q 15`，无 cut） | 99.61% | 149.0 | 94.62 | **0** | 2800 |
| **duzhong 标准**（本流程现用） | 99.54% | 149.0 | 94.66 | 0 | 2436 |
| duzhong + `-u 20` | 95.74% | 149.0 | 95.99 | 30,674 | — |

读法与坑：

- **fastp 默认就已经去接头了**（`--disable_adapter_trimming` 才是关闭）。所以"显式写接头"
  相对默认的增量很小：`--detect_adapter_for_pe` 只多探到 1,234 个 read 的接头、多裁 14,063 个碱基。
  显式写接头的价值主要是**跨平台可移植**（换数据源不必重测），不是清洁度。
- **`--cut_front --cut_tail` 是温和裁剪**：两端各自裁到合格为止，不会砍掉大半条 read。
  相比 `--cut_right`（遇一个低质量窗口就砍掉它**及其右侧全部**）损失小得多。
- **`-u 20` 才是让"低质量 read 被判死"的那一项**。duzhong 原配置没写 `-u`（默认 u40），
  在这批高质量数据上表现几乎等同默认（只丢 526 条 <36bp 短 read，质量剪切实际裁掉 0 个碱基）。
  **数据质量差的项目（Q30 < 90%）应显式加 `--unqualified_percent_limit 20`**，
  但要接受误杀（见下）。
- **`--length_required 60`（2026-10-23 收紧，原 36）= PE150 的 40%。** 实测（无花果
  SRR10992859，149 bp RNA-seq，20 万 read 对，用已清洗数据隔离长度门槛的净效应）：
  36 → 60 只多丢 902 条短 read（0.23%），bases 保留 99.99% → 99.91%，Q30 不变。
  挡掉的是接头二聚体与剪切后的残端，不是有用数据。**若文库插入片段普遍 < 120 bp
  （降解样本、3' 端严重偏倚），降到 50；不要低于 36。**
- `--compression 6`（duzhong 值，默认 4）：clean FASTQ 体积更小，CPU 略增，不影响结果。
- 接头序列若非 DNBSEQ/MGI（reads header 不是 `@ML…`），不要沿用本文档的显式接头串，
  让 `--detect_adapter_for_pe` 自动探测即可。

### 低质量数据上 `-u` 的真实效果与误杀代价（2026-10-23 实测，无花果 SRR5933780）

该样本原始 **Q30 仅 83.47%**（面板里最低的一档），用它实测 `-u` 的取舍：

| 配置 | bases 保留 | Q30 | 低质 read 被剔 | 短 read 被剔 |
|------|-----------|-----|--------------|------------|
| fastp 默认 | 100.00% | 87.27→87.27 | 0 | 0 |
| 本流程现用 (`-u 40`) | 99.58% | 87.27→87.34 | 2 | 5,160 |
| `-u 20` | 84.88% | 87.27→**90.91** | **130,478** | 2,340 |
| `-u 10` | 63.94% | 87.27→94.25 | 301,168 | 1,182 |

**但 `-u 20` 有明确的误杀代价**（这是选 `-u` 时必须知道的事）：

- `-u 20` 剔掉的 63,828 条 read 里，有 **15,970 条（25%）其实是好 read**
  （Q<20 碱基比例 ≤10%、平均质量 ≥Q33）——它们只是**恰好超过 20% 低质阈值**。
- **误杀是系统性的，不是随机的**：被剔 read 平均长度 111 bp，保留的 121 bp。
  **长 read 天然更容易累积到 20% 的低质比例**（同样的 5 个低质碱基，在 120 bp 的
  read 上占 4.2%，在 100 bp 的 read 上占 5%，在 60 bp 的 read 上占 8.3%），
  所以阈值越严，越优先杀掉长 read —— 长 read 恰恰是比对信息量最大的。
- 净效果：`-u 20` 换到 Q30 +3.6 个点（87.3→90.9），代价是丢掉 15% 的数据量、
  且偏向性丢掉长 read。**对下游 DEG 分析，过度过滤 + 长度偏倚的风险可能大于
  Q30 提升带来的收益。**

**实操建议**：

- **Q30 ≥ 92% 的数据：不要加 `-u 20`**，用本流程现用配置即可（保留率 ~99.6%，
  低质 read 几乎为零）。加 `-u 20` 只会误杀。
- **Q30 85–92% 的数据：先跑一批看 `filtering_result.low_quality_reads`**。
  若这个数 < 总 read 的 1%，说明低质 read 本来就少，不用加 `-u`；
  若 > 5%，再加 `-u 20`，并**在报告里写明"低质量 read 过滤阈值 -u 20"**，
  同时报告过滤前后数据量（因为会偏向性丢长 read）。
- **Q30 < 85% 的数据：先查原因**（建库/测序/仪器/接头），不要靠调 fastp 参数补救。
  质控参数救不回坏数据。
- 想知道"多少低质 read 被剔"：`filtering_result.low_quality_reads`；
  想区分误杀：比对被剔 read 的长度分布 vs 保留 read 的长度分布（长度偏倚 = 有偏过滤）。

## Stage 2: HISAT2 Alignment

### GFF3 → GTF conversion (critical!)

HISAT2 `extract_splice_sites.py` and `extract_exons.py` require GTF format with `gene_id`/`transcript_id`. GFF3 (with `ID`/`Parent`) silently produces 0-byte output.

```bash
gffread annotation.gff3 -T -o annotation.gtf
```

### Build index

```bash
hisat2_extract_splice_sites.py annotation.gtf > splice_sites.txt
hisat2_extract_exons.py annotation.gtf > exons.txt
wc -l splice_sites.txt exons.txt  # Verify non-empty

hisat2-build -p 32 \
    --ss splice_sites.txt \
    --exon exons.txt \
    genome.fa hisat2_index
```

### Pitfall: `--no-unal` makes flagstat report 100% mapped — use the HISAT2 log rate instead

With `--no-unal` unmapped reads never enter the BAM, so `samtools flagstat` prints `mapped (100.00%)`
and `idxstats` shows `unmapped=0`. Reporting that as the alignment rate is wrong — quote the
`overall alignment rate` from the per-sample `.hisat2.log`. HISAT2 also never marks duplicates, so the
flagstat `duplicates` line is always 0; do not present it as a duplication metric for RNA-seq.

### Pitfall: sample IDs containing parentheses (or spaces) silently produce 0-read BAMs

Vendor sample names like `H4(3)`, `H5(2)`, `WT (1)` break HISAT2: the perl wrapper builds unquoted
`sh -c` command strings from the `-1`/`-2` paths, the shell aborts on the parenthesis, and the pipeline
still **exits 0**.

Symptoms (all together):
- `xx.hisat2.log` contains `sh: 1: Syntax error: "(" unexpected`
- then `0 reads` / `0.00% overall alignment rate`
- BAM is a few hundred bytes, yet `samtools index`/`flagstat` succeed → a naive `[ -s bam ]`
  skip-guard marks the sample done and the batch "completes" with empty files
- fastp, samtools, Python/R are all unaffected, so only HISAT2 fails — easy to miss in a batch

Fix: sanitize IDs to `[A-Za-z0-9._-]` for the whole analysis (`H4(3)` → `H4-3`) **before** running fastp,
and keep `lists/sample_id_mapping.tsv` with analysis_id ↔ vendor_name. Leave delivered raw filenames
untouched. Always add a post-mapping guard so a silent 0-read BAM can never count as done:

```bash
pairs=$(grep -oP '^\d+(?= reads; of these:)' "$LOGS/${s}.hisat2.log" || true)
if [ -z "${pairs:-}" ] || [ "${pairs:-0}" -eq 0 ] || [ "$(stat -c%s "$bam")" -lt 1000000 ]; then
    echo "[FAIL_EMPTY] $s"; rm -f "$bam"; return 1
fi
```

### Pitfall: NFS index corruption

On NFS >95% full, index files silently corrupt. Build on local `/tmp`:

```bash
TMPDIR=/tmp/hisat2_build_$$
mkdir -p $TMPDIR
hisat2-build -p 12 --ss splice_sites.txt --exon exons.txt genome.fa $TMPDIR/index
cp $TMPDIR/index.*.ht2 /target/nfs/dir/
rm -rf $TMPDIR
```

### Pitfall: hisat2-build fails with `/usr/bin/env: 'python': No such file or directory`

`hisat2-build` is a Python wrapper and needs `python` in PATH. On CNS nodes (and minimal/minimal-ssh shells), the default PATH often has only `/usr/bin/python3` (no bare `python`). Symptom: index build "finishes" instantly with empty .ht2 files.

```bash
# ❌ Fails on CNS1/CNS2 default PATH (hisat2-build can't find python)
hisat2-build -p 32 genome.fa index

# ✅ Fix: prepend miniforge bin so `python` resolves
export PATH=/media/nfs1/hermes/miniforge3/bin:$PATH
hisat2-build -p 32 genome.fa index
```

Verify index actually built before proceeding:
```bash
ls -la index.*.ht2          # should show 8 files, ~100MB+ each for large genomes
# Empty/small .ht2 (few bytes) + "Index complete" in log = python PATH failure, rebuild with PATH fix
```

### Pitfall: haplotyped-genome index build is slow (large)

Haplotype-resolved genomes (e.g. 29 chr named chr01a-14b, two haplotypes in one FASTA) produce ~2x larger index and slow builds. Expect several minutes to ~2h for ~600Mb-1Gb haplotyped genomes, building the SA-index stage. Check progress via `.sa` file count in `/tmp` (target ≈ index size/10MB), not just the log (which may not flush). The `.1.ht2`/`.2.ht2`/`.4.ht2` files are written early; `.sa` files accumulate during the suffix-array stage then get consolidated.

### Pitfall: Batch mapping fails en-masse when NFS disk fills up (exit 134 SIGABRT)

On large batch mapping (hundreds of samples), if the NFS output disk fills mid-run, a large fraction of samples fail with **`hisat2-align exited with value 134`** (SIGABRT) and log `Error while flushing and closing output` / `Aborted (core dumped)`. This happens during SAM output write. It is NOT an alignment/index problem — the successfully-completed samples (run before disk filled) have normal alignment rates.

Symptoms & diagnosis:
- Log shows many `FAIL_ALIGN: <sid>` entries; per-sample `.hisat2.log` ends with the SIGABRT error.
- `df -h` shows the NFS output disk at/near 100%.
- The failure pattern is **not random** — the first samples succeed, later ones all fail as disk fills.

Prevention & recovery:
1. **Before launching batch mapping, check disk headroom**: `df -h /media/nfs2` AND estimate BAM output size (BAM ≈ 1–1.5× clean FASTQ size). With ~0.6TB clean data → need ~1TB free for BAMs.
2. **Keep disk >10% free** — NFS at 95%+ is dangerous for write-heavy jobs.
3. **Recover by re-running only the pending samples** (those without a `.sorted.bam`): extract pending list from the sample manifest, write a `pending_samples.txt`, and point the retry script at it. The retry loop already skips samples whose BAM exists.
4. `SAMPLES` manifest can get out of sync if a `1.rawdata/` (symlink) dir is deleted and raw fastq cleaned — derive the manifest from `2.clean/*_clean_R1.fastq.gz` instead (the mapping input), not from rawdata.

### ⚠️ CRITICAL: exit 134 has TWO distinct root causes — NFS-full vs /tmp-full

`hisat2-align exited with value 134` (SIGABRT, `Error while flushing and closing output`) is produced by **BOTH** "output disk full" AND "temp /tmp disk full". Diagnose which before fixing, else the wrong fix fails again (this happened in production: NFS-full was fixed by moving SAM to /tmp, but then /tmp-full caused a SECOND mass failure).

**Same symptom, two causes, opposite fixes:**

| Cause | Signature | Fix |
|-------|-----------|-----|
| **NFS output disk full** | `df -h /media/nfs2` near 100%; first samples succeed then later ones fail as disk fills | Write BAM sort to `/tmp` (local SSD) then `cp` to NFS |
| **/tmp too small for parallel SAMs** | `df -h /tmp` low (e.g. 37GB); EVERY sample fails after running ~3-4 min (SAM files each 10-20GB, `-p` parallel × N samples exhausts /tmp mid-write) | Write SAM to **NFS** (when NFS has plenty free) instead of /tmp |

**Key diagnostic — don't blindly put SAM in /tmp.** In this fig project, a 619Mb haplotyped genome produced **SAM files of 10-21 GB per sample**. CNS1 `/tmp` (root partition) had only ~37GB free → parallel-3 hisat2 (each 16 threads) writing 10-20GB SAMs exhausted /tmp within ~4 min → every sample failed `exit 134`. Single-sample manual runs succeeded (21GB SAM wrote fine) because only one SAM in /tmp at a time.

**Check BOTH disks before choosing where SAM/BAM go:**
```bash
df -h /media/nfs2   # NFS output disk
df -h /tmp          # local scratch — check it actually has room for PARALLEL×SAMsize
```

**Correct decision table:**
- NFS has room (e.g. >1TB), /tmp is small → **write SAM to NFS** (`$MAPPING/tmp_sam/`), sort on NFS, `mv` to bam dir. (This is what actually worked for 138 pending samples.)
- NFS is full, /tmp has room → write to /tmp then cp to NFS (original guidance).
- **Never assume /tmp has space** — on multi-user shared nodes, `/tmp` root partition is often small and crowded by other users' files (and stale SAMs from prior runs can linger for tens of GB).

Also: a leftover multi-GB SAM from a previous failed/manual run in `/tmp` permanently shrinks the available scratch. Clean your own test SAMs (`rm /tmp/*_hisat2.sam`), but coordinate with the user — they may deny deletions of anything in `/tmp` since it's shared.

```bash
# Robust pattern that writes SAM on NFS (works when NFS has space, /tmp doesn't):
SAMTMP=$MAPPING/tmp_sam
mkdir -p "$SAMTMP" "$BAMDIR"
hisat2 -p 16 -x $INDEX -1 $r1 -2 $r2 -S $SAMTMP/${sid}.sam 2> logs/${sid}.hisat2.log
samtools sort -@ 4 -m 4G -o $SAMTMP/${sid}.sorted.bam $SAMTMP/${sid}.sam
rm $SAMTMP/${sid}.sam
mv $SAMTMP/${sid}.sorted.bam $BAMDIR/${sid}.sorted.bam
samtools index -c $BAMDIR/${sid}.sorted.bam
```

### Pitfall: NFS pipe SIGPIPE (exit 141)

`hisat2 ... | samtools sort ...` pipe can fail with SIGPIPE on NFS. Fix: write SAM to temp file first, then sort.

```bash
# ❌ Risky on NFS
hisat2 -p 16 -x index -1 R1.fq -2 R2.fq | samtools sort -@ 4 -o sample.bam -

# ✅ Safe: two-step
TMP_SAM=/tmp/${sample}.sam
hisat2 -p 16 -x index -1 R1.fq -2 R2.fq -S "$TMP_SAM"
samtools sort -@ 4 -o sample.sorted.bam "$TMP_SAM"
rm -f "$TMP_SAM"
```

### Pitfall: /tmp exhaustion → hisat2 SIGABRT (exit 134) under parallel

**Symptom:** hisat2 fails with `Error while flushing and closing output` / `terminate called after throwing an instance of 'int'` / `Aborted (core dumped)` / `hisat2-align exited with value 134`. A **single** manual run succeeds, but when running N samples in parallel they ALL fail after ~4 min.

**Root cause:** each hisat2 writes a 10-20GB SAM (RNA-seq ~2-5Gb of reads → 10-21GB SAM). Running 4+ in parallel each writing SAM to `/tmp` fills `/tmp` (shared machine, often only 30-50GB free after other work) → output flush fails → SIGABRT. **NOT an NFS problem** — it's /tmp space.

**Diagnosis:** check `/tmp` free space (`df -h /tmp`), and confirm single-sample run succeeds while parallel fails. Also note leftover test SAMs from earlier manual runs can eat `/tmp` (check `du -sh /tmp/*.sam`).

**Fix:** write the SAM (and sorted BAM) temp files to NFS instead of `/tmp` — NFS typically has far more space. Parallel 3-4 with SAM on NFS works reliably:

```bash
# ✅ Safe under parallel: SAM → NFS tmp dir, not /tmp
SAMTMP=/media/nfs2/project/3.Mapping/tmp_sam
mkdir -p "$SAMTMP"
hisat2 -p 16 -x index -1 R1.fq -2 R2.fq -S "$SAMTMP/${sid}.sam"
samtools sort -@ 4 -o "$SAMTMP/${sid}_sorted.bam" "$SAMTMP/${sid}.sam"
mv "$SAMTMP/${sid}_sorted.bam" bam/${sid}.sorted.bam
samtools index -c bam/${sid}.sorted.bam
rm -f "$SAMTMP/${sid}.sam"
```

**Verification of a completed mapping run:** check the run log for FAIL_ALIGN/FAIL_SORT counts (not just `FAIL:` — the script logs `FAIL_ALIGN:`/`FAIL_SORT:`), and confirm every sample has both a `.sorted.bam` AND `.csi` index.

### Pitfall: `hisat2 -x` takes the INDEX PREFIX, not the fasta — error message is misleading

Passing the fasta path to `-x` fails with `(ERR): ".../genome.fasta" does not exist` — even though the fasta DOES exist. hisat2 is looking for `genome.fasta.1.ht2` etc. (index files with that string as prefix). Symptom: both parallel samples fail instantly, empty stderr logs otherwise.

```bash
hisat2 -x $REF/ref/genome.fasta ...   # ❌ looks for genome.fasta.1.ht2 — never exists
hisat2 -x $REF/ref/genome ...         # ✅ index prefix (files genome.1.ht2 ... genome.8.ht2)
```

Verify before batch: `ls ref/genome.{1..8}.ht2` all present and nonzero (`.3.ht2`/`.7.ht2`/`.8.ht2` are legitimately tiny — 12-152 bytes).

### Recipe: LOW-IO batch mapping (user-requested config, verified on wuhuaguo 197 samples)

When the user asks to keep IO low (shared NFS, other jobs running), use this instead of the 16-thread × parallel-3 pattern above:

```bash
map_one() {
  local S=$1
  [ -s "$BAM/${S}.bam" ] && { echo "[skip] $S"; return; }   # resumable
  hisat2 -p 4 --dta --no-unal -x $REF \
    -1 ${CLEAN}/${S}_clean_R1.fastq.gz -2 ${CLEAN}/${S}_clean_R2.fastq.gz \
    2> $LOGS/${S}.hisat2.log \
  | samtools sort -@ 2 -m 3G -T /tmp/whgsort_${S} -o $BAM/${S}.bam -
  samtools index -@ 2 $BAM/${S}.bam
  samtools flagstat $BAM/${S}.bam > $SUMMARY/${S}.flagstat
}
export -f map_one
cat all_samples.txt | xargs -P 2 -I{} bash -c 'map_one {}'
```

Key choices:
- **2 samples × hisat2 -p 4 = 8 threads total** — gentle on NFS. ~590Gb clean data → ~25-30h.
- **Pipe straight into sort (no SAM at all)** — IO-lightest. `-T /tmp/...` puts sort temp shards on LOCAL disk, NFS only receives the final BAM. (If exit 141 SIGPIPE appears on a given NFS, fall back to the SAM-on-NFS two-step pattern above.)
- `--no-unal` cuts unmapped reads out of the BAM → smaller output, faster.
- **Script must start with `set -e -o pipefail`** — without pipefail, a failed hisat2 inside the pipe still yields an empty "successful" BAM and the skip-if-exists guard will treat it as done on re-run.
- Skip guard `[ -s bam ]` makes the script resumable after interruption.
- Do NOT trust file timestamps for progress checks — the cluster clock can be off by weeks; count BAM files instead.

### Killing a batch mapping pipeline cleanly

Killing `xargs` does NOT kill the running children: the perl hisat2 wrappers (and `hisat2-align-s`) get **reparented to PID 1** and keep burning CPU. Also `pkill -f "pattern"` can match your own bash -c wrapper and kill itself (exit -15). Procedure:

```bash
ps -ef | grep -E "hisat2|map_one|run_mapping|xargs" | grep -v grep   # list real PIDs
kill -9 <pid1> <pid2> ...                                            # explicit PIDs, NOT kill $(pgrep -f ...) — user denies that form
# repeat the ps check until count = 0; orphaned hisat2 may surface one at a time
rm -f $BAM/*.bam $SUMMARY/*.flagstat   # discard partial outputs from the failed run
```

### Pitfall: CSI index for large chromosomes

BAI format limit = 2^29 = 536,870,912 bp. For chromosomes >512Mb, use CSI:

```bash
samtools index -c sample.sorted.bam    # CSI (no limit)
# vs
samtools index sample.sorted.bam       # BAI (limit 512Mb)
```

### Alignment command

```bash
hisat2 -p 16 -x ref/hisat2_index \
    -1 sample_clean_R1.fastq.gz -2 sample_clean_R2.fastq.gz \
    -S /tmp/${sample}.sam 2> logs/${sample}.hisat2.log

samtools sort -@ 4 -m 4G -o bam/${sample}.sorted.bam /tmp/${sample}.sam
samtools index -c bam/${sample}.sorted.bam
rm -f /tmp/${sample}.sam
```

### Key parameters

| Parameter | Value | Why |
|-----------|-------|-----|
| `--dta` | NOT used | For featureCounts/DESeq2 quantification, not needed. Saves ~10GB memory. |
| `--rna-strandness` | Omit (unstranded) | Use `FR`/`RF` for stranded libraries |
| Threads | 16 | 16 threads × ~3GB each = ~48GB for large genomes |
| `samtools sort -m` | 4G | Per-thread memory for sort. 4 threads × 4G = 16G max |

## Stage 3: featureCounts Quantification

### Pitfall: GFF3 gene_id attribute

featureCounts needs `gene_id` in GTF annotation. Use gffread-converted GTF:

```bash
# Requires GTF (not raw GFF3)
featureCounts -T 16 \
    -a annotation.gtf \
    -g gene_id \
    -t exon \
    -s 0 \
    -p --countReadPairs \
    -o counts_raw.txt \
    sample1.sorted.bam sample2.sorted.bam ...
```

### Error message if using raw GFF3

```
ERROR: failed to find the gene identifier attribute in the 9th column of the provided GTF file.
The specified gene identifier attribute is 'gene_id'
```

## Stage 4: TPM Calculation + Tissue-Level Filtering

### TPM formula

```
TPM_i = (count_i / length_i_kb) * 1e6 / sum(count_j / length_j_kb for all j)
```

### Filter rule: keep if ANY tissue has mean TPM >= 1

```python
# Define tissue from sample name
def get_tissue(sample_name):
    # e.g. Gli_root1_rep1 → root1
    return sample_name.replace('Gli_', '').split('_rep')[0]

# Group samples by tissue
tissue_groups = {}
for s in sample_names:
    t = get_tissue(s)
    tissue_groups.setdefault(t, []).append(s)

# Keep gene if any tissue mean TPM >= 1
filtered_genes = []
for gene_id in tpm:
    keep = False
    for tissue, tissue_samples in tissue_groups.items():
        indices = [sample_names.index(s) for s in tissue_samples]
        mean_tpm = sum(tpm[gene_id][i] for i in indices) / len(indices)
        if mean_tpm >= 1:
            keep = True
            break
    if keep:
        filtered_genes.append(gene_id)
```

### WRONG approach (do NOT use)

```python
# ❌ Global mean across all samples — removes tissue-specific genes
mean_tpm = sum(tpm[gene_id]) / len(samples)
```

### Output files

```
counts_matrix.tsv                — raw counts (all genes)
tpm_matrix.tsv                   — TPM values (all genes)
counts_matrix_filtered.tsv       — raw counts (tissue-level TPM >= 1)
tpm_matrix_filtered.tsv          — TPM values (tissue-level TPM >= 1)
```

### Complete Python script

```python
#!/usr/bin/env python3
"""Calculate TPM from featureCounts, filter by tissue-level mean TPM >= 1."""
import os

counts_file = "counts_raw.txt"
out_dir = "."

def get_tissue(sample_name):
    """Extract tissue from sample name. Customize for your naming."""
    return sample_name.replace('Gli_', '').split('_rep')[0]

with open(counts_file) as f:
    lines = [l for l in f if not l.startswith('#')]

header = lines[0].strip().split('\t')
sample_names = [os.path.basename(s).replace('.sorted.bam', '') for s in header[6:]]

gene_data = {}
gene_lengths = {}
for line in lines[1:]:
    cols = line.strip().split('\t')
    gene_id = cols[0]
    gene_length = int(cols[5])
    counts = [int(c) for c in cols[6:]]
    gene_data[gene_id] = counts
    gene_lengths[gene_id] = gene_length

# Calculate TPM
tpm = {}
for gene_id in gene_data:
    rpk = [c / (gene_lengths[gene_id] / 1000.0) for c in gene_data[gene_id]]
    tpm[gene_id] = rpk

for j in range(len(sample_names)):
    total_rpk = sum(tpm[g][j] for g in tpm)
    scale = 1e6 / total_rpk if total_rpk > 0 else 0
    for gene_id in tpm:
        tpm[gene_id][j] *= scale

# Group by tissue
sample_to_tissue = {s: get_tissue(s) for s in sample_names}
tissue_groups = {}
for s in sample_names:
    t = sample_to_tissue[s]
    tissue_groups.setdefault(t, []).append(s)

# Filter
filtered = [g for g in sorted(tpm.keys())
    if any(sum(tpm[g][sample_names.index(s)] for s in ts) / len(ts) >= 1
           for ts in tissue_groups.values())]

# Write matrices
for data, suffix in [(gene_data, "counts_matrix"), (tpm, "tpm_matrix")]:
    with open(os.path.join(out_dir, f"{suffix}.tsv"), 'w') as f:
        f.write("GeneID\t" + "\t".join(sample_names) + "\n")
        for gene_id in sorted(data.keys()):
            f.write(gene_id + "\t" + "\t".join(
                f"{data[gene_id][j]:.4f}" if 'tpm' in suffix else str(data[gene_id][j])
                for j in range(len(sample_names))) + "\n")
    print(f"Written: {suffix}.tsv")

for data, suffix in [(gene_data, "counts_matrix_filtered"), (tpm, "tpm_matrix_filtered")]:
    with open(os.path.join(out_dir, f"{suffix}.tsv"), 'w') as f:
        f.write("GeneID\t" + "\t".join(sample_names) + "\n")
        for gene_id in filtered:
            f.write(gene_id + "\t" + "\t".join(
                f"{data[gene_id][j]:.4f}" if 'tpm' in suffix else str(data[gene_id][j])
                for j in range(len(sample_names))) + "\n")
    print(f"Written: {suffix}.tsv")

print(f"Total: {len(gene_data)} -> Filtered: {len(filtered)} (removed {len(gene_data) - len(filtered)})")
```

## Pitfall: Gene IDs in TPM/count matrices carry transcript/version suffixes

When the user asks to plot a heatmap / extract rows for a **specific gene list** from a TPM or count matrix, the matrix row IDs often carry a transcript/version suffix (e.g. `Seita.6G138300.1` for foxtail millet / guzi, or `X.mRNA1`), while the user gives bare gene IDs (`Seita.6G138300`). This breaks exact matching and silently returns "missing genes".

Two traps:
1. **`grep -c "^$gene" matrix` is a false-positive sanity check** — `^Seita.6G138300` prefix-matches `Seita.6G138300.1`, so it looks like the gene exists when an exact R match fails. Don't rely on it.
2. In R, `read.table` row names keep the full `...1` ID, so `setdiff(target, rownames(tpm))` reports all genes "missing".

**Fix — strip the version suffix from row names before matching:**
```r
tpm <- read.table(tpm_path, header=TRUE, row.names=1, check.names=FALSE, sep="\t")
baseid <- sub("[.][0-9]+$", "", rownames(tpm))   # Seita.6G138300.1 -> Seita.6G138300
rownames(tpm) <- baseid
sub <- as.matrix(tpm[target_genes, , drop=FALSE])  # drop=FALSE keeps 1-gene case as matrix
```

Use `[.]` in the regex (not `\\.`) to avoid R escape hell when the script is written via write_file/patch tools (see the R regex pitfall below).

**Also note:** when computing column means over a matrix, `as.matrix()` the extracted sub-matrix first (a `data.frame` keeps the `drop=FALSE` dims but column-mean with character indices returns NA on non-numeric). Verified on guzi project: `05.TPM/genes.TPM.EXPR.matrix` has 72 columns (36 CK control + 36 `*Y` treated), sample grouping in `06.DEG/stats/sample_info.tsv` (treatment/genotype/group), 5 target genes all mapped to `.1` isoforms.

**Heatmap recipe that worked** (user's guzi RNA-seq, 5 genes × 72 samples): `pheatmap` on `log2(TPM+1)` and on row z-score (`t(scale(t(logmat)))`), with `annotation_col` = Treatment(Control/Treatment), Genotype, Group(tolerant/sensitive) from sample_info; color blue-white-red `#2166AC`→`#F7F7F7`→`#B2182B`; output PDF + 300 DPI PNG to `figures/gene_heatmap/`. pheatmap PNG output is rejected by the vision API (400 code 10305) — verify the data programmatically instead of by image inspection.

## Stage 5: DESeq2 Differential Expression

### Prerequisites

```bash
# Use deseq2 conda environment
/media/nfs1/hermes/miniforge3/envs/deseq2/bin/Rscript run_deseq2.R
```

### R script

```r
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Read count matrix
counts <- read.table("counts_matrix_filtered.tsv", header=TRUE, row.names=1, check.names=FALSE)

# Build metadata
samples <- colnames(counts)
condition <- gsub("_rep\\d+$", "", samples)  # e.g. Gli_root1_rep1 → Gli_root1
coldata <- data.frame(row.names=samples, condition=factor(condition))
coldata$condition <- relevel(coldata$condition, ref="root1")  # Set baseline

# DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)

# Results for each comparison
res <- results(dds, contrast=c("condition", "root2", "root1"))
resLFC <- lfcShrink(dds, contrast=c("condition", "root2", "root1"), type="ashr")
# If ashr fails: use unshrunken LFC (padj/pvalue unchanged)
```

### When lfcShrink with ashr fails

```r
res <- results(dds, contrast=c("condition", "root2", "root1"))
resLFC <- tryCatch({
    lfcShrink(dds, contrast=c("condition", "root2", "root1"), type="ashr")
}, error = function(e) {
    message("ashr unavailable, using unshrunken LFC")
    res
})
```

### Output files

```
DEG/DEG_summary.csv
DEG/DEG_root2_vs_root1.csv
DEG/volcano_root2_vs_root1.png
DEG/MA_root2_vs_root1.png
DEG/PCA_all_samples.png
DEG/heatmap_top_DEGs.png
```

## Stage 6: Functional Enrichment

### Prerequisites: eggNOG-mapper annotation

Run eggNOG-mapper on proteins to get GO/KEGG annotations:

```bash
emapper.py \
    -i proteins.fa \
    -o output_prefix \
    --data_dir /media/nfs1/hermes/db/eggnog \
    --cpu 40 \
    --tax_scope 33090 \
    --go_evidence non-electronic \
    --pfam_realign none \
    --temp_dir /tmp \
    --override
```

### Build OrgDb from eggNOG results

For non-model plants without org.*.eg.db, use `enricher()` with custom term2gene mapping:

```r
library(clusterProfiler)

# Parse eggNOG results
emapper <- read.delim("output_prefix.emapper.annotations", comment.char="#")
go_annot <- emapper[, c("query", "GOs")]
go_annot <- go_annot[go_annot$GOs != "" & go_annot$GOs != "-", ]

# Split multi-GO entries
go_list <- strsplit(go_annot$GOs, ",")
term2gene <- data.frame(
    GO = unlist(go_list),
    Gene = rep(go_annot$query, lengths(go_list))
)
term2name <- go_annot  # GO ID → Description mapping from go-basic.obo

# Run enrichment
ego <- enricher(
    gene = sig_genes,
    universe = background_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)
```

### Or build OrgDb from eggNOG annotations

For non-model plants without org.*.eg.db, use `enricher()` with custom term2gene mapping. **Critical: TERM2GENE column names must be uppercase `TERM` and `GENE`.**

### Pitfall: enricher() TERM2GENE column names

`clusterProfiler::enricher()` requires the TERM2GENE data frame to have EXACTLY two columns named `TERM` and `GENE` (uppercase). Any other names (e.g. `GO`, `Gene`, `KO`, `term`, `gene`) will silently produce 0 results — `enricher()` returns an empty object.

```r
# ❌ WRONG — column names don't match, enrichment returns 0 results
term2gene <- data.frame(
  GO = go_terms,    # should be "TERM"
  Gene = go_genes,  # should be "GENE"
  stringsAsFactors = FALSE
)
ego <- enricher(gene = sig_genes, TERM2GENE = term2gene)  # Returns 0 terms

# ✅ CORRECT — column names must be TERM and GENE
term2gene <- data.frame(
  TERM = go_terms,
  GENE = go_genes,
  stringsAsFactors = FALSE
)
ego <- enricher(gene = sig_genes, TERM2GENE = term2gene)  # Works correctly
```

Also, all downstream `intersect()` calls must use `$GENE` not `$Gene`:
```r
sig_go <- intersect(sig_genes, term2gene$GENE)   # ✅ correct
sig_go <- intersect(sig_genes, term2gene$Gene)    # ❌ NULL, returns 0
```

### Pitfall: eggNOG annotation parsing

eggNOG output files have multiple `#` comment lines. `read.delim(comment.char="#", header=TRUE)` will use the first data row as column names instead of the actual header. Always use:

```r
# ❌ WRONG — header=TRUE with comment.char causes column name corruption
ann <- read.delim("file.emapper.annotations", comment.char="#", header=TRUE, sep="\t")

# ✅ CORRECT — specify column names explicitly
col_names <- c("query","seed_ortholog","evalue","score","eggNOG_OGs",
               "max_annot_lvl","COG_category","Description","Preferred_name",
               "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
               "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy",
               "BiGG_Reaction","PFAMs")
ann <- read.delim("file.emapper.annotations",
                  comment.char="#", header=FALSE, col.names=col_names,
                  sep="\t", stringsAsFactors=FALSE)
```

### Pitfall: R regex escape in write_file scripts

When writing R scripts via `write_file` or `patch` tools, R string escaping is tricky. `sub("\\\\.mRNA", ...)` in R source means `sub("\\.mRNA, ...)` which matches literal `.mRNA`. But the tool's write pipeline may add extra backslashes.

**Best practice:** Use `[.]` (character class) instead of `\\.` to avoid escape issues:

```r
# ❌ Escape hell — backslash count varies between write_file, patch, and R
ann$gene <- sub("\\\\.mRNA[0-9]+$", "", ann$query)

# ✅ Safe — character class [.] matches literal dot, no escaping needed
ann$gene <- sub("[.]mRNA[0-9]+$", "", ann$query)
```

### Build term2gene from eggNOG annotations

```r
library(clusterProfiler)

# Read eggNOG annotations (skip comment lines, specify column names manually)
col_names <- c("query","seed_ortholog","evalue","score","eggNOG_OGs",
               "max_annot_lvl","COG_category","Description","Preferred_name",
               "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
               "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy",
               "BiGG_Reaction","PFAMs")

ann <- read.delim("output_prefix.emapper.annotations",
                  comment.char = "#", header = FALSE, col.names = col_names,
                  sep = "\t", stringsAsFactors = FALSE)

# Remove transcript suffix (e.g. .mRNA1) to match DESeq2 gene IDs
ann$gene <- sub("\\\\.mRNA[0-9]+$", "", ann$query)  # Glh1_01G000001.mRNA1 → Glh1_01G000001

# Build GO term2gene
go_ann <- ann[ann$GOs != "" & ann$GOs != "-", c("gene", "GOs")]
go_list <- strsplit(go_ann$GOs, ",")
term2gene_go <- data.frame(
  TERM = unlist(go_list),
  GENE = rep(go_ann$gene, lengths(go_list)),
  stringsAsFactors = FALSE
)

# Build KEGG Pathway term2gene (NOT KO level — see preferred enrichment below)
kegg_p_ann <- ann[ann$KEGG_Pathway != "" & ann$KEGG_Pathway != "-", c("gene", "KEGG_Pathway")]
kp_list  <- strsplit(kegg_p_ann$KEGG_Pathway, ",")
kp_terms <- unlist(kp_list); kp_genes <- rep(kegg_p_ann$gene, lengths(kp_list))
keep <- grepl("^ko[0-9]{5}$", kp_terms)   # ko##### only, drop map##### dups
term2gene_ko <- data.frame(TERM = kp_terms[keep], GENE = kp_genes[keep], stringsAsFactors = FALSE)
# Plant filter: exclude animal-specific pathway prefixes
term2gene_ko <- term2gene_ko[!grepl("^(ko05|ko047|ko049|ko046|ko043)", term2gene_ko$TERM), ]

# Run enrichment
ego <- enricher(
  gene = sig_genes,
  universe = bg_genes,
  TERM2GENE = term2gene_go,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500
)
```

### ✅ PREFERRED enrichment for multi-group comparison (user decided 2026-08)

The user's standard for DEG enrichment figures is the **multi-group dotplot with x = group** (NOT per-group facet/single plots). When you have multiple comparisons/time-points (e.g. root2/3/4/6 vs root1), bind all groups' significant results with a `Group` column and plot with x=Group, y=term, size=Count, color=`-log10(p.adjust)`. Full recipe is in `bio-pathway-enrichment` and the TO-GCN skill.

**Key improvements over the older KO-level/GO-ID approach** (this is what the user wanted changed — the original `run_enrichment.R` produced unreadable figures):
1. **GO: pass TERM2NAME** so Description shows readable names, not GO IDs. Get names from GO.db: `AnnotationDbi::select(GO.db, keys=go_ids, columns=c("TERM","ONTOLOGY"), keytype="GOID")` (the `TERM` column = description). **Plot ONLY the BP subontology** (filter results by `ONTOLOGY=="BP"`).
2. **KEGG: use `KEGG_Pathway` column (ko##### pathway level), NOT `KEGG_ko` (K-number level).** KO-level gives few (5-18) unreadable `K###` labels; Pathway-level shows readable names ("Photosynthesis", "Ribosome", "MAPK signaling pathway - plant"). Apply the **plant filter** (exclude prefixes `ko05/ko047/ko049/ko046/ko043` — animal pathways). Term2name from the reusable global DB `/media/nfs1/hermes/db/kegg_pathway_names.tsv` (col1=TERM, col2=NAME; 1170 rows). Pass `TERM2NAME` to enricher.
3. **Background** = all genes with non-NA pvalue in that comparison's DESeq2 result (`deg$gene[!is.na(deg$pvalue)]`), intersected with each term2gene's genes.
4. Plot color blue→red `#4575B4`→`#D73027` (low→high `-log10(p.adjust)`), size=Count, top-N (e.g. 10) terms per group merged & deduped, y-order by min p.adjust across all groups, save 8×7 PDF+PNG (600 dpi). Handle `-log10(0)=Inf` (cap at max finite + 1).

Working reference script: `RNA-seq/6.Enrichment/multigroup/scripts/run_multigroup_enrich.R` (珊瑚菜 project). Output dir `6.Enrichment/multigroup/{results,figures}`.

**珊瑚菜 DEG enrichment result (validates biology; DEG counts from DESeq2 padj<0.05):**
- DEG counts vs root1: root2=4,434 (up 2,620 / down 1,814), root3=1,979 (1,296/683), root4=6,713 (3,737/2,976), root6=6,674 (3,833/2,841).
- GO BP sig terms: root2=48 / root3=54 / root4=51 / root6=75 (73 unique). Key: response to chitin (all time points, most significant), water/light/cold stress responses (root2-4), defense responses to fungus/bacterium (root6).
- KEGG Pathway sig: root2=11 / root3=12 / root4=9 / root6=18 (45 unique). Key: plant hormone signal transduction + MAPK (root2), phenylpropanoid biosynthesis (root3+), phenylpropanoid + flavonoid biosynthesis (root4), alpha-linolenic acid metabolism (JA precursor) + hormone signaling + phenylpropanoid (root6).
- Biological trend: defense/environment responses active throughout; secondary metabolism (phenylpropanoid/flavonoid) intensifies in later stages; mature root (root6) dominated by biotic defense + secondary metabolism — consistent with medicinal-component accumulation in mature roots.
- Report generated for this: `RNA-seq/00.report/report_14sample_5stage_deg_shanhucai.html` + `_bundle.tar.gz`.

## Directory Structure

```
project/
├── 1.rawdata/
│   └── *.fastq.gz
├── 2.clean/
│   ├── *.fastq.gz
│   ├── logs/
│   └── run_fastp.sh
├── 3.Mapping/
│   ├── ref/
│   │   ├── genome.fasta
│   │   └── G.li.hap1.*.ht2
│   ├── bam/
│   │   ├── *.sorted.bam
│   │   └── *.sorted.bam.bai (or .csi)
│   ├── summary/
│   ├── logs/
│   └── scripts/
│       └── run_hisat2.sh
├── 4.Quantification/
│   └── featureCounts/
│       ├── counts_raw.txt
│       ├── counts_matrix.tsv
│       ├── tpm_matrix.tsv
│       ├── counts_matrix_filtered.tsv
│       └── tpm_matrix_filtered.tsv
├── 5.DEG/
│   ├── run_deseq2.R
│   ├── DEG_summary.csv
│   ├── DEG_*_vs_root1.csv
│   ├── volcano_*.png
│   ├── MA_*.png
│   ├── PCA_all_samples.png
│   └── heatmap_top_DEGs.png
└── 6.Enrichment/
    └── ...
```

## Related Skills

This is the **primary end-to-end bulk RNA-seq pipeline** (HISAT2 + featureCounts + DESeq2 + enrichment). For alternative alignment/quantification strategies and specialized sub-analyses, use the sibling skills — each covers one workflow in detail:

- `bio-read-qc-fastp-workflow` — detailed fastp options
- `bio-read-alignment-hisat2-alignment` — comprehensive HISAT2 reference
- `bio-read-alignment-star-alignment` — STAR alignment (alternative to HISAT2)
- `bio-rna-quantification` — **Salmon/kallisto** pseudo-alignment + tximport + count QC (alternative to featureCounts)
- `bio-differential-expression` — DESeq2/edgeR detailed reference
- `bio-differential-expression-timeseries-de` — time-series / longitudinal RNA-seq
- `bio-alternative-splicing` — splicing analysis (rMATS, leafcutter, MAJIQ, SUPPA2)
- `bio-small-rna-seq` — miRNA / small RNA pipeline
- `bio-pathway-enrichment` — clusterProfiler enrichment detailed reference
- `bioinformatics/go-enrichment-pipeline` — GO/KEGG enrichment with 7 plot types
- `bioinformatics/functional-annotation-pipeline` — eggNOG-mapper functional annotation
- `bio-differential-expression-batch-correction` — ComBat/SVA batch correction

## Citations

### Tools

- fastp: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 2018, 34(17):i884–i890. doi:10.1093/bioinformatics/bty560
- HISAT2: Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 2019, 37:907–915. doi:10.1038/s41587-019-0201-4
- featureCounts: Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 2014, 30(7):923–930. doi:10.1093/bioinformatics/btt656
- DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 2014, 15:550. doi:10.1186/s13059-014-0550-8
- clusterProfiler: Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2021, 2(3):100141. doi:10.1016/j.xinn.2021.100141
- SAMtools: Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience, 2021, 10(2):giab008. doi:10.1093/gigascience/giab008
- gffread (GFF->GTF): Pertea M, et al. GFF Utilities: GffRead and GffCompare. F1000Research, 2020, 9:304. doi:10.12688/f1000research.23297.1
- eggNOG-mapper: Cantalapiedra CP, et al. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 2021, 38(12):5825-5829. doi:10.1093/molbev/msab293
- eggNOG 5.0 database: Huerta-Cepas J, et al. Nucleic Acids Research, 2019, 47(D1):D309-D314. doi:10.1093/nar/gky1085

### Methods template

> RNA-seq reads were quality-filtered using fastp (Chen et al., 2018) with default parameters. Clean reads were aligned to the reference genome using HISAT2 (Kim et al., 2019). Gene-level quantification was performed with featureCounts v2.1.1 (Liao et al., 2014) using the GTF annotation with parameters `-g gene_id -t exon -s 0 -p --countReadPairs`. Transcripts per million (TPM) were calculated from raw counts normalized by gene length and library size; genes with mean TPM < 1 across all samples within each tissue group were filtered as low-expressed. Differential expression analysis was performed using DESeq2 (Love et al., 2014) with the design ~ condition, treating root1 as the baseline. Genes with adjusted p-value < 0.05 were considered significantly differentially expressed. Functional enrichment analysis was performed using clusterProfiler (Wu et al., 2021) with eggNOG-mapper-derived GO and KEGG annotations.

### BibTeX（可直接粘贴）

```bibtex
@article{chen2018fastp,
  title={fastp: an ultra-fast all-in-one FASTQ preprocessor},
  author={Chen, Shubin and Zhou, Yanfang and Chen, Ying and Gu, Jianhui},
  journal={Bioinformatics}, volume={34}, number={17}, pages={i884--i890},
  year={2018}, doi={10.1093/bioinformatics/bty560}
}
@article{kim2019hisat2,
  title={Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype},
  author={Kim, Daehwan and Paggi, Giovanni M. and Park, Chul and Bennett, Cynthia and
          Salzberg, Steven L.},
  journal={Nature Biotechnology}, volume={37}, pages={907--915}, year={2019},
  doi={10.1038/s41587-019-0201-4}
}
@article{pertea2020gffread,
  title={GFF Utilities: GffRead and GffCompare},
  author={Pertea, Mihaela and Pohl, AW and Jovic, N and Gifford, M and
          Berard, S and Flicek, P and Gordon, L and Graves, S and
          Schmutz, J and Kim, D and Heng, H},
  journal={F1000Research}, volume={9}, pages={304}, year={2020},
  doi={10.12688/f1000research.23297.1}
}
@article{liao2014featurecounts,
  title={featureCounts: an efficient general purpose program for assigning sequence
         reads to genomic features},
  author={Liao, Yifang and Smyth, Gordon K. and Shi, Wei},
  journal={Bioinformatics}, volume={30}, number={7}, pages={923--930}, year={2014},
  doi={10.1093/bioinformatics/btt656}
}
@article{danecek2021samtools,
  title={Twelve years of SAMtools and BCFtools},
  author={Danecek, Petr and Abeysinghe, Sudruga and Twardziok, Petr and
          Li, Heng and Marth, Gabor and Abecasis, Goncal and Durbin, Alistair},
  journal={GigaScience}, volume={10}, number={2}, pages={giab008}, year={2021},
  doi={10.1093/gigascience/giab008}
}
@article{love2014deseq2,
  title={Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
  author={Love, Michael I. and Huber, Wolfgang and Anders, Simon},
  journal={Genome Biology}, volume={15}, pages={550}, year={2014},
  doi={10.1186/s13059-014-0550-8}
}
@article{wu2021clusterprofiler,
  title={clusterProfiler 4.0: A universal enrichment tool for interpreting omics data},
  author={Wu, Tianyi and Hu, Erjie and Xu, Shuai and others},
  journal={The Innovation}, volume={2}, number={3}, pages={100141}, year={2021},
  doi={10.1016/j.xinn.2021.100141}
}
@article{cantalapiedra2021eggnog,
  title={eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and
         Domain Prediction at the Metagenomic Scale},
  author={Cantalapiedra, Carlos P. and Hernandez-Plaza, Ana and Letunic, Ivica and
          Bork, Peer and Huerta-Cepas, Jaime},
  journal={Molecular Biology and Evolution}, volume={38}, number={12},
  pages={5825--5829}, year={2021}, doi={10.1093/molbev/msab293}
}
@article{huertacepas2019eggnog5,
  title={eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
         orthology resource based on 5090 organisms and 2502 viruses},
  author={Huerta-Cepas, Jaime and Szklarczyk, Damian and Heller, Davide and
          Hernandez-Plaza, Ana and Forslund, Sofia K. and Cook, Helen and
          Mende, Daniel R. and Letunic, Ivica and Rattei, Thomas and
          Jensen, Lars J. and von Mering, Christian and Bork, Peer},
  journal={Nucleic Acids Research}, volume={47}, number={D1}, pages={D309--D314},
  year={2019}, doi={10.1093/nar/gky1085}
}
@article{kim2013stringtie,
  title={StringTie: a transcriptome assembler and its applications},
  author={Kim, Daehwan and Pertea, Mihaela and Trapnell, Cole and Salzberg, Steven L.},
  journal={Nature Biotechnology}, volume={31}, pages={940--942}, year={2013},
  doi={10.1038/nbt.1196},
  note={cited only if transcripts were assembled with StringTie}
}
```

### Methods 模板句（中文，对应上文英文版）

> 原始 reads 经 fastp（Chen et al., 2018）质控后，使用 HISAT2（Kim et al., 2019）比对至参考基因组。
> 基因水平定量采用 featureCounts v2.1.1（Liao et al., 2014），参数为 `-g gene_id -t exon -s 0 -p --countReadPairs`。
> TPM 由原始 counts 按基因长度与文库大小归一化得到；按组织分组过滤低表达基因（任一组织平均 TPM >= 1 则保留）。
> 差异表达用 DESeq2（Love et al., 2014）分析，设计为 `~ condition`，以 root1 为基准，padj < 0.05 视为显著差异表达。
> 功能富集用 clusterProfiler（Wu et al., 2021），注释来源为 eggNOG-mapper 导出的 GO 与 KEGG。

## 中文使用说明

### 本 skill 适用场景

植物批量 RNA-seq 分析的标准流程，从原始数据到差异表达基因和功能富集，一步到位。

### 流程步骤

1. **fastp 过滤**：双端数据质量过滤，自动检测接头
2. **HISAT2 比对**：建索引 → 比对到参考基因组
3. **featureCounts 定量**：基因水平表达量计数
4. **TPM 计算 + 过滤**：按组织级别过滤低表达基因（任一组织平均TPM≥1则保留）
5. **DESeq2 差异表达**：以root1为基准，比较其他组织
6. **富集分析**：eggNOG功能注释 → clusterProfiler GO/KEGG 富集

### 关键坑点

- **GFF3→GTF 转换**：HISAT2 和 featureCounts 都要求 GTF 格式，必须先用 gffread 转换
- **NFS 管道失败**：`hisat2 | samtools sort` 在 NFS 上可能 SIGPIPE，先写 SAM 再 sort
- **磁盘满批量失败（exit 134 两种原因）**：大批量 mapping 前必须 `df -h /media/nfs2` 确认磁盘余量（BAM≈1-1.5×clean大小）；NFS≥95%时写输出 SIGABRT（exit 134）。**但 /tmp 空间不足也会同样的 exit 134**——大基因组每样本 SAM 10-20GB，并行 N 个会瞬间耗尽 /tmp（如 CNS1 仅 37G 可用）。必须先 `df -h /tmp` 判断：NFS 有空间时把 SAM 写 NFS（tmp_sam 目录）最稳；NFS 满才写 /tmp 再 cp 回。绝不能盲写 /tmp。重跑用 pending_samples.txt 只处理缺失 BAM 的样本
- **`-x` 必须传 index 前缀**（不是 fasta 路径）：传错时报 `"...genome.fasta" does not exist`（文件其实存在，hisat2 在找 genome.fasta.1.ht2）。批量前先 `ls ref/genome.{1..8}.ht2` 验证 8 个文件齐全
- **低 IO 批量配置**（用户明确要求时）：2 样本并行 × hisat2 -p 4，管道直排 sort（无 SAM），`-T /tmp/xxx` 排序临时分片写本地盘，NFS 只落最终 BAM；`set -e -o pipefail` 防失败样本产出空 BAM 被静默跳过；skip-if-exists 可断点续跑。约 25-30h/590Gb
- **杀批量任务要清扫孤儿进程**：kill xargs 后 hisat2（perl wrapper）会被 init 收养继续跑；用 `ps -ef | grep` 列 PID 后逐个 `kill -9`，反复直到数量为 0；不要用 `kill $(pgrep -f ...)` 形式（会被拒）。失败运行的半成品 bam/flagstat 要物理删除
- **CSI 索引**：染色体 >512Mb 时 BAI 索引会失败，用 `samtools index -c`
- **组织级过滤**：按组织内平均 TPM 过滤，不是全局平均
- **ashr 不可用时**：`lfcShrink` 可能失败，回退到未收缩的 LFC（p值不变）