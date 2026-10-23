# RNA-seq：HISAT2 + featureCounts + DESeq2 流程

植物批量 RNA-seq 标准流程：raw FASTQ → 差异表达基因 → GO/KEGG 富集，含全部实测坑点。

## 流程总览

```
1.rawdata/        → fastp 质控
2.clean/          → HISAT2 比对
3.Mapping/        → featureCounts 定量
4.Quantification/ → TPM 计算 + 组织级低表达过滤
5.DEG/            → DESeq2 差异表达
6.Enrichment/     → clusterProfiler GO/KEGG 富集
```

## 环境

- fastp（HISAT2 需要 PATH 里有 `python`，CNS 上要 `export PATH=/media/nfs1/hermes/miniforge3/bin:$PATH`）
- HISAT2 + samtools（`samtools index -c` 用 CSI，支持 >512 Mb 染色体）
- featureCounts（subread）
- R：DESeq2、clusterProfiler、ggplot2、pheatmap、GO.db、AnnotationDbi
- eggNOG-mapper（功能注释，DB 在 `/media/nfs1/hermes/db/eggnog`）

## 关键参数

| 步骤 | 参数 | 理由 |
|------|------|------|
| fastp | `--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 -q 20 -n 5 --average_qual 20 --length_required 60 --compression 6` | duzhong 标准（2026-10-23 长度门槛 36→60 收紧） |
| hisat2-build | `--ss splice_sites.txt --exon exons.txt` | 必须先 GFF3→GTF |
| hisat2 | `-p 16`，**不加** `--dta` | featureCounts/DESeq2 不需要，省 ~10 GB 内存 |
| hisat2 | `--rna-strandness` 省略（unstranded） | 链特异文库才用 FR/RF |
| samtools sort | `-@ 4 -m 4G` | 4 线程 × 4 G = 16 G 上限 |
| samtools index | `-c`（CSI） | BAI 上限 2^29 = 536 Mb |
| featureCounts | `-g gene_id -t exon -s 0 -p --countReadPairs` | 注释必须是 GTF |
| TPM 过滤 | 任一组织平均 TPM ≥ 1 | **不是全局平均** |

## 目录结构

```
project/
├── 1.rawdata/            *.fastq.gz
├── 2.clean/              *_clean_R1.fastq.gz / *_clean_R2.fastq.gz + logs/ + run_fastp.sh
├── 3.Mapping/
│   ├── ref/              genome.fasta + genome.{1..8}.ht2
│   ├── bam/              *.sorted.bam + *.sorted.bam.csi
│   ├── summary/          *.flagstat + clean_and_mapping_summary.tsv
│   ├── logs/             *.hisat2.log
│   └── scripts/          run_hisat2.sh
├── 4.Quantification/featureCounts/
│   ├── counts_raw.txt
│   ├── counts_matrix.tsv / tpm_matrix.tsv
│   └── counts_matrix_filtered.tsv / tpm_matrix_filtered.tsv
├── 5.DEG/
│   ├── run_deseq2.R
│   ├── DEG_summary.csv / DEG_*_vs_*.csv
│   └── volcano_*.png / MA_*.png / PCA_all_samples.png / heatmap_top_DEGs.png
└── 6.Enrichment/         clusterProfiler 结果 + 图
```

## 关键坑点（完整清单见 SKILL.md）

### fastp / QC
- **参数用 duzhong 标准配置**：`--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --n_base_limit 5 --average_qual 20 --length_required 60 --compression 6`
- **`-q`/`-u` 只做整条 read 判定、不裁碱基**；裁碱基必须给 `--cut_front`/`--cut_tail`。fastp 默认 `-q 15 -u 40` 几乎不过滤（实测整批一条不丢）
- **fastp 默认就已去接头**（`--disable_adapter_trimming` 才是关闭），显式写接头的价值是跨平台可移植，不是清洁度
- **数据质量差（Q30 < 90%）的项目要加 `--unqualified_percent_limit 20`**，否则低质量 read 全部放行
- **数据量报 bases (Gb)，不报压缩文件 GB**——从 fastp JSON 的 `after_filtering.total_bases` 读，不要重新解压去数
- **fastp 的 `total_reads` 是单端 read 数**（PE 时 = pair 数 × 2）；PE pair 数要从 HISAT2 log 取
- **假 PE 检测**：ENA 标 PAIRED 但实际单端 → 看 FASTQ header 的 read-index 是 1,2,1,2（真 PE）还是 1,2,3,4（单端）
- ENA/NCBI API 在 CNS 上常被墙 → 用 `vdb-dump --info` 从 .sra 文件拿平台信息

### HISAT2 建索引
- **GFF3 必须先转 GTF**（`gffread annotation.gff3 -T -o annotation.gtf`），否则 splice_sites/exons 是 0 字节
- **`python` 必须在 PATH**——否则 hisat2-build "瞬间完成"但生成空 .ht2
- **NFS >95% 满时索引静默损坏** → 先在 /tmp 建好再 cp
- 单倍型分辨基因组索引构建慢（SA 阶段几分钟到 2 小时），看 .sa 文件数判断进度，别看日志

### HISAT2 比对
- **⚠️ exit 134 有两个完全不同的原因**：NFS 满 **或 /tmp 满**。大基因组每样本 SAM 10-21 GB，并行 N 个瞬间耗尽 /tmp（CNS1 只有 37 G）。**先 `df -h /tmp` 再决定 SAM 写哪**——NFS 有空间就写 NFS（`tmp_sam/`），别盲写 /tmp
- **`-x` 传 index 前缀不是 fasta 路径**：报错说 "genome.fasta does not exist"，其实他在找 `genome.fasta.1.ht2`
- **样本名带括号 `H4(3)` → 静默产出 0-read BAM 且 exit 0**。必须 sanitize 成 `[A-Za-z0-9._-]`，并加 post-mapping guard 检查 read 数与 BAM 大小
- **`--no-unal` 让 flagstat 报 100% mapped**，真实比对率看 `.hisat2.log` 的 overall alignment rate
- **NFS 管道 SIGPIPE（exit 141）**：先写 SAM 再 sort
- **CSI 索引**：染色体 >512 Mb 必须 `samtools index -c`
- **杀批量任务要清扫孤儿进程**：kill xargs 后 hisat2 perl wrapper 被 init 收养继续跑。用 `ps -ef | grep` 列 PID 后逐个 `kill -9`，别用 `kill $(pgrep -f ...)`（会被拒）
- **低 IO 配置**（用户要求时）：2 样本并行 × `-p 4`，管道直排 sort（不落 SAM），`-T /tmp/xxx` 排序临时片写本地盘；`set -e -o pipefail` 防失败样本产出空 BAM 被静默跳过

### featureCounts
- **必须 GTF**，用 GFF3 会报 `failed to find the gene identifier attribute`

### TPM / 过滤
- **组织级过滤**：任一组织平均 TPM ≥ 1 就保留。全局平均会删掉组织特异基因
- 矩阵行名带转录本/版本后缀（`Seita.6G138300.1`）→ 匹配前用 `sub("[.][0-9]+$", "", rownames(tpm))` 剥掉

### DESeq2
- `lfcShrink(type="ashr")` 可能失败 → `tryCatch` 回退到未收缩 LFC（p 值不变）
- `relevel(coldata$condition, ref="root1")` 设基准

### 富集
- **`enricher()` 的 TERM2GENE 列名必须是大写 `TERM` 和 `GENE`**，写成 GO/Gene 会静默返回 0 条结果
- **eggNOG 输出必须显式给 `col.names`**（`header=FALSE` + `comment.char="#"`），否则列名错位
- **R 正则写脚本用 `[.]` 而不是 `\\.`**，避免转义地狱
- **KEGG 用 `KEGG_Pathway`（ko#####）不是 `KEGG_ko`（K#####）**；排除动物通路前缀 `ko05/ko047/ko049/ko046/ko043`
- **背景集 = 该比较里 pvalue 非 NA 的全部基因**，不是全基因组
- **GO 只画 BP 子本体**，并传 TERM2NAME 显示可读名称

## 模板脚本

| 文件 | 作用 |
|------|------|
| `templates/run_fastp.sh` | PE 批量质控（幂等、并发可调） |
| `templates/build_hisat2_index.sh` | GFF3→GTF → splice_sites/exons → 建索引（含 python PATH 修复） |
| `templates/run_hisat2.sh` | 批量比对（两种模式：SAM-on-NFS / 低IO管道直排） |
| `templates/calc_tpm_filter.py` | TPM 计算 + 组织级过滤，输出 4 个矩阵 |
| `templates/run_deseq2.R` | DESeq2 多比较 + volcano/MA/PCA/heatmap |
| `templates/run_multigroup_enrichment.R` | 多组 clusterProfiler 富集（x=group 散点图） |

## 引用

软件与方法原文、BibTeX、中英文 Methods 模板句见 **SKILL.md 的 Citations 章节**。

### 方法模板句（中文）

> 原始 reads 经 fastp（Chen et al., 2018）质控后，使用 HISAT2（Kim et al., 2019）比对至参考基因组。
> 基因水平定量采用 featureCounts v2.1.1（Liao et al., 2014），参数为 `-g gene_id -t exon -s 0 -p --countReadPairs`。
> TPM 由原始 counts 按基因长度与文库大小归一化得到；按组织分组过滤低表达基因（任一组织平均 TPM ≥ 1 则保留）。
> 差异表达用 DESeq2（Love et al., 2014）分析，设计为 `~ condition`，以 root1 为基准，padj < 0.05 视为显著差异表达。
> 功能富集用 clusterProfiler（Wu et al., 2021），注释来源为 eggNOG-mapper（Cantalapiedra et al., 2021）导出的 GO 与 KEGG。

## 许可

内部流程沉淀，公开使用请注明来源。
