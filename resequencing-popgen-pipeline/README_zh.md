# 重测序群体遗传全流程

从原始 FASTQ → VCF → 群体结构五图 + Fst/π/TajD/LD decay → HTML 报告。

> **SNP calling 硬性约定**：calling 一律用 **freebayes**，不用 `bcftools mpileup/call`。
> bcftools 只允许出现在 calling 之后的 VCF 工具链（concat / norm / index / stats / +fixploidy）。

## 环境要求

- **fb_calling** conda 环境：freebayes 1.3.10, bcftools 1.23.1, vcftools 0.1.17
- **pop_genetics** conda 环境：plink2, admixture, fasttree, angsd
- **R** 包：ggplot2, ggtree, ggtreeExtra, ggnewscale, ape, patchwork
- **minibwa**：快速 BWA-MEM 后继（`~/bin/minibwa`），短读模式 `-x sr`
- **PopLDdecay**：`/media/nfs1/hermes/miniforge3/bin/PopLDdecay`（不在 pop_genetics 环境里，用全路径）

## 目录结构

```
项目/req/
├── 01.rawdata/          # 原始 FASTQ 文件
├── 02.cleandata/        # fastp 质控后的 reads
├── 03.Mapping/          # minibwa 比对 + samtools 排序
│   ├── ref/             # 参考基因组（genome.fasta 软链接）
│   ├── bam/             # *.sorted.bam + .bai
│   └── summary/         # 比对统计
├── 04.Variants/         # freebayes 联合 calling + VCF 过滤
│   ├── regions/chunks/  # 100 kb 分块区间
│   ├── lists/           # bamlist.txt
│   ├── vcf/parts_raw/   # 块级 VCF
│   └── scripts/         # make_regions.py / freebayes_call_unit.sh ...
├── 05.population/       # 主分析（全部样本）
│   ├── ld_prune/        # LD pruning 位点集
│   ├── pca/  admixture/  tree/
├── 06.pv25_only/        # 亚群分析（可选）
├── 07.selection/        # 选择分析（Fst, π, TajD, LD decay）
└── 00.report/           # HTML 报告 + 图表
```

## 分析流程

### 第1步：质控（fastp）

统一采用 **duzhong 标准配置**（`--cut_front --cut_tail` 双向质量剪切，不用 `cut_right`）：

```bash
fastp \
    --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --n_base_limit 5 \
    --average_qual 20 \
    --length_required 36 \
    --thread 8 \
    -i R1.fq.gz -I R2.fq.gz \
    -o clean_R1.fq.gz -O clean_R2.fq.gz \
    --html sample_fastp.html --json sample_fastp.json
```

> DNBSEQ/MGI 的接头不是 Illumina（reads header `@ML…`），先 grep 确认再显式指定。
> 非 DNBSEQ 平台可去掉 `--adapter_sequence*` 两行，靠 `--detect_adapter_for_pe` 自动探测。

**参数严格度实测对比**（duzhong 面板 SRR37245294 前 40 万 read 对，Q30 94.62%）：

| 配置 | bases 保留 | Q30 | 被质量剔除 | 3'端 Q<20 碱基/20万read |
|------|-----------|-----|-----------|----------------------|
| fastp 默认 | 99.62% | 94.62 | 0 | 2803 |
| 旧参数（`cut_right` + `-u 40`） | 93.08% | 96.91 | 366 | 1254 |
| **duzhong 标准（本流程现用）** | 99.54% | 94.66 | 0 | 2436 |
| duzhong + `-u 20` | 95.74% | 95.99 | 30,674 | — |

坑：`-q`/`-u` 只做整条 read 判定、不裁碱基；`--cut_*` 才裁碱基；fastp 默认 `-q 15 -u 40`
几乎不过滤；数据质量差（Q30 < 90%）的项目要显式加 `--unqualified_percent_limit 20`。

### 第2步：比对（minibwa + samtools）
```bash
minibwa map -t 6 -x sr -R "@RG\tID:$s\tSM:$s\tPL:DNBSEQ\tLB:$s"     ref/genome.fasta ${s}_clean_R1.fastq.gz ${s}_clean_R2.fastq.gz   | samtools sort -@ 2 -m 2G -T /tmp/mbtmp_$s -o bam/$s.sorted.bam -
samtools index bam/$s.sorted.bam
```
> ⚠️ **管道符必须是 `|`**。写成 `||` 会静默把 SAM 写进日志（SGE 日志涨到 150 GB 而非 12.6 GB BAM）。
> 染色体 > 500 Mb（BAI 上限 2^29 = 536 Mb）必须用 `samtools index -c`（CSI）。
> `minibwa index` 不接受软链接路径，必须用真实路径。

### 第3步：样本物种核验（混样/错样检测，必做）
比对后按 `mapped_pct` 排序，< 90% 的样本一律核验（不要当成"比对参数问题"）：
```bash
zcat S_clean_R1.fastq.gz | head -400000 | gzip -1 > S.sub.fq.gz
sourmash sketch dna -p k=31,scaled=5000 -o panel.sig *.sub.fq.gz
sourmash compare -o cmp panel.sig
```
再用叶绿体 reads（`samtools view -q 30 BAM Pt`）比 NM 标签 identity，最后跨物种 BLAST / 回贴定案。

### 第4步：样本重复 / 混样核验（必做）
```bash
plink2 --pfile ld_prune/<panel> --make-king-table --out kinship/<panel>_king --threads 8
plink2 --pfile <panel> --het --out kinship/<panel>_het
```
判读标尺：KINSHIP ≥ 0.45 + 位点一致率 99.6%+ 才是真重复；混样看 F 值（实测混样 F = −0.83、杂合率 62%）。

### 第5步：变异检测（freebayes 联合 calling）
```bash
# 1) 100 kb 分块
python3 make_regions.py 100000 && python3 make_work_units.py 40

# 2) 每块的命令（核心参数）
freebayes -f ref/genome.fasta -r chr01:1-100000     -q 20 -m 20 --use-mapping-quality --genotype-qualities --ploidy 2     --min-alternate-count 2 --min-alternate-fraction 0.05     --bam-list lists/bamlist.txt -v vcf/parts_raw/chr01-1-100000.vcf

# 3) SGE 数组提交（40 单元 × -pe smp 6 = 240 槽 → 恰好 40 个 freebayes 进程）
qsub -q all.q@CNS2 -pe smp 6 -l h_vmem=10G -cwd -o logs/fb_array -e logs/fb_array      -t 1-40 scripts/freebayes_call_unit_array.sh
```
> 耗时 ∝ **候选位点数 × 样品数**，不是数据量。346 样本 × 1.015 Gb 实测中位 **32 min/块**，
> 全量约 135 h（5.6 天）。不要用"复制 N 份 BAM"外推（会低估 19 倍）。
> 模板脚本见 `templates/`：`make_regions.py`、`make_work_units.py`、`freebayes_call_unit.sh`、
> `freebayes_call_unit_array.sub`。

### 第6步：VCF 过滤
| 参数 | 阈值 | 说明 |
|------|------|------|
| QUAL | ≥20 | 碱基质量过滤 |
| GQ | ≥20 | 基因型质量 |
| DP | [4, 1000] | 去除低/高覆盖度伪影 |
| F_MISSING | ≤0.2 | 最大 20% 缺失基因型 |
| MAF | ≥0.05 | 最小等位基因频率 |
| Alleles | 仅双等位 SNP | 去除 indel/多等位位点 |

```bash
bcftools concat -a -f parts_raw/*.vcf.gz | bcftools norm -d exact -Oz -o merged.vcf.gz
bcftools index -c merged.vcf.gz
bcftools view -i 'QUAL>=20 && FORMAT/GQ>=20 && FORMAT/DP>=4 && FORMAT/DP<=1000' merged.vcf.gz   | bcftools view -i 'F_MISSING<=0.2' | bcftools view -m2 -M2 -v snps   | bcftools view -i 'MAF>=0.05' > filtered.vcf.gz
```
> ⚠️ **Tajima's D 必须在未做 MAF 过滤的位点集上算**（MAF≥5% 会把 D 系统抬高）。
> ⚠️ **vcftools 只出每 SNP 的 `.weir.fst`**（0.1.17 无 `--fst-window`），窗口值自行聚合。
> ⚠️ **混合 `.` 与 `./.` 缺失编码会让 vcftools 静默丢位点** → 用 `bcftools +fixploidy -- -f 2` 归一化。

### 第7步：群体结构分析
**PCA（plink2）：**
```bash
plink2 --vcf filtered.vcf.gz --set-all-var-ids @:# --rm-dup   --pca 10 --out pca_results
```
**ADMIXTURE：**
```bash
plink2 --vcf filtered.vcf.gz --set-all-var-ids @:# --rm-dup --make-bed --out pop
for K in 2 3 4 5 6; do admixture --cv pop.bed $K -j8 | tee log${K}.out; done
```
**树（IBS + ape::nj / hclust UPGMA 出图，FastTree 仅存档）：**
```bash
plink2 --vcf filtered.vcf.gz --set-all-var-ids @:# --rm-dup --export A --out ibs
```
> plink2 **没有** `--distance`；`.eigenvec` 带表头要 `header=TRUE`；`.bed` 文件头 3 字节要 `raw[3:]`。

### 第8步：选择分析
| 指标 | 命令 | 窗口 |
|------|------|------|
| π | `vcftools --gzvcf VCF --keep G.txt --window-pi 50000 --window-pi-step 10000` | 50 kb/10 kb |
| Tajima's D | `vcftools --gzvcf VCF --keep G.txt --TajimaD 50000` | 50 kb |
| Fst | `vcftools --gzvcf VCF --weir-fst-pop A.txt --weir-fst-pop B.txt --fst-window-size 50000 --fst-window-step 10000` | 50 kb/10 kb |
| ROH | `plink --bfile bed --keep G_fid.txt --homozyg --homozyg-window-snp 50 --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-window-missing 5` | 需 PLINK 1.9 |
| LD 衰减 | PopLDdecay → log 分箱 + 池化 r² + `nls` 拟合 | 到 500 kb |

> **选择信号打分：缺数据的窗口必须剔除**（窗口内变异数 ≥10、两群体 π 均 >0 且非 nan、负 Fst 截断为 0、打分用秩变换），否则 top 榜全是空窗口。

## 可视化（五图模板）

**图1：** 环形系统发育树（整枝着色，邻近向上传播）
**图2：** 环形树 + K=2..N Admixture 外圈环（自算极坐标布局，不依赖 ggtreeExtra）
**图3：** PCA（PC1/2 + PC1/3，群体着色）
**图4：** Admixture K=2-4（组内按祖先比例排序）
**图CV：** 交叉验证误差曲线

绘图统一：**英文 + Arial + 600 dpi PNG + cairo_pdf 矢量**；报告里每张图配"怎么算好"的解读段。
模板：`templates/circular_tree_admixture_rings.R`、`templates/tree_structure_combined.R`。

## 关键避坑（完整 26 条见 SKILL.md）

1. **plink2 LD**：必须 `--set-all-var-ids @:# --rm-dup`
2. **ADMIXTURE**：染色体 ID 必须是整数；.Q 文件必须在当前工作目录
3. **ANGSD `-fold`**：仅在 `realSFS` 中使用
4. **ANGSD 2D SFS**：>140 样本不可行，改用 vcftools `--weir-fst-pop`（分钟级 vs 卡死）
5. **PCA 特征值**：除以所有特征值之和才是方差贡献率
6. **Q 矩阵对齐**：`Q[match(tree$tip.label, fam_order), ]`
7. **枝干着色**：邻近向上传播（子→父），非多数投票
8. **`xargs -P` 里函数用到的每个变量都必须 `export`**（漏一个 → 路径塌成 `/logs/...`，197 样本全 FAIL）
9. **长任务**：用 `nohup cmd > log 2>&1 &` 避免包装器超时
10. **PopLDdecay**：全路径调用；`comment.char=""` 读 `.stat.gz`；避免 `-SubPop`（segfault）
11. **`pkill -f` 会杀掉自己的 shell** → 用 `kill -9 <明确PID>`
12. **曼哈顿图 Bonferroni 分母是标记数**（不是基因组长度）；点 size 1.2，别抽样
13. **`minibwa | samtools sort` 的 `|` 写成 `||` 会静默烧盘**（SAM 落进日志）

## 模板脚本（templates/）

| 文件 | 作用 |
|------|------|
| `make_regions.py` | .fai → 100 kb 分块区间 `regions/chunks/<chr>.regions` |
| `make_work_units.py` | 分块轮转均分成 N 个工作单元（N = 想要并发数，惯例 40） |
| `freebayes_call_unit.sh` | 跑一个工作单元（块级独立 VCF、已存在则 SKIP、可续跑） |
| `freebayes_call_unit_array.sub` | SGE 数组 N 任务 × `-pe smp 6` |
| `circular_tree_admixture_rings.R` | 环形树 + K=2..N admixture 外圈（自算极坐标） |
| `tree_structure_combined.R` | 矩形 NJ 树 + K=2..N admixture 条形图合一 |

## 引用

软件与方法原文、BibTeX、Methods 模板句（中英文）见 **SKILL.md 的 Citations 章节**。

### 方法模板句（中文，摘要版）

> 原始 reads 经 fastp（Chen et al., 2018）质控后，用 minibwa（Li & Homer, 2026）比对至参考基因组。
> 变异检测采用 freebayes（Garrison & Marth, 2012）联合 calling，并使用 bcftools（Danecek et al., 2021）
> 进行硬过滤（QUAL ≥ 20，基因型深度 DP 4–1000，最小等位基因频率 MAF ≥ 0.05）。群体结构通过 PCA
> （plink2；Chang et al., 2015）、模型聚类（ADMIXTURE；Alexander et al., 2009）和基于 IBS 距离的邻接树
> （ape；Paradis & Schliep, 2018）评估。群体分化（Fst）、核苷酸多样性（π）和 Tajima's D 使用 vcftools
> （Danecek et al., 2011）计算，连锁不平衡衰减使用 PopLDdecay（Zhang et al., 2019）估算。
> 可视化使用 ggtree（Yu et al., 2017）和 ggplot2（Wickham, 2016）。

## 许可

内部流程沉淀，公开使用请注明来源。
