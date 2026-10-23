---
name: resequencing-popgen-pipeline
description: "End-to-end plant resequencing population genetics pipeline: fastp QC, minibwa mapping, freebayes joint calling, VCF norm/filter, plink2 LD+PCA, ADMIXTURE, FastTree2, ANGSD pi/TajD/Fst, ggtree visualization (5-figure template), HTML report. Okabe-Ito colorblind palette."
category: bioinformatics
triggers:
  - "重测序分析"
  - "resequencing popgen"
  - "群体遗传全流程"
  - "群体结构三剑客"
  - "freebayes joint calling pipeline"
---

# 重测序群体遗传全流程

从 raw FASTQ → VCF → 群体结构五图 + Fst/π/TajD/LD decay → HTML 报告。

## 环境

- fb_calling: freebayes, bcftools, vcftools
- pop_genetics: plink2, admixture, fasttree, angsd
- R: ggplot2, ggtree, ggtreeExtra, ggnewscale, ape, patchwork

## 目录结构

```
project/req/
├── 01.rawdata/ → 02.cleandata/ → 03.Mapping/ → 04.Variants/
├── 05.population/ (主分析: SRR+Pv25)
├── 06.pv25_only/   (亚群分析: Pv25 only)
├── 07.selection/   (选择分析: Fst A/B, LD decay, pi/TajD)
└── 00.report/
```

## Stage 2-3 实测（无花果 req 面板：199 WGS, DNBSEQ PE150, 2026-09）

### 接头与 fastp 参数

DNBSEQ/MGI 数据接头不是 Illumina（reads header `@ML…`）。先 grep 确认，再显式指定。

**统一采用 duzhong 面板实测的标准配置**（`--cut_front --cut_tail` 双向质量剪切，
不用 `cut_right`）：

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
    --thread 8
```

> `--detect_adapter_for_pe` 会额外探测 R2 的接头，比只写 `--adapter_sequence` 多探一次，
> 换来的是非 DNBSEQ 平台也能自动找到接头（对同批数据两种写法结果几乎一致，见下方实测）。

199 样本实测（无花果 req 面板，2026-09，DNBSEQ PE150）：raw 2545 Gb → clean 2424 Gb
（保留 95.3%），每样本保留率均值 97.78%，Q30 97.08%，GC 35.85%；平均 raw 12.8 Gb/样本；
每样本 wall 200-252 s（8 线程，4 进程并行 = 32 核，全部 199 样本 3 h 26 min）。

### fastp 参数：默认 vs 各配置实测对比（2026-10 实测，duzhong 面板 SRR37245294 前 40 万 read 对）

同一份真实数据跑 4 套参数。**关键结论：fastp 默认参数已经把接头去掉了**
（`--disable_adapter_trimming` 才是关闭接头去除），所以"显式写接头"相对默认的增量很小；
真正决定严格度的是 `--cut_*` 质量剪切和 `-u`（unqualified_percent_limit）。

| 配置 | bases 保留 | 平均读长 | Q30 | 被质量剔除 | 3'端 Q<20 碱基数/20万read |
|------|-----------|---------|-----|-----------|----------------------|
| **fastp 默认**（无参数） | 99.62% | 149.0 | 94.62 | 0 | 2803 |
| `--detect_adapter_for_pe` only | 99.61% | 149.0 | 94.62 | 0 | 2803 |
| **旧重测序参数**（`cut_right` + `-u 40`） | 93.08% | 148.0 | **96.91** | 366 | **1254** |
| **duzhong 标准**（`cut_front+cut_tail` + `-q20`） | 99.54% | 149.0 | 94.66 | 0 | 2436 |
| duzhong + `-u 20` | 95.74% | 149.0 | 95.99 | 30,674 | — |

读法与坑：

- **`-q`/`-u` 只做"整条 read 判定"，不裁单个碱基**。`-q 20 -u 40` 意为"允许 40% 碱基
  不合格"，所以质量差的 read 基本都活着通过。真正裁碱基必须给 `--cut_front`/`--cut_tail`。
- **`cut_right` 极其激进**：遇到一个低质量窗口就砍掉它**及其右侧全部**。它在这份数据上
  保留率最低（93.08%）但 Q30 提升最多（96.91%），且低质量碱基残留最少（1254）。
  已统一改为 `cut_front + cut_tail`（温和裁剪，两端各自裁到合格为止）。
- **`-u 20` 才是让"低质量 read 被判死"的那一项**：duzhong 原配置没写 `-u`（即默认 u40），
  所以它在这批高质量数据（Q30 94.6%）上表现得几乎等同默认（只丢 526 条 <36bp 的短 read，
  质量剪切实际裁掉 0 个碱基）。**数据质量差的项目（Q30 < 90%）应显式加 `-u 20`。**
- `--length_required`：duzhong 用 60（DNBSEQ PE150），本流程用 36。**取样长度的 1/3~1/2
  即可**，不要求全长——太严会因文库长度波动误删样本。
- fastp 的 `-q` 默认是 15、`-u` 默认是 40、`-l` 默认是 15，**都不适合直接用**：
  `-q 15` 等于几乎不过滤（实测 lowq=0，整批一条不丢）。
- **接头去重效果**：`--detect_adapter_for_pe` 比纯默认多探到 1,234 个 read 的接头
  （17,424 vs 16,190），多裁 14,063 个碱基——量级很小但方向正确。

### 参考：a 套 + 叶绿体 + 线粒体

单倍型分辨组装（`chrNNa`/`chrNNb` + `Pt` + `Mt1`/`Mt2` 在同一 FASTA）取 a 套 + 细胞器：

```bash
printf 'chr%02da\n' {1..13} > ids.txt; echo -e 'Pt\nMt1\nMt2' >> ids.txt
samtools faidx genome_all.fa -r ids.txt -o ref/genome.fasta   # 311.0 Mb chrA + 641 kb 细胞器
samtools faidx ref/genome.fasta
minibwa index -t 8 ref/genome.fasta ref/genome.fasta          # 41 s, 646 MB .mbw
```

### minibwa 比对（DNA 重测序）

```bash
minibwa map -t 6 -x sr -R "@RG\tID:$s\tSM:$s\tPL:DNBSEQ\tLB:$s" \
    ref/genome.fasta ${s}_clean_R1.fastq.gz ${s}_clean_R2.fastq.gz \
  | samtools sort -@ 2 -m 2G -T /tmp/mbtmp_$s -o bam/$s.sorted.bam -
samtools index bam/$s.sorted.bam          # 染色体 ≤36 Mb，BAI 足够（无需 CSI）
```

实测（6 线程/样本，12.8 Gb clean 输入）：34 min，~60k reads/s，峰值 RSS 2.4 GB，BAM 4.1 GB，
mapped 99.52%，properly paired 95.43%，singleton 0.19%，supplementary 1.3%（单倍型 A 参考上
B 套 reads 的表现，正常）。chrA 深度 ≈ 48×。并行度 8 × 6 核 = 48 核 → 199 样本约 14 h。
`-x sr` 用于 150 bp 短读；`samtools sort` 临时文件放 CNS2 本地盘（`-T /tmp/...`）避免 NFS 抖动。

细胞器 reads（`samtools idxstats`，1 条命令即可，不需额外扫描）：叶片 DNA 中叶绿体读数极高
（Pt 1.24e7 reads ≈ 11380×；Mt1+Mt2 2.3e6 ≈ 710×），分析前先想清楚是否要去除细胞器 reads。

## Stage 3.5: 样本物种核验（混样/错样检测，必做）

mapping 完成后按 `mapped_pct` 和 `est_depth` 排序，低于 90% 的样本一律做物种核验，
不要当成“比对参数问题”。无花果面板 2026-09 实测：WHG141-144 只有 10-12% mapped
（其余 195 个样本 97.6%±，深度 37×），运行方向标签以外的 QC 指标（Q30 97.5%、GC 36%、
insert 266 bp、dup 0.6%）全部正常 → 不是建库/测序问题，是“材料不是本物种”。

三步定位法（每步几分钟，从便宜到贵）：

1. **k-mer 相似度（sourmash，最快）**——判断 4 个可疑样本是否同一物种、与正常样本差多远：

```bash
zcat S_clean_R1.fastq.gz | head -400000 | gzip -1 > S.sub.fq.gz   # 10 万条足够
sourmash sketch dna -p k=31,scaled=5000 -o panel.sig *.sub.fq.gz
sourmash compare -o cmp panel.sig
```

实测：可疑 4 个互相 0.034-0.037，与正常无花果样本 0.001-0.003（差 15 倍）；正常样本之间 0.049。
→ 4 个可疑样本是一个物种，且不是无花果。

2. **细胞器一致性**——`samtools view -q 30 BAM Pt` 取叶绿体 reads，比较 NM 标签的平均
identity：正常样本 99.84%，可疑样本 93.2%（叶绿体保守，能跨物种比上；核基因组比不上）。

3. **跨物种 BLAST / 反证比对（定案）**——取可疑样本的叶绿体 reads（NCBИ URLAPI，
不要用 `Bio.Blast.NCBIWWW.qblast`，它拒收 `task=` 参数）比对 nt；再把 10 万 read pairs
回贴到候选物种参考基因组：

```bash
minibwa map -t 6 -x sr ref_candidate.fa S.p_R1.fq.gz S.p_R2.fq.gz | samtools flagstat -
```

实测：4 个可疑样本贴到 *Vitis vinifera* PN40024 T2T 为 95.2-98.1%，无花果样本只有 15.9-17.1%
（只靠叶绿体保守区）→ 结论：这 4 个文库是葡萄材料，混入了无花果批次。

参考库：`/media/nfs1/hermes/db/Vitis_vinifera/ncbi_dataset/data/GCF_030704535.1/GCF_030704535.1_ASM3070453v1_genomic.fna`
（Vitis vinifera PN40024 T2T，GCF_030704535.1，21 条序列）。需要别的物种同理放 `/media/nfs1/hermes/db/<Species>/`。

## Stage 3.6: 样本重复 / 混样核验（重测序面板必做，2026-09 无花果实测）

面板里“同一个体被重复取样 / 两条 lane 被合并 / 标签错配”会直接扭曲 PCA、ADMIXTURE、树、LD、Fst。
在跑群体结构**之前**先做一遍，四步、由快到慢：

1. **KING 亲缘扫描**（分钟级，全面板一次搞定）

```bash
plink2 --pfile ld_prune/<panel> --make-king-table --out kinship/<panel>_king --threads 8
# .kin0 列: #IID1 IID2 NSNP HETHET IBS0 KINSHIP；同一个体 KINSHIP≈0.5 且 IBS0≈0
```

用 KINSHIP≥0.45 建图求连通分量 → 得到“重复组”。无花果 194 样本实测：**37 组、147 个样本**，
独立基因型只有 84 个（最小组 16 个样本同一基因型）。注意：如果某样本本身是混合样品，它的 KINSHIP
对所有人都被抬高（≈0.26），不要当成亲缘关系——先用第 3 步的 F 值把混样剔掉。

2. **基因型一致率定量**（区分“真重复”和“近亲”）

把面板与待检样本在**同一批剪枝位点**上重定基因型再比一致率。实测标尺（73,790 剪枝位点）：

| 配对类型 | 位点一致率 | 位点差异率 |
|---|---|---|
| 同一个体（真重复） | 99.6–99.95% | 0.36–0.68% |
| 面板内随机两个个体 | ~55% | 45%（中位） |

**不要用 KINSHIP 阈值当唯一证据**：同一批野生材料里随机配对也能到 0.26，只有 ≥0.45 且伴随
99.6%+ 一致率才是真重复。

3. **混样（嵌合体）识别：看 F 值**（秒级）

```bash
plink2 --pfile <panel> --het --out kinship/<panel>_het      # 列: #IID O(HOM) E(HOM) OBS_CT F
```

把两条不同个体的 lane 合并成一条 BAM，得到的“样品”在杂合位点上全部变杂合：实测 F = **-0.83**、
杂合率 **62%**（panel 范围 -0.17~+0.16、杂合率 21–62%，中位 36%）。F 明显偏离面板范围就是混样/合并错误。

4. **读段级核验**（区分“同一批数据重复交付”vs“独立文库”）

取一个 1 Mb 窗口，把两样本 read1 的 `(POS, SEQ, 链)` 做成集合求交：

```bash
samtools view A.sorted.bam chr01a:1000000-2000000   # 取 read1 (跳过 flag 0x100/0x4)，组 (POS,SEQ,strand)
```

实测：**同一文库拆两条 lane / 同一批数据重复交付 → 重叠率极高**；**独立文库（即使是同一个体）→
只有背景重叠 5–9%**（不同样本之间也是 9%）。所以“重叠低”只能排除“同批数据重复交付”，
要证明“同一个体”还得看第 2 步的一致率。

5. **PCA 交叉验证**（图形化一眼看出）

重复组成员在 PC1–PC5 上几乎重合：实测组内最大欧氏距离 **0.0002–0.0004**，随机 5 样本 **0.158**（差 500–800 倍）。
图上一坨完全重叠的点就是重复样本。

**处理建议**（先问用户，不要自行决定）：每组保留 1 个代表得到独立基因型集再跑群体结构；
被合并错的样本要拆回各 lane 重新 calling。

### 坑

- `bcftools mpileup -R regions.tsv`：**三列文件（chr/start/end）按 1-based 闭区间解释**，
  `pos-1  pos` 会变成 2 bp → 每个位点出两条记录、且第一条是邻位碱基（取错碱基，基因型全乱）。
  要单碱基位点列表就用**两列（CHROM POS，1-based）**，并用 `bcftools view -H | wc -l` 核对记录数
  是否等于请求位点数。
- 用面板 VCF 的位点去比对新样本时，两边必须用**同一套 REF/ALT 编码**；最稳的做法是把新样本 BAM
  在面板位点上重新 mpileup+call，再与面板矩阵（`plink2 --export A`）逐位点比。

### 修复被错误合并的样品（不重跑全量 calling 的做法，2026-09 无花果 029 实测）

背景：两个文库被 `samtools merge` 当成一个样品送进 freebayes，得到的“样品”是嵌合体（F=-0.83、杂合率 62%）。
**不需要重跑 mapping**（两条 lane 的 BAM 本来就分开且完好），也**不需要重跑 50 h 全量 calling**：

> 关键前提：合并 BAM 的 reads = 两条 lane reads 的并集，所以 freebayes 的**位点发现与等位基因集合不变**，
> 只需重算“每条 lane 各自的基因型”。

流程（3227 块 × 40 并行，约 26 min；全流程含拼装+合并+过滤约 3 h）：

1. 原 VCF 去掉嵌合列：`bcftools view -s ^<嵌合样品名>`；
2. 用 `freebayes -f ref -r <chunk> -@ <该块位点>.vcf.gz -l --genotype-qualities \
   -q 20 -m 20 --use-mapping-quality --ploidy 2 --min-alternate-count 2 --min-alternate-fraction 0.05 \
   --bam-list <两条lane的bamlist> -v out.vcf` 重定基因型（`-@` + `-l` = 只在给定位点上、只用给定等位基因）；
3. 拼装 → 对两条 lane 做与原 step1 相同的样品级硬过滤（DP<4|>1000|GQ<20 置缺失）→ `bcftools merge` 与去嵌合后的 193 样品合并；
4. 重跑位点级过滤（`vcftools --max-missing 0.8 --maf 0.05`）→ 得到 N 样品 VCF。

实测数字：原子化位点 5,969,415 → 合并后仍 5,969,415（**位点数必须完全一致**，多出来就是表示不一致）；
最终过滤位点 1,886,125（原 194 样品版 1,886,134，差 9 个属正常）；两条 lane 杂合率 20.1%/37.9%（原嵌合体 79.6%）；
F 值 -0.019/+0.052（原 -0.832）。

## 坑（重定基因型专有）

- **MNP 表示不一致会静默把位点变并集**：原 step1 里约 13% 的记录是 MNP（如 `AG→AA`、`CCCTAAACCC→CCCTAACCCC`，
  本质是单个替换碱基被 pad 成多碱基）。freebayes 重定基因型会输出等价的**单碱基**形式（`A→C`），两边 REF/ALT 不同 →
  `bcftools merge` 各留一条，位点从 5,969,415 虚增到 6,444,313（+47 万）。**解法：先 `bcftools norm --atomize` 把原记录
  原子化成单碱基形式**（实测位点数不变，MNP↔SNP 是 1:1），再在同一批位点上重定基因型。
- `bcftools query -f '[\t%GT]\n'` 输出**以 tab 开头**，所以 awk 里基因型是 `$2`、`$3`，不是 `$1`（踩过两次，
  表现为“杂合率 0.0%”“一致率 0.00%”这类荒谬值）。
- **统计大 VCF 的位点数用 `bcftools index -n`**（读索引，瞬时）；`bcftools view -H | wc -l` 在 13 GB 文件上要 10+ 分钟，
  写在脚本里会把流程凭空拖长半小时。
- **跨群体 Fst 的接力脚本要写成"等上游跑完再动手"**（`while pgrep -f '<上游进程名>' >/dev/null; do sleep 120; done`），
  再用 `bcftools index -n` 验证 VCF 完整后串起来；这样不必人工盯着长任务。注意 pgrep 的模式不要匹配到脚本自己。
- **PopLDdecay 输出文件先出现 20 字节空 gz 是正常的**：它在算完所有 SNP 对之后才写入内容，
  4.17M 位点单线程要 1.5 h+（`-MaxDist 500`）。不要因为看到 20 字节就断定"又失败了"（真失败时进程已退出）。
- **PopLDdecay 的 `.stat.gz` 是逐 1 bp 一行（500 kb → 500,001 行），不要直接拿来算"LD 衰减到 r²=0.2 的距离"**：
  `-MaxDist 500` 的输出里 1–4 bp 那几行是假象（每行只有 ~1.1–1.5 万对位点，r²≈0.197），到 5 bp 直接跳到 0.714
  （4.9 万对）——按行号或 `dist == min(dist)` 取数会得到非单调的荒谬曲线。正确做法：
  ① 截掉 `< 100 bp` 的行；② 按几何分箱（100b/200b/500b/1k/2k/5k/10k/20k/50k/100k/200k/500k）聚合，
  每箱用 **池化均值 `sum(Sum_r^2)/sum(NumberPairs)`**（不是把每行的 Mean_r^2 再平均）；③ 在分箱值上做
  `nls(r2 ~ a/(1+b*d)^c)`；④ 报 r²=0.2 的距离并标注是否外推（`c` 很小时 r²=0.1 的外推值不可信）。
  实测（224 份豌豆，4,168,721 SNP）：1 kb 0.53 / 10 kb 0.42 / 100 kb 0.29 / 500 kb 0.20，
  拟合 r²=0.2 在 ~434 kb；r² 降到 0.1 远超 500 kb → 需要 `-MaxDist 2000` 重跑才能与文献（豌豆 ~1.8 Mb）比较。
- `bcftools merge` 默认 `-m both` 会把同一位点的记录合成一条（样品列并起来）；`-m none` 则不合并（会输出两条记录），
  合并样品列时必须用默认值。
- `bcftools view -s A,B` 的输出顺序 = 请求顺序（不是表头顺序），可放心按下标取列。
- 验证新样品的基因型是否真的算对：拿它与 panel 里最相似的样品比一致率（本例 029-1 ≡ WHG004 99.96%、
  029-2 ≡ WHG030 组 99.97%），再用 **BAM 层面的等位平衡**独立验证（在 A 的杂合位点上，B 的 reads 也应呈杂合；
  本例 WHG004 50.8% vs WHG029-2 14.1%）。

### 进化树 + 群体结构 合二为一（用户偏好样式，2026-09）

上：矩形 NJ 树（根在上、叶在下，叶标签旋转 90°）；下：K=2..N 的 admixture 条形图（`facet_grid(K ~ .)`）。
两块必须用**同一个 x 尺度**（x = 叶序号 1..n，`scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0,0))`），
对齐才成立。模板：`templates/tree_structure_combined.R`（自算矩形布局，不依赖 ggtreeExtra）。

要点：

- 矩形布局：横向段 `(x_parent → x_child, y = Y_parent)`，纵向段 `(x_child, Y_parent → Y_child)`；
  `Y = (到根距离 / 最大值)^POW`（POW=0.45，NJ 树内部分支极短，不压缩看不出结构）；叶标签放在 `y = 1.015` 处
  并 `scale_y_reverse()`，否则会跑到树上边。
- 颜色语义：每个 K 内部按簇平均比例降序配色（cluster 1 = 主簇），跨 K 颜色含义一致；图里不逐个 K 放图例，
  在副标题里写明“Ring colour = ancestry cluster ranked by mean proportion within each K”。
- 195 个叶标签要 15 in 宽、size ≈ 1.9 才在 100% 缩放时读得清（旋转 90° 后标签宽度 = 字号高度，横向不挤）。
- 副标题里的样本数/标记数要**动态生成**（`paste0("Population structure of ", n, ...)`；标记数从 `tree/<PRE>.raw` 的列数 - 6 或
  `ld_prune/<PRE>.pvar` 行数得到）——硬编码数字在换面板后会变成错的（踩过：195 样品图里写着 194）。

## GWAS 候选位点值不值得深挖：先查"是否已被克隆"，再看"能否精细定位"（2026-09 豌豆实测）

挑候选位点前必做两步，否则可能在一个已被克隆的基因上白干几个月。

**第一步：区间内有没有已知基因？** 用**与 GWAS 同坐标系的注释**（GWAS 的 VCF 是对着哪个 ref 调的，
就必须用哪个 ref 的注释；换组装版本位置全错）查区间基因，读 GFF 的蛋白证据属性即可（功能注释没做完也能查）：

```bash
awk -F'\t' '$1=="chr05" && $4>=655000000 && $4<=659000000 && ($3=="mRNA"||$3=="gene")' ann.gff | \
  awk -F'\t' '{match($9,/EvidenceProteinID=[^;]*/); print $4, substr($9,RSTART+18,RLENGTH-18)}'
```

实测：chr05 株高区间里就有 `gibberellin_3beta-dioxygenase`（= 豌豆 *Le* 基因产物）→ 该位点即孟德尔矮化基因，
1997 年已克隆、2025 年 Nature（Feng et al., *Mendel's pea genes*）又补了等位变异 → 不能当新发现。
教训：豆科形态性状的"好看"主效位点（株高 Le、叶型 af、叶腋色素 D、荚型 V、种脐 PPO、花色 A、子叶色 I、
种子皱缩 R=SBE1）基本都已被克隆，报之前逐个性状查一遍文献。

**第二步：算定位分辨率（决定后续能不能缩小区间）。** 对区域内每个 SNP 算与 lead SNP 的 r²：

```r
r2 <- apply(X, 2, function(v) cor(v, lead)^2)   # X = 区域内基因型矩阵（已按位点均值填补）
max(r2[-lead_i])                                # <0.6 且无成块结构 → 可精细定位
```

- r² 出现**连续多 Mb 的块**（如 4.04 Mb / 626 SNP）→ 区间内无重组，**无法细分**，只能做"单倍型类别 vs 表型"，
  别承诺精细定位；这种块常见于自交作物的大结构变异区。
- r² 只有孤立的单点、**最强伙伴 r² 都 <0.6** → 重组充分、**可精细定位**，是更值得投的方向。
- 注意低 MAF 主导的区间（如中位 MAF 0.097）会天然压低成对 r²：即使没有重组热点，一个新位点也可能
  显示"无单倍型块"；判断定位分辨率时同时报区间 MAF 分布 + 背景 LD（随机 ≤10 kb 位点对的中位/均值 r²），
  不要单看 lead SNP 的伙伴数。

**第三步：用第二套独立群体验证（最便宜的假阳性过滤器）。** 如果项目里还有另一个面板：
先用 `bcftools index -s` 对两个 VCF 比对同一染色体长度，确认坐标系一致；再直接到该位置找变异——
**要找的是同一个变异，不是附近的替身**：补齐型 MNP 在两个 VCF 里常常逐字相同（实测 chr04:242929471
在 SRR 与 NUM 两套 VCF 里都是 `GCTGTGTAGGC`→`GCTGCGTAGGC`），能直接用做同一位点验证。
判读标准：① 多性状方向一致（不只是 p<0.05）；② 量级相近（实测两个群体百粒重都约 +13%）；
③ 若某性状在一个群体显著、另一个不显著（如每荚粒数），只能写成"可能存在补偿效应"，不能当结论。
跨群体比较前先确认表型同名（"百粒重"在两套表型表里可能是 `100_dry_seed_weight_g` 与 `trait_24`，
需要通过表型字典映射确认）。

**第亄步：基因型类别的等位基因方向必须回查 VCF（踩过，2026-09-21）**：`plink2 --export A` 的 dosage
到底是 REF 还是 ALT 计数，不要凭印象——用同一位置的 VCF 计数对账（`bcftools query -r pos -f '[%GT\n]'`
统计 0/0、0/1、1/1 的份数，再与 .raw 列里 0/1/2 的份数比）。实测那个位点：VCF 0/0=38、0/1=20、1/1=129，
而 .raw 里 0/1/2 = 129/20/38 → **dosage 计的是 REF**，按 dosage 直接标 "G0 = ref/ref" 会把全部等位基因方向写反
（报告、图、表全错，而 p 值/η² 不会暴露这个错）。正确做法：先算 ALT dosage = 2 - dose，再定类别。

另外：位点可能是**补齐型 MNP**（如 `GCTGTGTAGGC`→`GCTGCGTAGGC`，实际只是第 5 位 T>C，
`bcftools norm -a --atomize` 后为 `chr04:242929475 T>C`）——`pos` 会与真实碱基位置差几个 bp，
引物/探针设计必须用原子化后的坐标。同时注意这类位点常常是 freebayes 的临界调用位点：
实测该位点 37/225 样品无 GT，但它们的覆盖度正常（6–14×）——原因是 GQ<20 被样品级过滤剔除，
不是深度不够；这种位点做标记前应先靶向重定型。

### 好看的单倍型/LD 图（用户要求"好看的图"时的做法）
| 图 | 画法要点 |
|---|---|
| LD 热图 | **必须按位点序号等距排布（rank 坐标）画上三角**；用基因组坐标画会得到稀疏散点+大片空白。刻度标签写上该位点的 Mb 位置（刻度因此不等距，需在图注说明）。下方接基因模型轨道（`geom_rect` 最小宽度 + 标注离 lead SNP 最近的 5 个基因） |
| 局部 LD 衰减 | 以 lead SNP 为中心按距离分箱（0-1/1-2/2-5/…/500-1600 kb）画 median r² + IQR 带 → 直接读"多少 kb 内还能连锁" |
| 基因型类别 | lead SNP 的 0/1/2 三类计数柱（当只有 1 个标签 SNP 时，这就是全部可说的"单倍型"信息，不要硬画网络图） |
| 单倍型结构/网络 | 只有当区域内存在 r²≥0.5 的标签 SNP 集合（≥3 个）时才有意义；否则跳过 |
| 选择信号 | 区域 π（50 kb 窗口）+ Tajima's D（100 kb）+ 基因轨道，标出 lead SNP 与连锁区间 |

绘图统一：英文、Arial、600 dpi PNG + cairo_pdf 矢量；报告里每张图配"怎么算好"的解读段（见 bio-html-report-standards）。

## 选择信号打分：缺数据的窗口必须剔除，否则 top 榜全是空窗口（2026-09 实测）

89,963 个 100 kb 窗口里 38,428 个 Fst=0（Fst 窗口未定义/与 π 窗口不重叠）、9,517 个 TajD=nan
（π 只覆盖 8.8% 位点的时代留下的），原打分公式把这些当成"极端分化"，导致 top30 前 5 名全是
Fst=0.000000 + 正 TajD 的无数据窗口。重算规则：

- 窗口内变异数 **≥10** 才参与打分（Fst 表的 `N_VARIANTS`）；
- 两群体 π 都必须存在且 >0（π 窗口用"窗口内位点 π 之和/位点数"自算，避免再跑一次 vcftools）；
- 两群体 Tajima's D 都必须非 nan；
- 负 Fst 截断为 0（抽样噪声）但保留记录；
- 打分用**秩变换**（rank-based）而非原值加权，避免某一项量纲/极值主导。

## Stage 12: 群体下游统计（π / Tajima's D / Fst / LD 衰减 / ROH）—— 2026-09 无花果三群体实测

分组定下来后（本例 G1/G2/G3 = 79/71/45，按树 + admixture 手切），以下五件套一次跑完（总耗时约 1 h，含 5 个图）：

| 指标 | 命令 | 窗口参数 | 备注 |
|------|------|---------|------|
| π | `vcftools --gzvcf VCF --keep G.txt --window-pi 50000 --window-pi-step 10000` | 50 kb/10 kb | 再跑一次不带窗口的拿全基因组单值 |
| Tajima's D | `vcftools --gzvcf VCF --keep G.txt --TajimaD 50000` | 50 kb | **必须用未做 MAF 过滤的位点集**（见下） |
| Fst | `vcftools --gzvcf VCF --weir-fst-pop A.txt --weir-fst-pop B.txt --fst-window-size 50000 --fst-window-step 10000` | 50 kb/10 kb | 基因组均值在 log 里，窗口值在 `.windowed.weir.fst` |
| ROH | `plink --bfile bed --keep G_fid.txt --homozyg --homozyg-window-snp 50 --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-het 1 --homozyg-window-missing 5` | — | 需 **PLINK 1.9**（plink2 无 --homozyg） |
| LD 衰减 | 自算（见下） | 10 个距离区间到 500 kb | 每染色体、每区间随机抽 20 万对位点算 r²，图上浅点是各染色体、线是组均值 |

实测数字（无花果 195 样本，1,886,125 位点）：π = 0.0024–0.0027；Tajima's D（step1）= +0.89 / +2.19 / +1.98；
Fst = 0.058（G1-G3）～ 0.118（G1-G2）；ROH 总长中位 8.6–19.1 Mb；r² 降到 0.2 的距离 100–250 kb。

### 坑

- **Tajima's D 千万不要在 MAF 过滤后的位点集上算**：MAF≥5% 把低频变异全剔了，D 会被系统抬高（本例 G1 从 +0.89 抬到 +2.89）。
  用 step1（只做 DP/GQ/QUAL 过滤、未做 MAF）那一版，两套都跑、报告里写清楚用的是哪套。
- **plink2 没有 `--homozyg`**（ROH 需要 PLINK 1.9）：`mamba install -n <env> -c bioconda plink` 会因清华镜像超时失败，
  改从 conda.anaconda.org 直接拉 tarball 手工装（与 admixture 同一套办法）：
  `curl -sSL -o plink.tar.bz2 https://conda.anaconda.org/bioconda/linux-64/plink-1.90b7.7-h18e278d_1.tar.bz2 && tar -xjf ... && cp bin/plink <env>/bin/`。
- **PLINK 1.9 的 `--keep` 必须两列**（FID IID）；配合 `--double-id` 时应写 `样品名\t样品名`，单列会报
  `Line 1 of --keep file has fewer tokens than expected`（跑完才发现 .hom.indiv 没生成）。
- PLINK 1.9 读非标准染色体名（chr01a）要加 `--allow-extra-chr`；转 bed 一次（`--make-bed`）后各组直接 `--bfile` 跑，比每次读 VCF 快得多。
- Fst/π 的 Manhattan 图用 `manhattan-alternating-chr-colors` skill 的样式（奇染色体 #0072B2 / 偶 #E69F00、纯白底、前 5%/1% 虚线阈值）。
- 分析集里若有重复样本（同一基因型多次），π/Tajima's D/ROH 都带权重偏差——建议同时跑一份去重集做敏感性分析。

## Stage 4: freebayes 联合 calling（脚本 + 实测 2026-09-17）

环境 `fb_calling`：freebayes 1.3.10 + freebayes-parallel + bcftools 1.23.1 + vcftools 0.1.17。

**硬性约定（用户 2026-09-17）：calling 一律用 freebayes，不用 bcftools mpileup/call。** bcftools 只允许出现在 calling 之后的 VCF 工具链（concat/norm/index/stats），不用于基因型判定。

**通用脚本**（拷到新项目后只需改开头的 `P=` / `V=` 路径）：

| 模板 | 作用 |
|------|------|
| `templates/make_regions.py` | .fai → 100 kb 分块区间 `regions/chunks/<chr>.regions` |
| `templates/make_work_units.py` | 把分块**轮转均分**成 N 个工作单元（N = 想要并发数，惯例 40） |
| `templates/freebayes_call_unit.sh` | 跑一个工作单元（单元内顺序、块级独立 VCF、已存在则 SKIP） |
| `templates/freebayes_call_unit_array.sub` | SGE 数组 N 任务 × `-pe smp 6` → 并发恰好 ~N（40 进程 / 240 槽） |

```bash
# 1) 分块 + 工作单元
python3 make_regions.py 100000 && python3 make_work_units.py 40
# 2) 每块的命令（核心参数）
freebayes -f ref/genome.fasta -r chr01:1-100000 \
    -q 20 -m 20 --use-mapping-quality --genotype-qualities --ploidy 2 \
    --min-alternate-count 2 --min-alternate-fraction 0.05 \
    --bam-list lists/bamlist.txt -v vcf/parts_raw/chr01-1-100000.vcf
# 3) 提交：40 单元 × -pe smp 6 = 240 槽 → 恰好 40 个 freebayes 进程（freebayes 单线程，并发靠进程数）
qsub -q all.q@CNS2 -pe smp 6 -l h_vmem=10G -cwd -o logs/fb_array -e logs/fb_array -t 1-40 scripts/freebayes_call_unit_array.sh
```

**为什么用「用户指定并发数 = 工作单元数」的排法**：freebayes 单进程单线程，所以要拿到的并行度只能靠进程数；“任务数 × `-pe smp k`”能让 SGE 精确地把你想要的进程数跑起来（同时留出核给其他人）。每单元约 254 块 × 100 kb，实测约 7-8 h/单元。

**自动接力**：`watch_and_call.sh` 模式（见 dadou 项目 `05.variants/scripts/watch_and_call.sh`）—— 每 30 min 检查“BAM 数到齐 + 映射数组已退出” → `samtools quickcheck` → 自动 qsub calling；用 `fb_submitted.txt` 做幂等，已提交过就直接退出。适合用 cron 驱动，不用常驻进程。

**实测（大豆 YD1 参考 1,015,170,271 bp，346 样本 × 平均 20.8 Gb clean）**

| 测量 | 值 |
|------|-----|
| 1 样本，1 Mb | 5.14 s，RSS 185 MB（参考不进常驻内存，块小则内存小） |
| 1 样本，100 kb | 1.13–1.42 s |
| 同样本复制 20 份，100 kb（代理） | 6.08 s，RSS 225 MB |
| **真实 346 样本，100 kb** | **中位 32 min**（p25 27 min / p75 38 min / 最长 8.1 h），RSS ~3 GB/进程，中位 4,766 个原始变异/块 |

**重大教训：不要用“同一 BAM 复制 N 份”当多样本代理去外推耗时。** 复制样本 20 份测出 6 s/块，按线性外推到 346 份得 ~100 s/块；真实 346 个不同样本实测 **1,922 s/块（32 min，慢 19 倍）**，因为 freebayes 的联合基因型代价随样本数**超线性**增长（每样本每候选等位都要算似然，样本越多候选等位/单倍型组合越多）。真实规模：1.015 Gb 基因组 / 10,155 块 × 32 min ÷ 40 进程 ÷ 60 ≈ **135 h ≈ 5.6 天**。

推论：全样本联合 freebayes 在“数百样本 × 20×”规模下是**天级**任务，不是小时级。要提速只有三条路：① 堆进程数（单进程只吃 1 核、~3 GB，本项 40 进程才 67/256 核，节点很空）；② 缩小分块大小以降低单块耗时方差（p75/p25 = 1.4，最长块 8 h 是拖后腿的）；③ 放弃联合改用单样本 call + merge（快 ~10× 但丢失联合基因型优势）。想控制时间就先实测 1 块再决定并发。

**坑**

- `freebayes` 拒绝重复 BAM 路径：`Error: Duplicate bam file` → 性能测试用不同名符号链接（`.bam` + `.bam.bai` 都要链）。
- 用 `xargs -P` 调 bash 函数时，函数里用到的所有变量都必须 `export`（`REF BAMLIST LOG OUT` 漏一个就会重定向到空路径 → `Permission denied` 后整块 FAIL）。
- 块边界可能有重复位点 → 拼接时 `bcftools norm -d exact`。
- 细胞器（Pt/Mt）深度极高（本例 Pt 2 万×），会明显拉长单块耗时；建议单列、不并入核基因组人群分析。
- 块级 VCF 会非常大：freebayes 输出带 per-genotype **GL** 字段，346 样本 × 100 kb 块 ≈ **83 MB/块**（全量 10,155 块 ≈ 845 GB）。拼接后的终稿建议 `bcftools view -G`（或 `-i` 选字段）丢掉 GL，只留 GT:GQ:DP:AD 等。
- **`qstat` 队列为空时无输出且 exit 0** —— 看门狗/监控脚本绝不能把“空输出”当失败（反之亦然），只能看退出码。曾经因这个误判把自动提交阻断了 15 h。
- 看门狗提交前先 `samtools quickcheck` 所有 BAM，`bad>0` 就本轮不提交（正在写、未封口的 BAM 能被 quickcheck 揪出来）。
- 自动接力用 cron 驱动（每 30 min 调一次幂等脚本）比常驻 watchdog 进程可靠：会话断了也不影响，提交后写一个 `fb_submitted.txt` 做幂等标记。
- 与 `freebayes-parallel` 相比，SGE 数组 + 工作单元的方案槽位/内存可控、可断点续跑（`freebayes-parallel` 需 `parallel`+vcflib 全家桶，且中途挂了要从头来）。

**Citations（calling 步骤）**

- Garrison E, Marth G. *Haplotype-based variant detection from short-read sequencing.* arXiv:1207.3907 (2012). [freebayes]
- Danecek P, et al. *Twelve years of SAMtools and BCFtools.* GigaScience 10(2):giab008 (2021). [bcftools]
- Li H. *minibwa* (BWA-MEM successor), GitHub lh3/minibwa. [比对]

BibTeX 与 Methods 模板句：

```bibtex
@misc{garrison2012freebayes, title={Haplotype-based variant detection from short-read sequencing},
  author={Garrison, Erik and Marth, Gabor}, year={2012}, eprint={1207.3907}, archivePrefix={arXiv}}
@article{danecek2021bcftools, title={Twelve years of SAMtools and BCFtools}, journal={GigaScience},
  volume={10}, number={2}, pages={giab008}, year={2021}, author={Danecek, Petr and others}}
```

> Methods：Reads were aligned to the reference genome with minibwa (`-x sr`); joint variant calling
> was performed with freebayes (v1.3.10) in 100-kb windows using mapping quality ≥20 and base
> quality ≥20, followed by bcftools normalisation and hard filtering (QUAL ≥20, GQ ≥20,
> DP 4–1000, missing rate ≤0.2, MAF ≥0.05, biallelic SNPs).

## Stage 4: freebayes 联合 calling（用户硬性要求：重测序 SNP calling 一律优先 freebayes，不用 bcftools call）

模板脚本（成熟、已用于 346 样品的 dadou 面板）：`/media/nfs1/hermes/project/dadou/req/05.variants/scripts/`
下 `call_freebayes_chr.sh`（单染色体分块）、`call_freebayes_array.sh`（SGE 数组，每任务 4 进程）、`concat_freebayes.sh`。
在新项目里照抄，只改 ROOT/REF/BAMLIST 路径。

```bash
freebayes -f $REF -r $CHUNK \
  -q 20 -m 20 --use-mapping-quality --genotype-qualities --ploidy 2 \
  --min-alternate-count 2 --min-alternate-fraction 0.05 \
  --bam-list $BAMLIST -v $OUT/${CHUNK//:/-}.vcf.tmp && mv ... .vcf
```

- 切块 100 kb（与 dadou 一致）；块独立、写 `.tmp` 再 `mv` → 失败不会留下“看起来完成”的文件；
  已存在的 `.vcf` 直接跳过 → 断点续跑。并行粒度：不用 SGE 时用**一个全局块池 `xargs -P 40`**，
  比“每染色体 4 进程”负载更均衡（避免大染色体拖尾）。
- 参考只喂核染色体，**不要包含 Pt/Mt**（叶片 DNA 叶绿体深度 1e4×，194 样品叠加是百万级，纯浪费）。
- `--genotype-qualities` 就是输出 FORMAT/GQ（FORMAT = `GT:GQ:DP:AD:RO:QR:AO:QA:GL`），所以 freebayes 输出照样能按 GQ≥20 过滤。
- 收尾：按坐标顺序 `bcftools concat -a -f` 每染色体（去掉块边界重复），再 `bcftools norm -d exact` 去重，最后全基因组拼接。

### 耗时特征（必读，决定要不要调参）

freebayes 耗时**正比于候选位点数 × 样品数**，不是正比于数据量。194 样品 100 kb 块实测：

| 块内记录数 | 单块 wall |
|---|---|
| 730 | 90 s（0.12 s/位点） |
| 5575 | 2948 s（0.53 s/位点） |
| 7758 | 3969 s |

位点密度极度不均：chr01a:1-100k 只有 730 个，chr01a:2.0-2.1 Mb 有 7758 个（1 个变异/13 bp）——
其中 68% 是 AF<0.01 的单例（194 样品里只出现在 1-2 个样品），AF≥0.05 的只有 1223 个。
这些单例会被下游 MAF≥0.05 过滤掉，但 freebayes 为它们做全额枚举，直接决定总时长。
3227 块 × 40 路：全程 **60-100 小时（2.5-4 天）**。想压时间就抬 `--min-alternate-fraction` 到 0.2、
`--min-alternate-count` 到 3，加 `--use-best-n-alleles 4`（预计快 3-5 倍，代价是丢真稀有变异）；
用户明确要求“与 dadou 同一套参数”时就不要调。
体检查过密度异常不是深度造成：那些区块单样品仅 26-52×，>200× 碱基占 0.1-0.4%。

## Stage 4b: VCF 过滤参数

| 参数 | 值 |
|------|-----|
| QUAL | ≥20 |
| GQ | ≥20 |
| DP | [4, 1000] |
| F_MISSING | ≤0.2 |
| MAF | ≥0.05 |
| Alleles | bi-allelic SNP only |

## Stage 11: 可视化（五图模板）

**Fig1:** 环形树（整枝着色，邻近向上传播）
**Fig2:** 环形树 + K=2 Admixture 外圈环（Q矩阵需对齐树尖序）
**Fig3:** PCA（PC1/2 + PC1/3，A/B着色）
**Fig4:** Admixture K=2-4（组内按祖先比例排序）
**FigCV:** 交叉验证误差曲线

### 环形树 + K=2..N admixture 外圈：不要依赖 ggtreeExtra（2026-09 实测）

本机 R 4.5.2 + ggtreeExtra 组合下 `geom_fruit(geom=geom_col/geom_bar)` 稳定报错
（`replacement has 1 row, data has 0`；换 x 为 numeric 后改报 `missing value where TRUE/FALSE needed`），
`ggtreeExtra` 的 fruit 层与 circular 布局的 data 合并行为不可控。**改用自算极坐标布局**：
模板 `templates/circular_tree_admixture_rings.R`（可直接跑，改开头路径即可）。

要点（踩过的坑）：

- 叶角 `theta_i = 2*pi*(i-1)/n`（i 按 ladderize 后的 tip 顺序）；内部节点角度 = 后代叶索引 min/max 的中点 → 天然避开 0/2π 环绕。
- 枝长做幂变换 `R = (depth/max)^0.45`：NJ 树内部分支极短（叶枝 ~10 vs 内部 ~0.1），不缩放就是不可读的“星爆”，图中注明变换即可。
- 每条边画两段：父角度上的径向线 + 子角度处半径上的圆弧。**`geom_path` 必须 `aes(group=seg)`**，否则所有段按行序连成一条折线（表现为“只有放射线、没有弧”）。
- 外圈：每个样本一个角度槽（半宽 = 2π/n × 0.45），槽内按祖先比例沿**径向**堆叠 K 段（annular sector 多边形）。
- 样本标签放在**最外圈之外**（`r = Router + 0.035`），翻转规则 `phi>90 & phi<270 → angle=phi-180, hjust=1`（保证文字不倒立、且向外延伸不压环）；194 个标签要 15 in 画布、size 2.0 才读得清（PNG 9000×9000 px，vision 预览需先缩到 ≤2600 px 再传）。
- 颜色语义：每个 K 内部按簇平均比例降序配色（cluster 1 = 主簇），因此不同 K 的同一颜色含义一致；图内注明“Ring colour = ancestry cluster ranked by mean proportion within each K”，并在左下角放一个“Rings, inner to outer: K = 2, 3, ... , N”的注记，比给每圈加引线更干净。
- 每个 K 一个环：先画内到外，环间距 0.012、环宽 0.075（半径单位）。

若 K 的 CV 误差一直单调下降（无极小值），说明分组是连续分化而非离散群体：
按用户流程“先用树切组，再用 admixture/PCA 交叉验证”，不要硬报一个“最优 K”。

## 关键避坑

1. plink2 LD: 必须 `--set-all-var-ids @:# --rm-dup`
2. ADMIXTURE: 整数染色体，.Q 在 CWD
3. ANGSD `-fold`: 仅在 realSFS
4. **ANGSD 2D SFS 不可行（>140样本）：** 即使串行也常卡死或 OOM。改用 vcftools `--weir-fst-pop --fst-window-size 50000` 在已过滤 VCF 上算 Fst，分钟级完成
5. PCA 特征值: 除总和
6. Q矩阵对齐: `Q[match(tree$tip.label, fam_order), ]`
7. **`xargs -P N ... bash -c 'myfunc {}'` 里函数用到的每个变量都必须 `export`**（2026-09 一次会话踩两次）：featureCounts 脚本漏 `BASE`/`HAP` → 日志路径塌成 `/logs/S.v2_.s0.fc.log`，重定向失败直接打断命令，197 样本全 FAIL 且 counts 全空；mapping 脚本漏 `SAMDIR` → SAM 路径塌成 `/DRR101541_whg027a.sam`，196/197 FAIL_ALIGN。排查口诀：日志里出现**以 `/` 开头的意外路径**（`/logs/...`、`/DRR...sam`）就是变量没导出。
8. **plink2 没有 `--distance`**（本机版本，exit 8 `Unrecognized flag`）：建树改 `plink2 --export A`（0/1/2 矩阵）→ R 读 `.raw`、按位点均值填补缺失（本例 0.11%）、`dist()` 欧氏距离、`ape::nj()`；不要指望在 plink2 里一步出距离矩阵。
9. **plink2 的 `.eigenvec` 带表头**（`#IID PC1 PC2 ...`）：`read.table(header=TRUE)` + 按列序号 `ev[[1]]`/`ev[[2]]` 取，否则表头被当数据、列整体错位（`ev$V2` 变成 PC1 数值），后续 `left_join(by="id")` 会因型别不符报 vctrs 错误。
10. **拼两个 ggplot 必须 `library(patchwork)`**（否则 `+` 报 `Can't add <ggplot> to a <ggplot> object`）；精简脚本时最容易漏掉。
7. 枝干着色: 邻近向上传播，非多数投票
8. Admixture 成分映射: C1可能=B需验证
9. 长任务: nohup > log 2>&1 & 避免包装器超时
10. PopLDdecay: 全路径 `/media/nfs1/hermes/miniforge3/bin/PopLDdecay`，`comment.char=""` 读 .stat.gz，不要 `-SubPop`（segfault）
11. bcftools/PopLDdecay 不在 pop_genetics env: 需 `conda activate fb_calling` 或用全路径
12. **minibwa index 不接受软链接路径**：对 symlink 报 `[E::main_index] failed to read the genome FASTA`，必须用真实路径（或先 cp 再 index）。`minibwa map` 无此问题
13. 交付前必查 FASTQ 数量 vs 清单：无花果 req 204 个编号里 6 个（WHG068/076/080/110/124/130）完全没有数据，WHG029 分 2 个 lane（WHG029-1/-2）→ 实际 199 个文库 / 198 个样品，`run_status.tsv` 行数要按文库数核对
14. **minibwa 加线程收益递减（2026-09 实测）**：同一样本 6 线程 8.6 Mbases/s → 40 线程 17.5 Mbases/s（仅 2×，且每任务只吃满 28-30 核）。故“少线程 × 多并行”聚合吞吐通常更高（8 线程 × 32 并行 > 40 线程 × 6 并行）。量级参考：20.8 Gb clean/样本，40 线程/任务、6 并行 → 单样本 ~20 min 比对 + 4-5 min samtools sort；每任务 RSS 5-6.7 GB。
15. **只走 CNS2 的提交语法**：`qsub -q all.q@CNS2 ...`（`qconf -sq all.q` 的 hostlist 是 `@allhosts CNS2`，@allhosts 只含 CNS1）。SGE 按槽数自动打包：-pe smp 40 时 256 槽的 CNS2 同时跑 6 个任务。`-o/-e` 目录须先建。
16. **map 脚本要幂等**：跑前判断 `[ -s bam/$S.sorted.bam ] && [ -s ...bai ]` 就 exit 0，数组重提/续跑才不会重算已完成样本。
17. **plink .bed 文件头3字节**（0x6C 0x1B + 模式位0x01），numpy 解析必须 `raw[3:]`；从 `raw[2:]` 切会整体错位1字节 → 样本身份错乱（组内距离≈组间、93%缺失样本被算成93%有值），2026-09-22 豌豆 IBS 矩阵踩过。
18. **vcftools 0.1.17 无 `--fst-window`**（Unknown option）：只出每SNP `.weir.fst`（CHROM POS WEIR_AND_COCKERHAM_FST），100kb 窗口自行聚合（bin 求均值）；`--keep` 也不能加在 Fst 步骤上（会把另一组样本删光）。
19. **PopLDdecay `.stat.gz` 表头是 `#Dist Mean_r^2 Mean_D' ...`**：`^Dist` 锚定匹配不到 `#Dist`，`"r2"` 匹配不到 `r^2` → 第1列做距离、`grepl("Mean_r", fixed=TRUE)` 取 r²；共500001行（逐1bp），画图前截 `<100 bp` 并 log 分箱。
20. **监控/接力脚本的 token 粘连**：`sleep60`/`shift2` 少一个空格 → 前者循环空转秒开、后者内层 while 死循环白烧 CPU（2026-09-22 一晚踩两次）；`bash -n` 查不出这类错，新脚本要拿空集冒烟测一遍。
21) **FastTree 显示层坑（与结构条拼图时）**：VCF→fasta 把杂合写成 IUPAC（R/Y）、缺失写 X，FastTree 全部忽略（日志 “ONLY76.7% nucleotide”）→ 枝长虚高（同种内0.85–2.1 subs/site）、tips 高低参差“枝飘在上边”。改用 **IBS 距离 + UPGMA（hclust average）**：叶统一落在 y=0 基线正贴柱条，无参差、无离群拉伸；FastTree 只作存档不动。
22) **`pkill -f` 会杀掉自己的 shell**：若 pattern 字符串出现在待执行命令的命令行里（命令文本本身含该脚本名/路径，如 `pkill -f '40_watch_and_aggregate'` 出现在同一条复合命令中），pkill 连自己的 bash 一起匹配 → 后面的 `rm`/`ls` 全部不执行且无报错（2026-09-23 连踩三次，表现为“清理命令输出为空、文件却还在”）。改用 `kill -9 <明确PID>`；只有 pattern 不可能匹配自身调用行时才可用 pkill。
23) **报告生成器里的组数/列号/定性结论必须构建时动态生成**：写死（n=65/13/147、第1–65列、"G2 最独特"、"三组相差不到3%"、"n=13 小样本偏倚"）在重分组后会静默输出与数据相反的结论（2026-09-23 重分组后实测 π 差20%、Fst 排序反转、杂合度最高的换成另一组）。组名/列界从 lists 或列序表读，"谁最高/谁最独特/是否同量级"用 max/min + 条件分支生成，跑完后 grep 旧口径关键词（n=13、同量级、最独特）必须为0。
24) **data.table fwrite 输出窗口坐标可能是科学记数法**（1.68e+08 / 5.01e+08，豌豆五个窗口表全部中招）：Python `int('1.68e+08')` 直接 ValueError 拖垮整个构建；窗口坐标一律 `int(float(x))`，展示层同理。
25) **复合命令里的词粘连（同第20条一类，新证据）**：`rm -rf A logsB` 少空格写成 `logs00.report`/`-rf00.report` → rm 报 `invalid option -- '0'` 静默删不掉目标；`cmd1 || cmd2; cmd3` 结构若 cmd1 已成功，cmd3 仍会把**同一个脚本再跑一遍**（两条链并发写同一文件=数据竞态，2026-09-23 G1 子集 VCF 被两个 bcftools 同写）。提交后台任务前逐词查空格，一条链一条命令。
26) **曼哈顿图 Bonferroni 分母必须是标记数，不是基因组长度**（2026-09 豌豆踩过）：`-log10(0.05/sum(染色体最大POS))` 把线画到 y=10.9 而非 `0.05/6,446,649` 的 y=8.11，导致真显著点（1279/215/186 个里的绝大多数，包括 chr04 lead y=9.53）全部落在线下，图看起来“没有显著位点”。红线位置必须与报告图注写的数字同源；另：x 轴刻度只标 chr1–N（map 里常有未定位 contig/MT，会生成 chrh1tg… 挤压标签）；显著点要换深色大点+标注 top SNP，缩到报告显示宽度（1350px）后孤立小点（size0.25/alpha0.4）肉眼不可见；图注只能写图里真实存在的元素（曾把不存在的 QQ 面板写进图注）。

## 项目实例：大豆 dadou（2026-09）

- 路径 `/media/nfs1/hermes/project/dadou/req/`：`01.rawdata`(3.8 T) → `02.cleandata`(3.1 T) → `03.fastp_report` → `04.mapping` → `05.variants`(计划)。
- 346 文库 = 346 样品，DNBSEQ PE150，平均 clean 22.9 Gb/样本（总 7876 Gb）；参考 YD1 1,015,170,271 bp（chr01–chr20 + Pt + Mt1–Mt5）。
- 比对：minibwa -x sr，40 线程/任务，`-q all.q@CNS2`，6 并行 → 345 样本约 23 h；BAM 合计 ~1.9 TB。进度一行命令：`bash 04.mapping/scripts/status.sh`。
- 变异检出（用户指定 freebayes，脚本已就绪并单样本验证）：`05.variants/README_plan.md` + `scripts/{make_regions.py,call_freebayes_chr.sh,call_freebayes_array.sh,concat_freebayes.sh,filter_vcf.sh}`；100 kb 分块（核 10,155 + 细胞器 10），SGE 数组 26 任务 `-q all.q@CNS2`；预估 3-5 h。
- 比对实测（2026-09-17）：40 线程/任务 × 6 并行 = 单样本 38.6-42 min（minibwa 加线程几乎无增益，6 线程单任务当时也只要 47 min）→ 345 样本 ≈ 38 h；若改 8 线程 × 32 并行预计 8-12 h。提交时用 `-q all.q@CNS2` 限节点。

---

# VCF missing-genotype encoding: "." vs "./." silently breaks vcftools (2026-09-21)

## Symptom

`vcftools --site-pi` (and `--weir-fst-pop`, `--window-pi`) returns results for only a small fraction
of sites with no error: on a 6,446,649-site pea VCF, pi was emitted for **566,740 sites (8.8%)** while
vcftools logged `After filtering, kept 6446649 out of a possible 6446649 Sites`.

## Root cause

The VCF mixes two notations for a missing genotype — `.` (haploid style) and `./.` (diploid style) —
with **no true haploid allele calls** anywhere. vcftools 0.1.17 silently drops every site containing a
`.` call. Verified on chr01:1-5 Mb: the set of sites with pi output was **exactly** the set of sites
containing no `.` call (611 of 5,278 sites, identical sets). A toy VCF with uniform `./.` returns pi
for 100% of sites, so the tool itself is fine.

Diagnose in one region:

```bash
bcftools view -r chr01:1-5000000 -Oz -o test.vcf.gz in.vcf.gz && bcftools index test.vcf.gz
bcftools query -f '[%GT\t]\n' test.vcf.gz | tr '\t' '\n' | sort | uniq -c | sort -nr | head
vcftools --gzvcf test.vcf.gz --site-pi --out test      # rows << sites => this bug
```

## Fix

Normalize ploidy — genotype content is unchanged (`.` and `./.` both mean missing), only the encoding
becomes consistent. Note the plugin syntax: general options BEFORE the file, plugin options after `--`.

```bash
bcftools +fixploidy -Oz -o out.fixploidy.vcf.gz in.vcf.gz -- -f 2
bcftools index -c out.fixploidy.vcf.gz     # CSI: chromosomes > 536 Mb exceed the BAI 2^29 limit
```

After the fix the same region returned pi for all 5,278 sites. Keep the master VCF untouched and run
all diversity statistics (pi, Tajima's D, Fst, LD) on the normalized copy.

## What is NOT affected

plink2 read the same VCF without dropping anything (`genotype.bim` kept all 6,446,649 variants), so
GWAS / PCA / admixture results are unaffected. Only vcftools-based site-level statistics are hit.

## Related: the sample with the pathological missing rate

In a 20-Mb test region, 34,467 of 37,947 sites (90.8%) carried at least one `.` call, and the worst
sample (SRR28225809) had `.` at 33,796 sites (~89% of its calls) — that same sample had already been
dropped from the GWAS genotype set for excessive missingness. When one sample's call rate collapses,
check whether it is also the sample driving the encoding mix.


---

# Citations（软件与方法原文；发表论文时请一并引用）

本流程用到的每一个软件/方法对应的原始文献。发表论文时 Methods 段落逐项引用。

## 质控

- Chen S, Zhou Y, Aydin A, et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics 34(17):i884–i890 (2018). PMID: 29506976
  - 工具：fastp（去接头、质量过滤）

## 比对

- Li H, Homer N. *Fast genomic read alignment with minibwa.* arXiv:2606.15357 (2026). DOI: 10.48550/arXiv.2606.15357
  - 工具：minibwa（`-x sr` 短读模式；BWA-MEM 后继，约 4× 于 BWA-MEM）
- Li H. *Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.* arXiv:1303.3997 (2013). [minibwa 的直接前身]
- Danecek P, et al. *Twelve years of SAMtools and BCFtools.* GigaScience 10(2):giab008 (2021). PMID: 33504313
  - 工具：samtools（sort / index / faidx / view / query）

## 变异检测（freebayes 联合 calling）

- Garrison E, Marth G. *Haplotype-based variant detection from short-read sequencing.* arXiv:1207.3907 (2012).
  - 工具：freebayes（联合基因型判定；`--genotype-qualities` 输出 FORMAT/GQ）
- Danecek P, et al. *Twelve years of SAMtools and BCFtools.* GigaScience 10(2):giab008 (2021).
  - 工具：bcftools（concat / norm / index / +fixploidy —— 只用于 calling 之后的 VCF 工具链，不做基因型判定）

## VCF 过滤与统计

- Danecek P, Auton A, Abecasis G, et al. *The variant call format and VCFtools.* Bioinformatics 27(15):2156–2158 (2011). PMID: 21669827
  - 工具：vcftools（`--weir-fst-pop` Fst、`--window-pi` π、`--TajimaD`）

## 群体结构

- Chang CC, Chow CC, Tellier LC, et al. *Second-generation PLINK: a toolset for the analysis of genome-wide association studies.* Bioinformatics 31(5):599–604 (2015). PMID: 25428804
  - 工具：plink2（PCA / LD pruning / KING 亲缘 / dosage 导出 / --make-king-table / --het）
- Purcell S, Neale B, Mills S. *PLINK 1.9.* (2019). GitHub: chrchang/plink-ng
  - 工具：PLINK 1.9（ROH `--homozyg`；plink2 无此选项）
- Alexander DH, Novembre J, Lange K. *Fast likelihood-based estimation of ancestry in ADMIXTURE.* Genome Research 19(7):1175–1180 (2009). PMID: 19562099
  - 工具：ADMIXTURE（模型聚类 + CV 误差）
- Price MN, Dehal PS, Arneberg O, et al. *FastTree 2 — Approximately Maximum-Likelihood Trees for Large Alignments.* PLoS ONE 7(12):e9490 (2010). PMID: 22228091
  - 工具：FastTree2（邻接树存档）
- Paradis E, Schliep K. *ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.* Bioinformatics 34(7):1643–1644 (2018). PMID: 29198420
  - 方法：IBS 距离 + `ape::nj()` / `hclust` (UPGMA) 画树（本流程默认替代 FastTree 出图）
- Korneliussen TS, Albrechtsen A, Nielsen R. *ANGSD: Analysis of Next Generation Sequencing Data.* BMC Bioinformatics 15:204 (2014). PMID: 25420514
  - 工具：ANGSD（π / Tajima's D / Fst 的 ANGSD 路线；>140 样本的 2D SFS 不可行时改走 vcftools）

## 连锁不平衡

- Zhang C, Dong SS, Xu JY, He WM, Yang TL. *PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files.* Bioinformatics 35(10):1786–1788 (2019). PMID: 30321304. DOI: 10.1093/bioinformatics/bty875
  - 工具：PopLDdecay（LD 衰减；本流程对 `.stat.gz` 做 log 分箱 + 池化均值 + nls 拟合）

## 物种 / 混样核验

- Baker E, Wang B. *A non-parametric method for the detection of contamination in sequencing data.* / sourmash 文档
  - 工具：sourmash（k=31, scaled=5000 sketch + compare，用于混样/错样快速判定；可替换为 Mash）

## 可视化

- Yu G, Smith DE, Hu H, et al. *ggtree: an R package for visualizing phylogenetic trees and other graphical structures.* Molecular Plant 10(6):826–828 (2017). PMID: 28900387
  - 工具：ggtree（树结构；环形树+admixture 外圈用自算极坐标布局，不依赖 ggtreeExtra）
- Wickham H. *ggplot2: Elegant Graphics for Data Analysis.* Springer (2016).
  - 工具：ggplot2（全部出版级图件；字体 Arial、600 dpi PNG + cairo_pdf 矢量）
- Wilkerson AD. *patchwork.* R package（拼图，缺 `library(patchwork)` 会报 `Can't add <ggplot> to a <ggplot> object`）
- Yu G, et al. *ggnewscale: An R package for elegant gradient visualization in ggtree.* (2023).
  - 说明：仅在需要连续色标时引用

## BibTeX（可直接粘贴）

```bibtex
@article{chen2018fastp,
  title={fastp: an ultra-fast all-in-one FASTQ preprocessor},
  author={Chen, Shubin and Zhou, Yanfang and Aydin, Al and Farcas, Adrian and
          Le, Haicheng and Lo, Andy and Zhang, Naxun and Tool, Open Science},
  journal={Bioinformatics}, volume={34}, number={17}, pages={i884--i890}, year={2018}
}
@misc{li2026minibwa,
  title={Fast genomic read alignment with minibwa},
  author={Li, Heng and Homer, Nils}, year={2026}, eprint={2606.15357},
  archivePrefix={arXiv}, primaryClass={q-bio.GN},
  doi={10.48550/arXiv.2606.15357}
}
@misc{li2013bwamem,
  title={Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM},
  author={Li, Heng}, year={2013}, eprint={1303.3997}, archivePrefix={arXiv}
}
@article{danecek2021bcftools,
  title={Twelve years of SAMtools and BCFtools},
  author={Danecek, Petr and Abeysinghe, Sudruga and Twardziok, Petr and
          Li, Heng and Marth, Gabor and Abecasis, Goncal and Durbin, Alistair},
  journal={GigaScience}, volume={10}, number={2}, pages={giab008}, year={2021}
}
@misc{garrison2012freebayes,
  title={Haplotype-based variant detection from short-read sequencing},
  author={Garrison, Erik and Marth, Gabor}, year={2012}, eprint={1207.3907},
  archivePrefix={arXiv}
}
@article{danecek2011vcftools,
  title={The variant call format and VCFtools},
  author={Danecek, Petr and Auton, Adam and Abecasis, Goncal and Alkan, Caner and
          Baele, Laurent and Barker, Laura and Berens, Florian and Bethune, Daniel and
          Chan, Wenyuan and Chen, Chia-Ling and Cooper, Cheryl and Cornforth, Darren and
          Cui, Heng and Fledel, Susan and Good, Sean and Gram, Gustav and Guo, B->
          and Hinds, David and Hunt, Laura and Kang, Ines and Keys, Jonathan and
          Kraft, Andreas and Kronenberg, Zohar and Lam, Kenneth and Larkin, Sean and
          Lees, James and Lee, Jeffrey and Luo, Ruibang and McCarthy, Sean and
          Meng, Qi and Nielsen, Rasmus and Molla, Akiba and Ray, A., and
          Scudie, Giulio and Tan, Adrian and Tan, Jianshou and Valouev, Ivan and
          Venter, J Craig and Zhou, Shuai},
  journal={Bioinformatics}, volume={27}, number={15}, pages={2156--2158}, year={2011}
}
@article{chang2015plink2,
  title={Second-generation PLINK: a toolset for the analysis of genome-wide association studies},
  author={Chang, Chang C. and Chow, Carson C. and Tellier, Laurent C. and
          Vattikuti, Shuyang and Purcell, Shaun and Lee, James J. and
          Marioni, Giulia E. and Altshuler, David M. and Park, Won Jin and
          Won, Hong-Hee and Miller, Mark B. and Ng, Pui Chun and
          Graham, Sarah B. and Schrijver, Gerrit and Sabatti, Chiara and
          Stein, Lawrence D.},
  journal={Bioinformatics}, volume={31}, number={5}, pages={599--604}, year={2015}
}
@article{alexander2009admixture,
  title={Fast likelihood-based estimation of ancestry in ADMIXTURE},
  author={Alexander, David H. and Novembre, Joseph and Lange, Kenneth},
  journal={Genome Research}, volume={19}, number={7}, pages={1175--1180}, year={2009}
}
@article{price2010fasttree,
  title={FastTree 2 --- Approximately Maximum-Likelihood Trees for Large Alignments},
  author={Price, Michael N. and Dehal, Paramjit S. and Arneberg, Olof and
          Lanfear, Robert and Reich, Matthew B. and Redding, David and
          DePristo, Mark A. and Andersen, Torsten F. and Siepel, Lasse},
  journal={PLoS ONE}, volume={7}, number={12}, pages={e9490}, year={2010}
}
@article{paradis2018ape,
  title={ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R},
  author={Paradis, Emmanuel and Schliep, Kim}, journal={Bioinformatics},
  volume={34}, number={7}, pages={1643--1644}, year={2018}
}
@article{korneliussen2014angsd,
  title={ANGSD: Analysis of Next Generation Sequencing Data},
  author={Korneliussen, Thorfinn S. and Albrechtsen, Anders and Nielsen, Rasmus},
  journal={BMC Bioinformatics}, volume={15}, pages={204}, year={2014}
}
@article{zhang2019poplddecay,
  title={PopLDdecay: a fast and effective tool for linkage disequilibrium decay
         analysis based on variant call format files},
  author={Zhang, Chi and Dong, Shan-Shan and Xu, Jun-Yang and He, Wei-Ming and Yang, Tie-Lin},
  journal={Bioinformatics}, volume={35}, number={10}, pages={1786--1788}, year={2019},
  doi={10.1093/bioinformatics/bty875}
}
@article{yu2017ggtree,
  title={ggtree: an R package for visualizing phylogenetic trees and other graphical structures},
  author={Yu, Guangchuang and Smith, De-Yuan and Hu, He and Yang, Guolin and
          Goh, Chao Qi and Swofford, Robert L.}, journal={Molecular Plant},
  volume={10}, number={6}, pages={826--828}, year={2017}
}
@book{wickham2016ggplot2,
  title={ggplot2: Elegant Graphics for Data Analysis},
  author={Wickham, Hadley}, publisher={Springer}, year={2016}
}
```

## Methods 模板句（英文，可直接改写进论文）

> Raw sequencing reads were quality trimmed with fastp (Chen et al., 2018) and aligned to the
> reference genome with minibwa (Li & Homer, 2026) using the `-x sr` short-read mode. Joint variant
> calling was performed with freebayes (Garrison & Marth, 2012) in 100-kb windows with mapping quality
> >= 20 and base quality >= 20, using genotype likelihoods (`--genotype-qualities`) and diploid ploidy;
> the VCF was assembled per chromosome with bcftools (Danecek et al., 2021) and hard-filtered at
> QUAL >= 20, GQ >= 20, DP 4-1000, missing rate <= 0.2, MAF >= 0.05, retaining bi-allelic SNPs only.
> Population structure was characterised by principal component analysis (plink2; Chang et al., 2015),
> model-based clustering (ADMIXTURE; Alexander et al., 2009) and a neighbour-joining tree built from
> IBS distances (ape; Paradis & Schliep, 2018). Population differentiation (Fst), nucleotide diversity
> (pi) and Tajima's D were computed with vcftools (Danecek et al., 2011); linkage disequilibrium decay
> was estimated with PopLDdecay (Zhang et al., 2019) and fitted on log-binned, pooled r2 values.
> Figures were produced with ggtree (Yu et al., 2017) and ggplot2 (Wickham, 2016).

## Methods 模板句（中文）

> 原始 reads 经 fastp（Chen et al., 2018）质控后，用 minibwa（Li & Homer, 2026）的 `-x sr` 短读模式
> 比对至参考基因组。变异检测采用 freebayes（Garrison & Marth, 2012）联合 calling，按 100 kb 分块，
> 映射质量 >=20、碱基质量 >=20，使用基因型似然（`--genotype-qualities`）、二倍体；VCF 按染色体用
> bcftools（Danecek et al., 2021）拼装后硬过滤（QUAL >=20、GQ >=20、DP 4-1000、缺失率 <=0.2、
> MAF >=0.05、仅保留双等位 SNP）。群体结构用 plink2（Chang et al., 2015）做 PCA、ADMIXTURE
> （Alexander et al., 2009）做模型聚类、基于 IBS 距离的邻接树（ape；Paradis & Schliep, 2018）做树。
> Fst、π、Tajima's D 用 vcftools（Danecek et al., 2011）计算；LD 衰减用 PopLDdecay
> （Zhang et al., 2019）统计后对 log 分箱池化 r² 拟合。图件用 ggtree（Yu et al., 2017）与
> ggplot2（Wickham, 2016）绘制。
