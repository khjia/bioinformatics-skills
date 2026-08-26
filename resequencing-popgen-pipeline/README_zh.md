# 重测序群体遗传全流程

从原始 FASTQ → VCF → 群体结构五图 + Fst/π/TajD/LD decay → HTML 报告。

## 环境要求

- **fb_calling** conda 环境：freebayes, bcftools, vcftools
- **pop_genetics** conda 环境：plink2, admixture, fasttree, angsd
- **R** 包：ggplot2, ggtree, ggtreeExtra, ggnewscale, ape, patchwork
- **minibwa**：快速 BWA-MEM 替代工具（~/bin/minibwa）

## 目录结构

```
项目/req/
├── 01.rawdata/          # 原始 FASTQ 文件
├── 02.cleandata/        # fastp 质控后的 reads
├── 03.Mapping/          # minibwa 比对 + samtools 排序
│   ├── ref/             # 参考基因组（genome.fasta 软链接）
│   ├── *.bam            # 排序后的 BAM 文件
│   └── summary/         # 比对统计
├── 04.Variants/         # freebayes 联合变异检测 + VCF 过滤
├── 05.population/       # 主分析（全部样本）
├── 06.pv25_only/        # 亚群分析（可选）
├── 07.selection/        # 选择分析（Fst, π, TajD, LD decay）
└── 00.report/           # HTML 报告 + 图表
```

## 分析流程

### 第1步：质控（fastp）
```bash
fastp -i R1.fq.gz -I R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
  --thread 8 --html sample_fastp.html --json sample_fastp.json
```

### 第2步：比对（minibwa + samtools）
```bash
minibwa map -t 8 -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
  ref/genome.fasta clean_R1.fq.gz clean_R2.fq.gz | \
  samtools sort -@ 4 -o sample.bam -
samtools index -c sample.bam   # 大基因组（>500Mb）用 CSI 索引
```

**注意**：染色体长度超过 500Mb 的基因组必须使用 CSI 索引（非 BAI）。BAI 上限为 2^29 = 536Mb。

### 第3步：变异检测（freebayes 联合 calling）
```bash
# 生成 BAM 列表
ls *.bam > bamlist.txt

# 联合 calling
freebayes -f ref/genome.fasta -L bamlist.txt \
  -p 1 -P 0 -E 3 --min-coverage 4 \
  --use-best-n-alleles 5 \
  > raw_variants.vcf.gz
```

### 第4步：VCF 过滤

| 参数 | 阈值 | 说明 |
|------|------|------|
| QUAL | ≥20 | 碱基质量过滤 |
| GQ | ≥20 | 基因型质量 |
| DP | [4, 1000] | 去除低/高覆盖度伪影 |
| F_MISSING | ≤0.2 | 最大 20% 缺失基因型 |
| MAF | ≥0.05 | 最小等位基因频率 |
| Alleles | 仅双等位 SNP | 去除 indel/多等位位点 |

```bash
bcftools view -i 'QUAL>=20' raw_variants.vcf.gz | \
bcftools view -i 'F_MISSING<=0.2' | \
bcftools view -m2 -M2 -v snps | \
bcftools view -i 'FORMAT/DP>=4 && FORMAT/DP<=1000 && FORMAT/GQ>=20' | \
bcftools view -i 'MAF>=0.05' > filtered.vcf
```

### 第5步：群体结构分析

**PCA（plink2）：**
```bash
plink2 --vcf filtered.vcf --set-all-var-ids @:# --rm-dup \
  --pca 10 --out pca_results
```

**ADMIXTURE：**
```bash
# 先转换为 bed/bim/fam 格式
plink2 --vcf filtered.vcf --set-all-var-ids @:# --rm-dup --make-bed --out pop

# 运行 K=2 到 K=6
for K in 2 3 4 5 6; do
  admixture --cv pop.bed $K -j8 | tee log${K}.out
done
```

**邻接树（FastTree2）：**
```bash
# 从 VCF 生成 SNP 比对
vcf2phylip --input filtered.vcf --output aligned.phy
FastTree2 -gtr -nt aligned.phy > tree.nwk
```

### 第6步：选择分析

**Fst（vcftools，滑窗法）：**
```bash
vcftools --vcf filtered.vcf --weir-fst-pop popA.txt --weir-fst-pop popB.txt \
  --fst-window-size 50000 --fst-window-step 10000 --out fst_results
```

**核苷酸多样性 π（vcftools）：**
```bash
vcftools --vcf filtered.vcf --keep popA.txt \
  --window-pi 50000 --window-pi-step 10000 --out pi_popA
```

**Tajima's D（vcftools）：**
```bash
vcftools --vcf filtered.vcf --keep popA.txt \
  --TajimaD 50000 --out tajd_popA
```

**LD 衰减（PopLDdecay）：**
```bash
PopLDdecay -InVCF filtered.vcf -OutStat ld_stat \
  -SubPop popA.txt -MAF 0.05 -Miss 0.2
perl Plot_OnePop.pl -inFile ld_stat.stat.gz -output ld_decay_plot
```

## 可视化（五图模板）

**图1：** 环形系统发育树（整枝着色，邻近向上传播）
**图2：** 环形树 + K=2 Admixture 外圈环（Q 矩阵需对齐树尖序）
**图3：** PCA（PC1/2 + PC1/3，群体着色）
**图4：** Admixture K=2-4（组内按祖先比例排序）
**图CV：** 交叉验证误差曲线

## 关键避坑

1. **plink2 LD**：必须用 `--set-all-var-ids @:# --rm-dup`，否则报重复 ID 错误
2. **ADMIXTURE**：染色体 ID 必须是整数；.Q 文件必须在当前工作目录
3. **ANGSD `-fold`**：仅在 `realSFS` 中使用
4. **ANGSD 2D SFS**：>140 样本不可行，改用 vcftools `--weir-fst-pop`（分钟级 vs 卡死）
5. **PCA 特征值**：除以所有特征值之和才是方差贡献率
6. **Q 矩阵对齐**：`Q[match(tree$tip.label, fam_order), ]`
7. **枝干着色**：邻近向上传播（子→父），非多数投票
8. **长任务**：用 `nohup cmd > log 2>&1 &` 避免包装器超时
9. **PopLDdecay**：完整路径 `/media/nfs1/hermes/miniforge3/bin/PopLDdecay`；用 `comment.char=""` 读 `.stat.gz`；避免 `-SubPop`（segfault）
10. **bcftools/PopLDdecay**：不在 pop_genetics 环境中，需 `conda activate fb_calling` 或用全路径

## 引用

如果使用本流程，请引用相关工具的文献。详见 SKILL.md 的 Citations 章节。

### 方法模板句

> "原始测序 reads 经 fastp（Chen et al. 2018）质控后，使用 minibwa（Li 2024）比对至参考基因组。变异检测采用 freebayes（Garrison & Marth 2012）联合 calling 模式，并使用 bcftools（Danecek et al. 2021）进行硬过滤（QUAL ≥ 20，基因型深度 DP 4–1000，最小等位基因频率 MAF ≥ 0.05）。群体结构通过 PCA（plink2；Chang et al. 2015）、模型聚类（ADMIXTURE；Alexander et al. 2009）和邻接树（FastTree2；Price et al. 2010）评估。群体分化（Fst）、核苷酸多样性（π）和 Tajima's D 使用 vcftools（Danecek et al. 2011）计算，连锁不平衡衰减使用 PopLDdecay（Zhang et al. 2019）估算。可视化使用 ggtree（Yu et al. 2017）和 ggplot2。"
