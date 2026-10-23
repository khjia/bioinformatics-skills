# bioinformatics-skills

Bioinformatics analysis skills collection.

每个 skill 一个目录，内含 `SKILL.md`（主文档）、`README_zh.md`（中文使用说明）、可选
`templates/`（可直接复用的脚本模板）、`references/`、`scripts/`。所有 skill 的 SKILL.md
末尾都有 **Citations** 章节（软件与方法原文 + BibTeX + Methods 模板句），发表论文时直接引用。

## Skills

### `rnaseq-hisat2-featurecounts-deseq2` — 植物批量 RNA-seq 主线流程

raw FASTQ → fastp 质控 → HISAT2 比对 → featureCounts 定量 → TPM + 组织级低表达过滤
→ DESeq2 差异表达 → clusterProfiler GO/KEGG 多组富集。

```
rnaseq-hisat2-featurecounts-deseq2/
├── SKILL.md              # 主文档（中文，含全部实测坑点 + Citations）
├── README_zh.md          # 中文速查说明
└── templates/
    ├── run_fastp.sh                    # PE 批量质控（幂等）
    ├── build_hisat2_index.sh           # GFF3→GTF → splice_sites/exons → 建索引（含 python PATH 修复）
    ├── run_hisat2.sh                   # 批量比对（sam_nfs / lowio 两种模式，含 0-read guard）
    ├── calc_tpm_filter.py              # TPM 计算 + 组织级过滤，输出 4 个矩阵
    ├── run_deseq2.R                    # DESeq2 多比较 + volcano/MA/PCA/heatmap
    └── run_multigroup_enrichment.R     # 多组 clusterProfiler 富集（x=group 散点图）
```

核心坑点（详见 SKILL.md）：
- `exit 134` 有两个完全不同的原因——**NFS 满** 或 **/tmp 满**（大基因组每样本 SAM 10-21 GB，
  并行 N 个瞬间吃满 /tmp）。必须先 `df -h /tmp` 再决定 SAM 写哪
- `hisat2 -x` 传 **index 前缀**，不是 fasta 路径（报错信息极具误导性）
- 样本名带括号 `H4(3)` → **静默产出 0-read BAM 且 exit 0**，必须加 post-mapping guard
- `--no-unal` 让 flagstat 报 100% mapped，真实比对率看 `.hisat2.log`
- GFF3 必须先转 GTF，否则 splice_sites/exons 是 0 字节
- `enricher()` 的 TERM2GENE 列名必须大写 `TERM`/`GENE`，否则静默返回 0 条结果
- TPM 低表达过滤是**组织级**（任一组织平均 TPM ≥ 1），全局平均会删掉组织特异基因

### `resequencing-popgen-pipeline` — 植物重测序群体遗传主线流程

raw FASTQ → fastp 质控 → minibwa 比对 → **freebayes 联合 calling**（100 kb 分块、SGE 数组、
可断点续跑）→ VCF concat/norm/硬过滤 → plink2 LD+PCA / KING 亲缘 / ADMIXTURE
→ IBS 邻接树 → π / Tajima's D / Fst / LD decay / ROH → 群体结构五图 + HTML 报告。

```
resequencing-popgen-pipeline/
├── SKILL.md              # 主文档（中文，26 条避坑 + Citations）
├── README_zh.md          # 中文速查说明
└── templates/
    ├── make_regions.py               # .fai → 100 kb 分块区间
    ├── make_work_units.py            # 分块 → N 个工作单元（轮转均分）
    ├── freebayes_call_unit.sh        # 跑一个工作单元（可断点续跑）
    ├── freebayes_call_unit_array.sub # SGE 数组提交模板
    ├── circular_tree_admixture_rings.R   # 环形树 + K=2..N admixture 外圈
    └── tree_structure_combined.R         # 矩形 NJ 树 + admixture 条形图合一
```

**硬性约定**：SNP calling 一律使用 **freebayes**，不使用 `bcftools mpileup/call`；bcftools 只用于
calling 之后的 VCF 工具链（concat / norm / index / stats / +fixploidy）。
