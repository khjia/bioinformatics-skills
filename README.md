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

### `vcf-to-gwas` — VCF → GWAS → 单性状报告

VCF → plink2 质控 / LD 剪枝 / PCA → Fast3VmrMLM 多位点 GWAS → 最近基因注释
→ 每个性状一份 HTML 报告（Manhattan/QQ、主峰区域图、±1 Mb 全部基因、区域 QTN 信号、
候选基因单倍型）→ 可选 SuSiE 精细定位（含 LD 诊断）→ 可选 GWAS Catalog 已知位点查询。

```
vcf-to-gwas/
├── SKILL.md              # 主文档（流程、参数、方法、实测坑点 + Citations）
├── README_zh.md          # 中文速查说明
├── templates/
│   ├── 01_vcf2plink.sh              # VCF→bed、统一 ID、LD 剪枝、PCA
│   ├── 02_make_inputs.R             # trait.csv / PCA.csv（按 .fam 顺序）
│   ├── 03_run_fast3vmrmlm.R         # Fast3VmrMLM v2.0
│   ├── 04_annotate_qtl.py           # QTN 最近基因注释（±50 kb）
│   ├── 05_gene_coords.R             # GFF3 → gene_coords.rds
│   ├── 20_trait_report_figs.R       # 单性状图表
│   ├── 21_trait_report.Rmd          # 报告模板（含 SuSiE/变异位置补充章节）
│   ├── 22_render_all.R              # 批量渲染
│   ├── 23_finemap_susie.R           # SuSiE 精细定位（个体基因型、in-sample LD）
│   ├── 26_known_loci.py             # GWAS Catalog 已知位点查询
│   ├── 00_config.example.R          # 报告/精细定位配置模板
│   ├── pipeline.env.example         # 流程环境变量模板
│   └── run_all.sh                   # 一键串联（RUN_REPORTS=1 时含报告）
└── references/
    ├── fine-mapping.md              # SuSiE/FINEMAP 选择、LD 诊断、坑点
    ├── known-loci.md                # GWAS Catalog / 跨物种查询注意事项
    └── pitfalls.md                  # 全流程实测坑点
```

核心坑点（详见 SKILL.md）：
- Fast3VmrMLM 的 **LOD 与单标记 −log₁₀P 不可混用**（前者来自多位点模型）
- `trait.csv` / `PCA.csv` 按行对齐 `.fam`，表头必须是 `<Phenotype>` / `<pca>`
- 非标准染色体需 `--allow-extra-chr`；变异 ID 统一为 `chr:pos:ref:alt`
- 区域缺失率 10–17% 时，提取基因型要放宽 `--geno`（0.1 → 0.25 → 0.5）
- 表格单元格若只有 `+`/`-` 会被 pandoc 当成列表，破坏整张表

### `togcn` — 时序共表达网络（TO-GCN）与分层调控

单一时序或多条件时序 RNA-seq 的层级调控网络构建：

```
togcn/
├── SKILL.md              # 主文档（何时用、8 步流程、出图与 Methods 要点）
├── references/
│   ├── algorithm.md      # 官方三程序算法（C1/C2 语义、八类共表达、BFS 源码级整理）
│   └── pitfalls.md       # 19 条实测坑点 + 数字自检表
└── templates/
    ├── Cutoff_single.cpp / TO-GCN_single.cpp   # 官方源码
    ├── togcn_fast.cpp                          # 等价优化实现（邻接表）
    ├── build_matrices.R                        # TPM 建矩阵 + 组水平过滤
    ├── build_tree.py                           # 重建 BFS 树（带父节点）
    ├── merge_layers.py                         # 层合并（min_n / max_n）
    ├── layer_heatmaps.R                        # 层热图（基因级 z-score）
    ├── network_circles.R                       # 逐层调控圈
    └── seed_expression.R                       # 种子基因表达
```

核心要点（详见 SKILL.md）：

- 阈值必须由官方 `Cutoff` 给出，并扫描确认网络没被压成“小世界”（最大层占比 > 50% 就要警惕）
- 自写实现必须与官方二进制在同一矩阵上 `diff` 层级文件一致才能用
- 层级是**图距离**，报告时必须同时给出边数与覆盖节点数
- 跨层边必须保留：把基因按“最强 TF 所在层”放置会让层间调控消失

### `analysis-report` — 分析报告写作（材料方法 → 结果 → 解读）

把一次完整分析（原始数据到最终图表）写成自包含 HTML 报告：

```
analysis-report/
├── SKILL.md                    # 骨架 + 三条硬规则 + 技术实现
├── references/style-guide.md   # 逐段模板、措辞规则、自检清单、常见错误
└── templates/
    ├── report_template.Rmd     # 可直接改用的骨架
    └── render.R                # 渲染脚本（含 pandoc PATH 处理）
```

核心要点（详见 SKILL.md）：

- 每一步固定三段：**材料方法**（软件/版本/参数/阈值来源）→ **结果**（只放事实）→ **解读**（2–4 句）
- 解读回答三件事：数字说明什么、能推出什么、有什么保留；不许拿后文结论倒推
- 交付前跑自检清单：正文数字与表一致、每图有编号与引用、无服务器路径与账号
