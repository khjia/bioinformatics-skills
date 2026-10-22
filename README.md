# bioinformatics-skills

Bioinformatics analysis skills collection.

每个 skill 一个目录，内含 `SKILL.md`（主文档）、可选 `README_zh.md`（中文使用说明）、可选
`templates/`（可直接复用的脚本模板）与 `references/`、`scripts/`。所有 skill 的 SKILL.md
末尾都有 **Citations** 章节（软件与方法原文 + BibTeX + Methods 模板句），发表论文时直接引用。

## Skills

- **gwas** — complete GWAS pipeline (Fast3VmrMLM + post-GWAS + enrichment + report)
- **rna-seq-pipeline** — bulk RNA-seq from ENA download to DEG interaction analysis
- **togcn-time-ordered-coexpression** — time-ordered gene coexpression network (TO-GCN) construction from developmental time-series RNA-seq, with seed selection, BFS level assignment, two-circle network visualization, and stage-wise GO/KEGG enrichment (中文说明见各 skill 的 README_zh.md)
- **resequencing-popgen-pipeline** — end-to-end plant resequencing population genetics: fastp QC → minibwa mapping → **freebayes joint calling** (100-kb windows, SGE array, resumable) → VCF concat/norm/hard-filter → plink2 LD+PCA / KING kinship / ADMIXTURE → IBS neighbour-joining tree → π / Tajima's D / Fst / LD decay / ROH → 5-figure population structure + HTML report. Includes runnable templates for freebayes chunking, SGE array submission, and the circular-tree/admixture-ring plots (中文说明见 `resequencing-popgen-pipeline/README_zh.md`)

### resequencing-popgen-pipeline 目录

```
resequencing-popgen-pipeline/
├── SKILL.md              # 主文档（中文，含 26 条避坑 + Citations）
├── README_zh.md          # 中文使用说明（速查版）
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
