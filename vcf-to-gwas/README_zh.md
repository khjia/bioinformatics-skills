# vcf-to-gwas — 从 VCF 到 GWAS 与单性状报告

一条命令链：VCF → plink2 质控 / LD 剪枝 / PCA → Fast3VmrMLM 多位点 GWAS →
最近基因注释 → 每个性状一份 7 图 12 表 HTML 报告。

## 目录

```
vcf-to-gwas/
├── SKILL.md                  # 主文档（流程、参数、方法、Citations）
├── README_zh.md              # 本文件
├── templates/
│   ├── 01_vcf2plink.sh       # VCF→bed、统一 ID、LD 剪枝、PCA
│   ├── 02_make_inputs.R      # trait.csv / PCA.csv（按 .fam 顺序）
│   ├── 03_run_fast3vmrmlm.R  # Fast3VmrMLM v2.0
│   ├── 04_annotate_qtl.py    # QTN 最近基因注释（±50 kb）
│   ├── 05_gene_coords.R      # GFF3 → gene_coords.rds
│   ├── 20_trait_report_figs.R # 单性状图表
│   ├── 21_trait_report.Rmd   # 报告模板
│   ├── 22_render_all.R       # 批量渲染
│   ├── 00_config.example.R   # 报告配置模板
│   ├── pipeline.env.example  # 流程环境变量模板
│   └── run_all.sh            # 一键串联（可加 RUN_REPORTS=1）
└── references/pitfalls.md    # 实测坑点
```

## 快速开始

```bash
cp templates/pipeline.env.example pipeline.env   # 改 PROJECT、工具路径、阈值
set -a; source pipeline.env; set +a

bash templates/run_all.sh                        # Phase A：VCF → GWAS → 注释
RUN_REPORTS=1 bash templates/run_all.sh          # Phase A + B（需先建 00_config.R）
```

单步运行：

```bash
bash    templates/01_vcf2plink.sh
Rscript templates/02_make_inputs.R
Rscript templates/03_run_fast3vmrmlm.R
python3 templates/04_annotate_qtl.py
Rscript templates/05_gene_coords.R
```

## 关键约定

- `trait.csv` 首列表头必须是 `<Phenotype>`，`PCA.csv` 首列表头必须是 `<pca>`；
  行顺序必须与基因型 `.fam` 完全一致。
- 染色体名非标准时（如 `chr01a`）必须加 `--allow-extra-chr`。
- 变异 ID 统一为 `chr:pos:ref:alt`，后续所有连接都依赖它。
- 显著性：Bonferroni p = 0.05 / SNP 数，2,494,579 个 SNP 对应 **LOD = 7.7**。
- 区域信号取自 Fast3VmrMLM 的 QTN 输出（按 50 kb 去冗余），不要用单标记条件模型
  （多数性状功效不足，会得到 0 信号）。

## 输出

- `06.Fast3VmrMLM/results/<trait>_result.xlsx`：QTN 表（LOD、add/dom、r2、P、显著性）
- `06.Fast3VmrMLM/results/<trait>_midresult.csv`：单标记扫描
- `06.Fast3VmrMLM/report/qtl_annotated.tsv`：QTN-基因注释
- `09.postgwas2/report/*.html`：每性状一份报告（7 图 12 表）

## 常见问题

见 `references/pitfalls.md`：BOM/编码、Windows 文件名乱码、`sprintf` 的 `%`、
表格单元格 `+`/`-` 被当成列表、区域缺失率偏高时的 `--geno` 回退、LOD 与
−log₁₀P 不可混用等。
