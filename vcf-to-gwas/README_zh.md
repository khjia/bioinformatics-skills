# vcf-to-gwas — 从 VCF 到 GWAS、精细定位与单性状报告

VCF → plink2 质控 / LD 剪枝 / PCA → Fast3VmrMLM 多位点 GWAS → QTN 注释
→ 每个性状一份 HTML 报告（Manhattan/QQ、主峰区域、±1 Mb 全部基因、区域 QTN、单倍型）
→ 可选 SuSiE 精细定位 → 可选 GWAS Catalog 已知位点查询。

## 目录

```
vcf-to-gwas/
├── SKILL.md
├── README_zh.md
├── templates/
│   ├── 01_vcf2plink.sh         # VCF→bed、统一 ID、LD 剪枝、PCA
│   ├── 02_make_inputs.R        # trait.csv / PCA.csv（按 .fam 顺序）
│   ├── 03_run_fast3vmrmlm.R    # Fast3VmrMLM v2.0
│   ├── 04_annotate_qtl.py      # QTN 最近基因注释（±50 kb）
│   ├── 05_gene_coords.R        # GFF3 → gene_coords.rds
│   ├── 20_trait_report_figs.R  # 单性状图表
│   ├── 21_trait_report.Rmd     # 报告模板（含 SuSiE/变异位置补充章节）
│   ├── 22_render_all.R         # 批量渲染
│   ├── 23_finemap_susie.R      # SuSiE 精细定位（个体基因型、in-sample LD）
│   ├── 26_known_loci.py        # GWAS Catalog 已知位点查询
│   ├── 00_config.example.R     # 报告/精细定位配置模板
│   ├── pipeline.env.example    # 流程环境变量模板
│   └── run_all.sh              # 一键串联（RUN_REPORTS=1 时含报告）
└── references/
    ├── fine-mapping.md         # SuSiE/FINEMAP 选择、LD 诊断、坑点
    ├── known-loci.md           # GWAS Catalog / 跨物种查询注意事项
    └── pitfalls.md             # 全流程实测坑点
```

## 快速开始

```bash
cp templates/pipeline.env.example pipeline.env   # 改 PROJECT、工具、阈值
set -a; source pipeline.env; set +a

bash templates/run_all.sh                        # Phase A
RUN_REPORTS=1 bash templates/run_all.sh          # Phase A+B（需先建 00_config.R）

export GWAS_WDIR=/path/to/09.postgwas2
Rscript templates/23_finemap_susie.R fa fa_C24_1 # Phase C
python3 templates/26_known_loci.py --gene AT5G43760 --out known.tsv  # Phase D
```

## 关键约定

- `trait.csv` 首列表头 `<Phenotype>`，`PCA.csv` 首列表头 `<pca>`；行顺序与 `.fam` 一致。
- 非标准染色体加 `--allow-extra-chr`；变异 ID 统一 `chr:pos:ref:alt`。
- 显著性：Bonferroni p = 0.05 / SNP 数（2,494,579 个 SNP ↔ LOD 7.7）。
- 区域信号取 Fast3VmrMLM QTN 按 50 kb 去冗余；不要用单标记条件模型代替。
- 低 N（≈200）时 SuSiE 报告 90% 可信集，PIP 仅作排序；用 L=5/10/15 检查稳定性。

## 常见问题

见 `references/pitfalls.md`、`references/fine-mapping.md`、`references/known-loci.md`。
