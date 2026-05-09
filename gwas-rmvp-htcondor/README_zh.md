# gwas-rmvp-htcondor — 中文使用说明

在 HTCondor 集群上跑多性状 × 多模型 GWAS（GLM + MLM + FarmCPU + BLINK）。
英文完整文档见 `SKILL.md`，本文件只覆盖上手要点和常踩坑。

## 何时用 / 何时不用

适用：
- 性状多（≥ 5）、想一次跑完 4 个模型并做共识筛选
- 集群有 HTCondor（CNS1 = Schedd）或 SLURM（少量改动）
- 已有 PLINK BED + 表型 TSV + PC 协变量 TSV

不适用：
- 数据小到笔记本能跑 → 直接 Rscript 跑 rMVP
- 只要一个模型且追求极致速度 → 用 GEMMA
- 二分类 + 罕见变异 → 用 REGENIE

## 输入要求

三份文件，IID 列必须完全一致（包括前导零）：

```
genotype_full_bed/full_snps_NNN.{bed,bim,fam}   # GWAS 用全量 SNP，不要 LD 过滤
inputs/pheno_filtered.tsv                       # FID IID trait_01 trait_02 ...
regenie_covariates/quant_pc3.tsv                # FID IID PC1 PC2 PC3
```

一致性校验（必跑一遍，对不上就改源头补零或去零）：

```bash
awk '{print $2}' genotype_full_bed/full_snps_*.fam | sort > /tmp/fam.ids
tail -n +2 inputs/pheno_filtered.tsv     | awk '{print $2}' | sort > /tmp/pheno.ids
tail -n +2 regenie_covariates/quant_pc3.tsv | awk '{print $2}' | sort > /tmp/pc.ids
comm -3 /tmp/fam.ids /tmp/pheno.ids   # 应该为空
comm -3 /tmp/fam.ids /tmp/pc.ids      # 应该为空
```

## 目录约定

```
<project>/05.gwas/
├── inputs/                         # 源数据
├── scripts/                        # 一次性准备脚本
├── genotype_full_bed/              # 全量 BED
├── regenie_covariates/             # PC 协变量
├── rMVP/                           # GWAS 主目录
│   ├── 01_prepare_data.R           # BED → big.matrix（只跑一次）
│   ├── 02_run_rMVP.R               # 主驱动（带重试 + K 指纹缓存）
│   ├── run_batch.sh                # HTCondor 任务包装
│   ├── rmvp.condor                 # 提交文件
│   ├── kinship.rds + kinship.fp    # K 矩阵 + 基因型指纹（自动失效）
│   ├── tune/                       # 可选：自适应 PC 扫描
│   └── final_results/              # 最终 CSV（硬链接到 tune/runs/）
├── post_gwas/                      # 5 步 Python 后处理
└── report/                         # 阶段报告
```

注意：报告写在阶段目录内 `05.gwas/report/`，不要写到全局 `00.report/`。

## 五步走

### 步骤 1：BED → rMVP big.matrix（只跑一次）

```bash
cd <proj>/05.gwas/rMVP
/media/nfs1/hermes/miniforge3/bin/Rscript 01_prepare_data.R
```

产出 `mvp.geno.{bin,desc,ind,map}`。15M SNP × 140 样本约 5–20 分钟。

### 步骤 2：主驱动 02_run_rMVP.R

直接用 `templates/02_run_rMVP.R`，已经内置三大特性，别再用裸版本：

A. **K 矩阵指纹缓存校验**：把 `mvp.geno.bin` 的 `<size>_<mtime>_<md5>` 写入 `kinship.fp`。
   下次跑如果指纹一致就读 `kinship.rds`，否则自动重算。基因型一变 K 自动失效，不会用错旧矩阵。

B. **性状级重试（在 R 里，不在 shell 里）**：每个性状外面套 `tryCatch` + 重试循环（默认 3 次，
   失败后 sleep 30s + `gc()`）。**为什么不在 shell 层 retry**：shell 层 `&& break` 会让一个性状
   失败就把整批重跑，浪费几小时。性状级重试能保住已成功的工作。

C. **四模型一次出**：`method = c("GLM","MLM","FarmCPU","BLINK")`。BLINK 是 rMVP 自带的，
   多花 10–20% 时间，但给共识筛选多一票独立证据。

环境变量（`run_batch.sh` 里设置）：

| 变量 | 默认值 | 作用 |
|---|---|---|
| `OMP_NUM_THREADS` | 16 | 传给 MVP() 的 ncpus |
| `RMVP_NPC` | 3 | PC 协变量数量（0 = 不用 PC）；自适应调优会覆盖 |
| `RMVP_RETRY` | 3 | 每个性状最多重试次数 |
| `RMVP_OUTDIR` | results | 输出目录；调优模式会改成 `results/trait_XX_pcNN` |

### 步骤 3：HTCondor 提交

`rmvp.condor` 关键字段：

```
universe       = vanilla
executable     = run_batch.sh
getenv         = True
should_transfer_files = NO       # 直接走 NFS，省得拷数据
request_cpus   = 16
request_memory = 64GB            # >200 样本或并行 ≥3 批改 96GB
arguments      = 1 5             # 跑性状 1–5
queue
arguments      = 6 10
queue
# ... 按总性状数分批
```

提交：

```bash
cd <proj>/05.gwas/rMVP
mkdir -p logs results
condor_submit rmvp.condor
condor_q -submitter $USER
```

### 步骤 4：监控

```bash
condor_q
condor_tail <cluster>.<process>      # 看实时 stderr
ls results/ | wc -l                  # 应该涨到 性状数 × 4
```

吞吐参考：16 核单性状单模型约 2–5 分钟。29 性状 × 4 模型分 6 批 ≈ 2–4 小时。

### 步骤 5：Python 后处理（5 个脚本，按顺序）

```bash
cd <proj>/05.gwas/post_gwas
mkdir -p plots summary
PY=/media/nfs1/hermes/miniforge3/bin/python
BASE=/abs/path/to/05.gwas

# 1. 提取 Bonferroni / suggestive 显著位点
$PY extract_all_bonf.py --base "$BASE"

# 2. 多模型共识 SNP + 位点聚类
$PY extract_high_confidence.py --base "$BASE" --min-models 2 --window 100000 --cutoff bonferroni

# 3. 单性状 Manhattan + QQ（4 模型 × N 性状 × 2 图）
$PY plot_rmvp_all.py

# 4. 多性状汇总图：堆叠 Manhattan + 性状×染色体热点热图
$PY plot_multitrait_summary.py --base "$BASE" --model FarmCPU --threshold suggestive

# 5. 单文件 HTML 报告（PNG 全部 base64 内嵌，可直接邮件附件发出去）
$PY generate_html_report.py --base "$BASE"
```

依赖：numpy / pandas / matplotlib / scipy（miniforge3 base 都有）。
HTML 报告纯 stdlib + pandas，没用 Jinja/Rmarkdown，浏览器直接打开就能看。

## 自适应 PC 调优（λ_GC 不正常时再用）

初版用 PC=3 跑完后，看 `per_trait_model_summary.tsv` 的 `lambda_gc` 列。
如果很多 (性状, 模型) 组合的 λ 落在 [0.85, 1.15] 之外，就上调优：

- PC 扫描点：`{0, 1, 2, 3, 5, 7, 10}`（跳过 4/6/8/9，收益递减）
- 选择规则：每个 (性状, 模型) 选 λ 最接近 1.0 且在 [0.85, 1.15] 之内的 PC 数
- 用**硬链接**把选中的 CSV 链到 `final_results/`，省几百 GB 中间文件

硬链接验证（必须做，删之前确认）：

```bash
# 所有 final CSV 必须 link count ≥ 2
stat -c '%h %n' final_results/*.csv | awk '$1<2'   # 输出为空才安全
rm -rf tune/runs/                                   # 验证后才能删
```

`tune/state/` 里的 JSON 是审计日志，体积小，留着。

## 多模型共识（high-confidence 位点）

`extract_high_confidence.py` 做的事：

1. 每个 (性状, 模型) 取通过 `--cutoff` 阈值的 SNP（默认 Bonferroni）
2. 一个 SNP 被几个模型同时打中？计数 ≥ `--min-models`（默认 2）的留下
3. 同染色体 `--window` bp（默认 100kb）内的 SNP 聚成一个位点

输出（在 `summary/` 下）：

- `high_confidence_snps.tsv` — 每行一个 (性状, SNP)，记录 `models_hit`、`n_models`、最小 p
- `high_confidence_loci.tsv` — 每行一个位点（trait, chrom, start, end, lead_snp, lead_p, models_union, n_snps）
- `model_agreement_matrix.tsv` — 模型两两重叠数（GLM∩MLM 极低而 GLM∩FarmCPU 高 → GLM 膨胀）

调参建议：

- 严格发文章：`--min-models 3 --cutoff bonferroni`
- 发现性 / QTL 优先：`--min-models 2 --cutoff suggestive`
- `--window` 跟你的 LD 衰减尺度匹配；自交植物 100kb 够用，异交可放到 500kb

## 三个常见坑

### 坑 1：IID 前导零被吃掉

症状：`MVP()` 报 "sample names do not match" 或输出全 NA。
原因：fam 用 3 位补零（001, 002...），rMVP 默默把前导零去掉了。
修复：驱动里已写 `ind <- sprintf("%03d", as.integer(readLines("mvp.geno.ind")))`。
你的项目如果用别的 ID 风格，改这一行。

### 坑 2：CSV 的 p 值列名不是 "p-value"

症状：`header.index("p-value")` 报 ValueError。
原因：rMVP 写出的 p 值列名是 `trait_XX.MODEL`（最后一列）。
修复：按位置取 `df.columns[-1]`，别按名字。

### 坑 3：matplotlib 渲染中文性状名变方块

症状：图标题里中文显示成 □□□。
**首选方案**：英文标题 + TSV 留中文映射，不要折腾字体。

```python
title = f'{trait_id} - {model}'   # 例如 "trait_01 - MLM"
```

中文名留在 `pheno_category_map.json` 和 `per_trait_model_summary.tsv` 的 `trait_zh` 列里。
零字体依赖，图文件可移植。

## HTCondor hold（内存超限）的处理

```bash
condor_qedit <JobId> RequestMemory 80000   # MB
condor_release <JobId>
```

经验值：150 样本 × 16M SNP 峰值 30–50GB；>300 样本或单节点并行 ≥3 批，一开始就开 96GB。

## 验证清单

跑完后按这个列表挨个确认：

- [ ] `ls final_results/*.csv | wc -l` == 性状数 × 4（GLM/MLM/FarmCPU/BLINK，跳过的不算）
- [ ] 用了自适应调优时：`stat -c '%h %n' final_results/*.csv | awk '$1<2'` 输出为空
- [ ] `per_trait_model_summary.tsv` 每个 (性状, 模型) 都有一行，包括 BLINK
- [ ] FarmCPU / MLM / BLINK 的 λ_GC 中位数接近 1.0
- [ ] `kinship.fp` 和 `kinship.rds` 同时存在
- [ ] plots/ 下 PNG 数 = 2 × 性状数 × 4，标题全英文
- [ ] `summary/high_confidence_loci.tsv` 在已知有 QTL 的性状上非空
- [ ] `report.html` 浏览器直接打开能看（PNG 已 base64 内嵌）
- [ ] HTCondor 日志没有遗留的 `[FAIL]` 行（重试用尽）
- [ ] 报告写在 `<proj>/05.gwas/report/`，不在 `00.report/`
- [ ] 自适应调优用过的话 `tune/runs/` 已删
- [ ] post_gwas/ 下没残留 `__pycache__/`

## 引用要求（发文章必看）

发文章必须引：

- **rMVP 工具本身** — Yin et al., 2021, GPB（doi:10.1016/j.gpb.2020.10.007）
- **报告的每个模型的原文**：
  - MLM → Yu et al., 2006, Nat Genet
  - FarmCPU → Liu et al., 2016, PLoS Genet
  - BLINK → Huang et al., 2019, GigaScience
  - GLM 没单一规范文献，引 rMVP 的实现即可
- **K 矩阵方法** — VanRaden, 2008, J Dairy Sci

完整文献列表 + BibTeX + Methods 模板句见 `SKILL.md` 的 **Citations** 章节。

不引会被审稿人逮到，植物/动物基因组方向是高频拒稿点。

## 相关 skill

- `gwas-gemma-slurm` — 单模型 GWAS 备选
- `population-genetics-vcf-analysis` — GWAS 上游 PCA / 结构分析
- `htcondor-cross-node-bioinformatics` — 集群故障排查
