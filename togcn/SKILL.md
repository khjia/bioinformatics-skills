---
name: togcn
description: >-
  Build a time-ordered gene co-expression network (TO-GCN) from time-series
  RNA-seq and turn it into publication figures: iTAK transcription-factor
  identification, the official Cutoff/GCN/TO-GCN pipeline (or the verified
  equivalent implementation), BFS layer assignment, size-and-trend based layer
  merging, layer expression heatmaps, concentric level-by-level regulatory
  circles, and seed-gene expression panels.
  Use for TO-GCN、时序共表达网络、分层调控网络、转录因子层级图、TF 调控圈图、
  time-ordered co-expression、hierarchical TF network.
license: CC-BY-4.0
---

# TO-GCN 时序共表达网络：构建与出图

面向单一时序或多条件时序 RNA-seq，复现 Nature Communications 一类论文中的
TO-GCN（TO-GCN: time-ordered gene co-expression network，BFS 分层）。

> **上游衔接**：本 skill 是 `rnaseq-hisat2-featurecounts-deseq2` 流水线的 **Stage 7**。
> fastp → HISAT2 → featureCounts → TPM → DESeq2 完成后，把「差异基因 ∪ 全部 TF」的
> TPM 矩阵交给这里建网。报告里必须同时给出上游 QC（原始 reads、比对率、
> featureCounts 分配率、表达过滤前后基因数），否则网络规模无法解释。


## 何时用 / 何时不用

用：把时序转录组做成“层级 + 调控网络 + 层级热图/圈图”；需要 TF 层级级联。
不用：只做差异表达或富集；数据没有 ≥4 个时间点（BFS 分层没有意义）。

## 0. 先确认三件事

1. **几个条件**：两个条件 → 用 comparative 版（C1/C2 双条件）；单一时间序列 →
   用 `Single_Time-series_data` 版。判定错了整条链都要重做。
2. **几个时间点 × 几个重复**：n 太小（如 5 个点、每组 3 重复 = 15 样本）时，
   官方 5% 尾阈值会偏松，需要看“小世界”诊断（见第 5 步）。
3. **TF 列表从哪来**：必须用 iTAK（`tf_classification.txt`），不要用 PFAM 代理。

## 1. 输入准备

```r
# 计数矩阵 → TPM（基因长度取 GTF 外显子并集）
counts_all <- read.delim("raw_counts_matrix.tsv", row.names = 1, check.names = FALSE)
sel  <- grep("^(0|30|60|90|120)_[123]$", colnames(counts_all), value = TRUE)
sel  <- sel[order(as.numeric(sub("_.*$","",sel)), sub("^.*_","",sel))]   # 顺序必须对
counts <- as.matrix(counts_all[, sel]); rm(counts_all)
len  <- read.delim("gene_length.tsv", header = FALSE)     # 每个基因外显子并集长度
L    <- setNames(len$V2, len$V1)
counts <- counts[rownames(counts) %in% names(L), ]
rpk  <- counts / (L[rownames(counts)] / 1000)
tpm  <- sweep(rpk, 2L, colSums(rpk), "/") * 1e6

# 组水平过滤：任一组（时间点/组织/处理）平均 TPM > 1 即保留
grp  <- sub("_.*$", "", sel)                              # 组标签 = 时间点
mean_by_grp <- sapply(split(seq_along(sel), grp),
                      function(i) rowMeans(tpm[, i, drop = FALSE]))
keep <- apply(mean_by_grp, 1, max) > 1
tpm  <- tpm[keep, , drop = FALSE]
```

**过滤规则统一为「任一组平均 TPM > 1 则保留」**：

- 组 = 时间点（或组织/处理），先在组内对重复求平均，再取各组最大值
- 不要用全局平均、也不要用计数阈值（`rowMeans(counts) >= 5` 一类）——
  那样会把"只在某一个时间点高表达"的关键基因滤掉，而这类基因恰恰是时序网络的
  起点（种子）和早期层
- 过滤必须在算完 TPM 之后做，不能先过滤计数再算 TPM（会改变分母）

关键：**列顺序必须按"时间 × 重复"排好**，`times*3` 之类的写法是错的，必须
`np.repeat(times, 3)` / R 里显式排序。这个坑毁过一次种子选择。

写矩阵给 C++ 用：`基因ID \t 样本1 ... 样本n`，**无表头、无行名**（官方 `fscanf` 要求）。

可选的基因集：全表达基因；DEG（任一比较 vs 0 min，padj<0.05 且 |log2FC|>1）∪ 全部 TF。
DEG 集会显著缩小网络并改善层间均衡。

## 2. iTAK 定 TF

```bash
iTAK -p -i proteins.fa -o itak_out -m plants -b 2>&1 | tee itak.log
# 输出 tf_classification.txt：protein_id \t family \t TF|TR \t path
```

蛋白 ID 常带 isoform 后缀（`.1`），比对表达矩阵前要 `sub("\\..*$","")`。

## 3. 阈值：用官方 Cutoff，不要拍脑袋

```bash
g++ -O2 -o Cutoff_single Cutoff_single.cpp
g++ -O2 -o TO-GCN_single  TO-GCN_single.cpp
Cutoff_single <n_samples> matrix.tsv        # 打印推荐正/负阈值 + PCC_histogram.tsv
```

规则（源码）：PCC 分 201 个 0.01 bin；上尾 CDF>0.05 的 bin +2 即为正阈值；
下尾同理为负阈值。
**同一个基因集，阈值只取决于数据**；15 个样本、单一组织时典型值 0.86–0.93。

## 4. 建网 + BFS 分层

```bash
TO-GCN_single <n_samples> matrix.tsv seeds.txt <pos_cutoff>
# → Node_level.csv（TF gene ID, level in GCN）、Node_relation.csv
```

- 层级 = 从种子出发的**无向图最短路距离**（+1 后输出为 0-based）
- 官方 BFS 是 O(节点×边)，>5000 节点会跑几小时 → 换 `templates/togcn_fast.cpp`
  （算法等价，邻接表实现；**必须先用官方二进制在小数据集上校验输出逐行一致**）

种子选择（官方用 MFSelector，但只支持 R<4.0，通常装不上）：

```python
T0 == max(5 个时间点均值)          # 第一时期是峰值
T0 >= 50 CPM                       # 不是噪声
T0 / mean(T30..T120) >= 1.5        # 后面确实低
```

通常 10–30 个种子，全部落在 level 0。

## 5. 诊断：阈值是否把网络压成“小世界”

看层大小分布。如果 **80%+ 节点挤在前 5 层**，说明网络太稠、层级无意义，需要提高阈值。
扫一遍：0.88 / 0.92 / 0.94 / 0.95 / 0.96 / 0.97，记录
「边数 / 覆盖节点 / 最大层 / 层分布」，选层级能展开又能保持连通的最小阈值。
过高会突然碎裂（覆盖节点断崖式下降）。

## 6. 合并层

```bash
python3 merge_layers.py --run <run_dir> --tag <tag> --min_n 100 --max_n 2000
# → figs/groups_<tag>.csv（gene, group, group_order, layer）
```

规则：
1. `n < min_n` 的层并入**相邻两层中趋势更相似**（均值 z 曲线 Pearson r 更大）的一侧；
2. 任何合并结果不得超过 `max_n`（两端都超时选更小的）；
3. 原生层超过 `max_n` **不拆**（拆分只会切出趋势相同的子块，除非用户明确要求）。

不要用“趋势相似”单独作为合并理由：大层即使 r=0.99 也各自独立。

## 7. 出图

| 图 | 脚本 | 要点 |
|---|---|---|
| 层热图（z-score + log2 CPM，mean/median 双面板） | `templates/layer_heatmaps.R` | 基因级 z-score 后再取层均值，避免被高表达基因主导 |
| 分层调控圈（左→右逐层，TF 内圈/基因外圈） | `templates/network_circles.R` | 圈径 ∝ √n；TF–基因 + TF–TF；去掉基因–基因；方向过滤：TF(k)→基因(k+1) 与同层保留，反向删除 |
| 种子基因表达 | `templates/seed_expression.R` | z-score 热图 + log2(CPM+1) 折线（mean±SD） |

出图规范：图例必须解释颜色/线型映射；坐标轴带单位；同层用 `coord_fixed()`。

## 8. 必须写进 Methods 的参数

- 样本数、时间点、重复数；基因过滤阈值
- TF 来源与版本（iTAK 版本、TF 数）
- 阈值 + 它是怎么来的（Cutoff 的 5% 尾规则）
- 种子基因列表与选择规则
- 层合并参数（min_n / max_n）与合并结果（几层 → 几组）
- 边数、覆盖节点数、跨层边数

## 参考

- `references/algorithm.md`：官方三个程序的完整算法（含 C1/C2 语义、八类共表达）
- `references/pitfalls.md`：踩过的坑与数字自检清单
- `templates/`：可直接改用的脚本
