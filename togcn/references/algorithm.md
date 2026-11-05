# 官方 TO-GCN 算法（源码级整理）

源码位置（示例服务器）：`/media/40T/public/software/TO-GCN-master`
包含：`Cutoff.cpp`、`GCN.cpp`、`TO-GCN.cpp`（comparative，双条件）
以及 `Single_Time-series_data/source_code/` 下的 `*_single.cpp`（单条件）。

## 1. Cutoff —— 定阈值

```
Cutoff #Cond1_samples #Cond2_samples TF.tsv gene.tsv
```

- 对每个 TF×gene 对，在条件 1、条件 2 下各算一个 Pearson r
- 把 r 打进 201 个 bin（-1.00 ~ 1.00，步长 0.01），算 PDF / CDF
- 正阈值：从高到低扫，第一个 `CDF_desc > 0.05` 的 bin，再 +2 bin
- 负阈值：从低到高扫，第一个 `CDF_asc > 0.05` 的 bin
- 输出 `PCC_histogram.tsv` + 推荐的 `GCN` 命令行

**注意**：这个 5% 尾是**经验分布**，不是理论 null 分布。真实共表达数据里
r 的分布远比 N(0,1/√n) 宽，所以推荐值（0.85~0.93）远高于理论 5%（在 n=15 时约 0.44）。

## 2. GCN —— 建八类共表达网络

```
GCN #Cond1 #Cond2 TF.tsv gene.tsv posC1 posC2 negC1 negC2
```

输出 8 个 csv，列固定为 `TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2`：

| 类型 | 条件 | 用途 |
|---|---|---|
| C1+C2+ | r1 ≥ p1 且 r2 ≥ p2 | **TO-GCN 分层默认用这一张** |
| C1+C20 | r1 ≥ p1 且 r2 ∈ (-0.5, 0.5) | 条件间不一致 |
| C10C2+ | r2 ≥ p2 且 r1 ∈ (-0.5, 0.5) | 同上 |
| C1-C2- | r1 < n1 且 r2 < n2 | 抑制关系，**不参与分层** |
| C1-C20 / C10C2- | 单侧负 | |
| C1+C2- / C1-C2+ | 方向相反 | |

其中 ±0.5 是源码里的硬编码 `pos_no_cutoff / neg_no_cutoff`。

## 3. TO-GCN —— BFS 分层

```
TO-GCN #Cond1 #Cond2 TF.tsv gene.tsv posC1 posC2 seeds.txt coex_type
# coex_type: 0 = C1+C2+, 1 = C1+C20, 2 = C10C2+
```

```c
set_neighbor_pos(check_ID, L):
  for 每条未消费的边 e:
     if e.query == check_ID: e.checked = 1; if (e.target 是 TF 且 level == -1) level = L
  for 每条未消费的边 e:
     if e.target == check_ID: e.checked = 2; if (e.query 是 TF 且 level == -1) level = L
```

- 种子 level = 1，输出时 `level - 1`（0-based）
- **只有 TF 表里的节点会被赋层**：`TF.tsv` 传 TF 列表 → 只有 TF 有层级；
  传全部基因 → 全部基因有层级
- 边是“查询→目标”的有向记录，但 BFS 两个方向都走，等价于**无向图最短路**
- `checked` 只保证每条边被消费一次，不改变最短路语义

等价实现（`templates/togcn_fast.cpp`）就是：建无向邻接表 → BFS 记录距离。

## 4. single 版差异

```
Cutoff_single #features matrix.tsv
TO-GCN_single #features matrix.tsv seeds.txt pos_cutoff
```

- 只有一个条件，只接受**正阈值**
- 同一个文件同时作为 TF 表和基因表
- 官方 `function_three()` 有个 bug：`level==1` 原样输出、其余 `level-1`，
  导致种子和一级邻居混在一起；**改成统一 `level-1`** 才与 comparative 版一致

## 5. 层级的含义

层号不是“调控代数”，而是**图距离**：网络越稠，直径越小，层数越少；
超过某一阈值网络断裂，层数又会骤降。报告时必须同时给出边数与覆盖节点数。
