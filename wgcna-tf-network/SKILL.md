---
name: wgcna-tf-network
description: WGCNA module TF-structural gene co-expression network visualization. iTAK TF annotation + per-pathway ring-layout network plots with NPG colors.
tags: [wgcna, transcription-factor, network, itak, visualization, co-expression]
---

# WGCNA TF-Structural Gene Co-expression Network

## When to Use
- WGCNA 模块分析后，需要可视化 TF 与结构基因的共表达调控网络
- 需要对模块基因进行 iTAK 转录因子注释
- 需要按 GO 通路拆分绘制 per-pathway 网络图

## Pipeline Overview

```
Module genes → Extract protein → hmmscan (iTAK DB) → TF classification
                                                          ↓
GO enrichment → Select pathways → Split TF/Structural → Ring-layout network plot
```

## Step 1: iTAK TF Annotation

### Prerequisites
- iTAK database: `Tfam_domain.hmm`, `TF_selfbuild.hmm`, `TF_Rule.txt`
- Protein sequences for module genes (FASTA)
- hmmscan (HMMER3)

### Run hmmscan
```bash
# Against Tfam_domain.hmm
hmmscan --domtblout ${OUT}/module_tfam.domtbl --noali -E 1e-5 \
  ${ITAK_DB}/Tfam_domain.hmm ${PROTEINS} > /dev/null

# Against TF_selfbuild.hmm
hmmscan --domtblout ${OUT}/module_selfbuild.domtbl --noali -E 1e-5 \
  ${ITAK_DB}/TF_selfbuild.hmm ${PROTEINS} > /dev/null
```

### Classification Logic
1. Build NAME→ACC mapping from HMM files (NAME line → ACC line, strip version)
2. Parse domtblout: map domain hits to Pfam ACC per gene
3. Parse TF_Rule.txt (multi-line record format, `//` separator):
   - `Required`: domain requirements (AND=`--`, OR=`:`, count=`#N`)
   - `Forbidden`: domains that exclude classification
   - `Name`: TF subfamily name
   - `Type`: TF or TR
4. First-match wins (rules are priority-ordered)

### Pitfalls
- TF_Rule.txt uses Pfam ACC (e.g., `PF00847#2`) for Tfam domains
- TF_selfbuild domains use their NAME directly (no Pfam ACC)
- Must build NAME→ACC mapping from HMM file headers
- Gene IDs in domtblout include `.1` suffix; strip for matching with WGCNA gene lists
- `G2-like` in selfbuild maps to itself (no Pfam ACC)

## Step 2: Network Visualization

### Design
- **Inner ring (red #E64B35)**: TF genes — pathway TFs + top-N module TFs by co-expression weight to structural genes
- **Outer ring (cyan #4DBBD5)**: Structural genes (non-TF pathway genes)
- **Edge width**: proportional to WGCNA co-expression weight
- **TF→Structural edges**: bold red (#E64B35, alpha=0.6, width×1.8)
- **Structural↔Structural edges**: thin cyan (#91D1C2, alpha=0.3, width×0.8)
- **Labels**: black, TF labels include family name in parentheses

### NPG Color Palette
```python
C_TF = '#E64B35'       # TF nodes + TF-struct edges
C_STRUCT = '#4DBBD5'   # Structural nodes
C_EDGE_SS = '#91D1C2'  # Struct-struct edges
# Additional NPG: '#3C5488', '#00A087', '#F39B7F', '#8491B4'
```

### Layout Algorithm
```python
# Inner ring (TF)
for i, node in enumerate(tfs):
    angle = 2 * pi * i / n_tf - pi/2
    pos[node] = (r_inner * cos(angle), r_inner * sin(angle))

# Outer ring (structural)
for i, node in enumerate(structs):
    angle = 2 * pi * i / n_struct - pi/2
    pos[node] = (r_outer * cos(angle), r_outer * sin(angle))
```

### Top-N TF Selection (when pathway has few TFs)
For each pathway, find module TFs with highest total co-expression weight to the pathway's structural genes:
```python
tf_weights = defaultdict(float)
for edge in module_edges:
    if edge.node1 in tf_set and edge.node2 in struct_set:
        tf_weights[edge.node1] += edge.weight
top_tfs = sorted(tf_weights, key=tf_weights.get, reverse=True)[:N]
```

### Label Positioning
- Radially outward from node center
- Horizontal alignment based on angle (left/right/center)
- TF labels: bold, include `\n(family_name)`

## Dependencies
- Python: pandas, numpy, networkx, matplotlib
- System: hmmscan (HMMER3)
- Database: iTAK-db-v1

## Output
- Per-pathway PNG + PDF (300 dpi)
- TF annotation TSV: gene, gene_id, TF_family, type, full_family

## Citations

- **iTAK**: Zheng Y, et al. (2016) iTAK: A Program for Genome-wide Prediction and Classification of Plant Transcription Factors, Transcriptional Regulators, and Protein Kinases. *Molecular Plant* 9(12):1667-1670. doi:10.1016/j.molp.2016.09.014
- **WGCNA**: Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559. doi:10.1186/1471-2105-9-559
- **HMMER**: Eddy SR (2011) Accelerated Profile HMM Searches. *PLoS Computational Biology* 7(10):e1002195. doi:10.1371/journal.pcbi.1002195

```bibtex
@article{zheng2016itak,
  title={iTAK: a program for genome-wide prediction and classification of plant transcription factors, transcriptional regulators, and protein kinases},
  author={Zheng, Yuan and Jiao, Chen and Sun, Honghe and Rosli, Hernan G and Pombo, Marina A and Zhang, Pengfei and Banf, Michael and Dai, Xinbin and Martin, Gregory B and Giovannoni, James J and others},
  journal={Molecular Plant},
  volume={9},
  number={12},
  pages={1667--1670},
  year={2016}
}
@article{langfelder2008wgcna,
  title={WGCNA: an R package for weighted correlation network analysis},
  author={Langfelder, Peter and Horvath, Steve},
  journal={BMC Bioinformatics},
  volume={9},
  pages={559},
  year={2008}
}
```

## 中文使用说明

### 功能
对 WGCNA 模块基因进行 iTAK 转录因子注释，并按 GO 富集通路绘制 TF-结构基因共表达网络图。

### 使用步骤
1. 准备模块基因蛋白序列（从基因组注释提取）
2. 运行 hmmscan 扫描 iTAK 数据库（Tfam_domain + TF_selfbuild）
3. 运行分类脚本（`itak_classify.py`）得到 TF 注释表
4. 配置通路列表和文件路径，运行绘图脚本（`plot_tf_network.py`）

### 参数调整
- `N_TOP_TF`: 每个通路补充的模块 TF 数量（默认 5）
- `r_inner/r_outer`: 内外圈半径比（默认 1.0/2.2）
- `PATHWAYS`: 通路名称字典（key=输出文件名后缀, value=GO Description 精确匹配）

### 注意事项
- GO 文件中 geneID 列用 `/` 分隔
- 边文件格式：fromNode\ttoNode\tweight\tdirection\tfromAltName\ttoAltName
- 蛋白 ID 带 `.1` 后缀，基因 ID 不带
