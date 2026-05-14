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

## Step 3: TFBS Validation (FIMO)

### Purpose
区分 TF→结构基因边是否有转录因子结合位点（TFBS）证据。有证据=粗实线，无证据=细虚线。

### Prerequisites
- MEME Suite (conda env `meme`, version 5.5.9+)
- JASPAR 植物核心 motif 库: `/media/nfs1/hermes/db/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt`
- 基因组 FASTA + FAI 索引
- GFF 注释文件
- bedtools, samtools

### Workflow
```bash
# 1. 提取结构基因上游 2kb 启动子区域 (BED)
# GFF gene行 → 根据strand取上游2kb → BED文件

# 2. 用 bedtools 提取序列
bedtools getfasta -fi genome.fa -bed promoters.bed -s -name+ -fo /dev/stdout | \
  sed 's/::.*$//' > promoters_2kb_clean.fa

# 3. FIMO 扫描 (p < 1e-4)
mamba run -n meme fimo --thresh 1e-4 --oc fimo_output \
  /media/nfs1/hermes/db/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt \
  promoters_2kb_clean.fa
```

### TF Family → Motif Mapping
JASPAR motif alt_id 包含 TF 蛋白名，通过关键词匹配到 iTAK 家族：
```python
family_patterns = {
    'AP2/ERF': ['ERF', 'RAP2', 'DREB', 'CBF'],
    'WRKY': ['WRKY'],
    'MYB': ['MYB', 'GAMYB'],
    'NAC': ['NAC', 'ANAC', 'CUC', 'NAM', 'ATAF', 'VND', 'SND'],
    'bZIP': ['bZIP', 'TGA', 'ABF', 'AREB', 'GBF', 'HY5', 'ABI5', 'EmBP'],
    'bHLH': ['bHLH', 'MYC', 'PIF', 'BEE', 'ICE', 'FBH', 'EIL'],
    'C2C2-Dof': ['DOF', 'Dof', 'PBF', 'OBP', 'CDF'],
    'C2H2': ['ZAT', 'AZF', 'STZ'],
    'HSF': ['HSF', 'HSFA', 'HSFB'],
    'TCP': ['TCP'],
    'HD-ZIP': ['ATHB', 'HAT', 'PHV', 'REV'],
    'ARF': ['ARF'],
    'MADS': ['MADS', 'AGL', 'SOC', 'SVP', 'AP1', 'AP3', 'PI', 'AG', 'SEP', 'squamosa'],
    'SBP': ['SPL', 'SBP'],
    'GATA': ['GATA'],
}
```

### Edge Classification Logic
```python
# 对每条 TF→结构基因边:
tf_fam = tf_family[tf_gene]           # 该TF的iTAK家族
struct_families = tfbs_evidence[struct_gene]  # 该结构基因启动子上命中的TF家族列表
has_tfbs = tf_fam in struct_families   # 家族匹配 = 有结合位点证据
```

### Visualization
- **有 TFBS 证据**: 粗红实线 (#E64B35, width×2.0, alpha=0.7)
- **无 TFBS 证据**: 细灰虚线 (#999999, width×0.6, alpha=0.4, dashed)
- **结构基因间**: 细青实线 (#91D1C2, width×0.7, alpha=0.3)

### Pitfalls
- bedtools getfasta `-name+` 输出 header 含 `::chr:start-end(strand)`，FIMO 会截断；需 `sed 's/::.*$//'` 清理
- FIMO 对 39 个 2kb 序列 × 805 motif 约需 90s
- p-value 阈值 1e-4 是 FIMO 默认，足够严格

## Dependencies
- Python: pandas, numpy, networkx, matplotlib
- System: hmmscan (HMMER3), bedtools, samtools, fimo (MEME Suite)
- Database: iTAK-db-v1, JASPAR2024 plant core motifs
- Conda env: `meme` (MEME Suite 5.5.9)

## Output
- Per-pathway PNG + PDF (300 dpi)
- TF annotation TSV: gene, gene_id, TF_family, type, full_family

## Citations

- **iTAK**: Zheng Y, et al. (2016) iTAK: A Program for Genome-wide Prediction and Classification of Plant Transcription Factors, Transcriptional Regulators, and Protein Kinases. *Molecular Plant* 9(12):1667-1670. doi:10.1016/j.molp.2016.09.014
- **WGCNA**: Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559. doi:10.1186/1471-2105-9-559
- **HMMER**: Eddy SR (2011) Accelerated Profile HMM Searches. *PLoS Computational Biology* 7(10):e1002195. doi:10.1371/journal.pcbi.1002195
- **MEME Suite / FIMO**: Grant CE, Bailey TL, Noble WS (2011) FIMO: scanning for occurrences of a given motif. *Bioinformatics* 27(7):1017-1018. doi:10.1093/bioinformatics/btr064
- **JASPAR**: Castro-Mondragon JA, et al. (2022) JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. *Nucleic Acids Research* 50(D1):D165-D173. doi:10.1093/nar/gkab1113
- **bedtools**: Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26(6):841-842. doi:10.1093/bioinformatics/btq033

Methods template:
> TF-target regulatory relationships were validated by scanning 2-kb upstream promoter regions of structural genes for TF binding sites using FIMO (v5.5.9, p < 1e-4) against the JASPAR 2024 plant core motif database (805 non-redundant motifs). TF family-motif correspondence was established by keyword matching between iTAK-classified TF families and JASPAR motif identifiers.

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
@article{grant2011fimo,
  title={FIMO: scanning for occurrences of a given motif},
  author={Grant, Charles E and Bailey, Timothy L and Noble, William Stafford},
  journal={Bioinformatics},
  volume={27},
  number={7},
  pages={1017--1018},
  year={2011}
}
@article{castro2022jaspar,
  title={JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles},
  author={Castro-Mondragon, Jaime A and others},
  journal={Nucleic Acids Research},
  volume={50},
  number={D1},
  pages={D165--D173},
  year={2022}
}
@article{quinlan2010bedtools,
  title={BEDTools: a flexible suite of utilities for comparing genomic features},
  author={Quinlan, Aaron R and Hall, Ira M},
  journal={Bioinformatics},
  volume={26},
  number={6},
  pages={841--842},
  year={2010}
}
```

## 中文使用说明

### 功能
对 WGCNA 模块基因进行 iTAK 转录因子注释，FIMO 扫描启动子 TFBS，并按 GO 富集通路绘制带 TFBS 证据的 TF-结构基因共表达网络图。

### 完整流程
1. 准备模块基因蛋白序列（从基因组注释提取）
2. 运行 hmmscan 扫描 iTAK 数据库（Tfam_domain + TF_selfbuild）
3. 运行分类脚本（`itak_classify.py`）得到 TF 注释表
4. 提取结构基因上游 2kb 启动子（`extract_promoters.py` + bedtools）
5. FIMO 扫描 JASPAR 植物 motif 库（conda env `meme`）
6. 构建 TFBS 证据 JSON（`build_tfbs_evidence.py`）
7. 配置通路 JSON，运行绘图脚本（`plot_tf_network_tfbs.py`）

### 快速命令
```bash
# Step 4-6: TFBS 验证
python extract_promoters.py --gff genome.gff --gene-list struct_genes.txt --out promoters.bed
bedtools getfasta -fi genome.fa -bed promoters.bed -s -name+ | sed 's/::.*$//' > promoters.fa
mamba run -n meme fimo --thresh 1e-4 --oc fimo_out JASPAR_plants.meme promoters.fa
python build_tfbs_evidence.py --fimo fimo_out/fimo.tsv --out tfbs_evidence.json

# Step 7: 绘图
echo '{"salt_response":"response to salt stress"}' > pathways.json
python plot_tf_network_tfbs.py --edges edges.txt --go go.txt --tf tf.tsv \
       --tfbs tfbs_evidence.json --pathways pathways.json --out-dir figures/
```

### 参数调整
- `--top-n`: 每个通路补充的模块 TF 数量（默认 5）
- `--pvalue`: FIMO p-value 阈值（默认 1e-4）
- `--upstream`: 启动子长度（默认 2000bp）
- `--pathways`: JSON 字典 {输出文件名后缀: GO Description 精确匹配}

### 注意事项
- GO 文件中 geneID 列用 `/` 分隔
- 边文件格式：fromNode\ttoNode\tweight\tdirection\tfromAltName\ttoAltName
- 蛋白 ID 带 `.1` 后缀，基因 ID 不带
- bedtools getfasta `-name+` 输出含 `::chr:start-end`，需 sed 清理
- JASPAR motif 库路径：`/media/nfs1/hermes/db/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt`
