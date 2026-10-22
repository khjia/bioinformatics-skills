#!/usr/bin/env python3
"""把 100 kb 分块按「轮转」均分成 N 个工作单元（N = 想要的并发进程数）。

轮转分配的理由：相邻块耗时接近（细胞器除外），轮转比连续切分更均衡。
核基因组 -> unit_01..unit_N；细胞器单独成 unit_organelle.txt（深度极高，另行 --skip-coverage 处理）。

用法: python3 make_work_units.py [n_units]     # 惯例 40
"""
import os, sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
P = "<PROJECT_ROOT>"
V = f"{P}/req/05.variants"
CH = f"{V}/regions/chunks"
OUT = f"{V}/lists/work_units"
os.makedirs(OUT, exist_ok=True)

chroms = [l.split("\t")[0] for l in open(f"{CH}/index.tsv").readlines()[1:]]
nuc, org = [], []
for c in chroms:
    chunks = [l.strip() for l in open(f"{CH}/{c}.regions") if l.strip()]
    (org if not c.startswith("chr") else nuc).extend(chunks)

units = [[] for _ in range(N)]
for i, r in enumerate(nuc):
    units[i % N].append(r)
lines = []
for i, u in enumerate(units, 1):
    p = f"{OUT}/unit_{i:02d}.txt"
    open(p, "w").write("\n".join(u) + "\n")
    lines.append(p)
open(f"{OUT}/unit_organelle.txt", "w").write("\n".join(org) + "\n")
open(f"{V}/lists/workunits.txt", "w").write("\n".join(lines) + "\n")
sizes = [len(u) for u in units]
print(f"核块={len(nuc)} 单元={N} 每单元 {min(sizes)}-{max(sizes)} 块；细胞器 {len(org)} 块单列")
