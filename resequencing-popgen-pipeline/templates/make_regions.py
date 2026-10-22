#!/usr/bin/env python3
"""生成 freebayes 分块区间文件：每个序列一个 .regions（每行一块 chr:start-end）。

用法: python3 make_regions.py [chunk_bp]     # 默认 100000
输出: <V>/regions/chunks/<chr>.regions + index.tsv
随后: awk '{print $1}' regions/chunks/index.tsv | tail -n +2 > lists/chroms.txt
"""
import os, sys

CHUNK = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
P = "<PROJECT_ROOT>"          # 例: /media/nfs1/hermes/project/dadou
V = f"{P}/req/05.variants"
FAI = f"{P}/ref/genome.fasta.fai"
OUT = f"{V}/regions/chunks"
os.makedirs(OUT, exist_ok=True)

rows = []
for line in open(FAI):
    f = line.split("\t")
    chrom, ln = f[0], int(f[1])
    chunks, st = [], 1
    while st <= ln:
        en = min(st + CHUNK - 1, ln)
        chunks.append(f"{chrom}:{st}-{en}")
        st = en + 1
    with open(f"{OUT}/{chrom}.regions", "w") as fh:
        fh.write("\n".join(chunks) + "\n")
    rows.append((chrom, ln, len(chunks)))

with open(f"{V}/regions/chunks/index.tsv", "w") as fh:
    fh.write("chrom\tlength\tchunks\n")
    for c, l, n in rows:
        fh.write(f"{c}\t{l}\t{n}\n")
print(f"chunk={CHUNK} bp ; {len(rows)} seqs ; {sum(n for _,_,n in rows)} chunks -> {OUT}")
