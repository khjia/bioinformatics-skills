#!/bin/bash
# HISAT2 建索引：GFF3 -> GTF -> splice_sites/exons -> index
# 用法: bash build_hisat2_index.sh <genome.fa> <annotation.gff3> [outdir]
# 关键点:
#   1) GFF3 必须先转 GTF（否则 splice_sites/exons 是 0 字节）
#   2) python 必须在 PATH（CNS 上默认只有 python3，否则索引"瞬间完成"但是空文件）
#   3) NFS >95% 满时索引静默损坏 -> 先在 /tmp 建好再 cp
set -euo pipefail

export PATH=/media/nfs1/hermes/miniforge3/bin:$PATH   # 保证 python 可解析

FASTA=${1:?usage: build_hisat2_index.sh <genome.fa> <annotation.gff3> [outdir]}
GFF=${2:?usage: build_hisat2_index.sh <genome.fa> <annotation.gff3> [outdir]}
OUT=${3:-$(dirname "$FASTA")}
PREFIX=${PREFIX:-$(basename "$FASTA" .fa)}
THREADS=${THREADS:-32}

mkdir -p "$OUT"
GTF=$OUT/${PREFIX}.gtf

# 1) GFF3 -> GTF
if [ ! -s "$GTF" ]; then
    echo "[gffread] $GFF -> $GTF"
    gffread "$GFF" -T -o "$GTF"
fi

# 2) 提取剪接位点与外显子
hisat2_extract_splice_sites.py "$GTF" > "$OUT/${PREFIX}.splice_sites.txt"
hisat2_extract_exons.py      "$GTF" > "$OUT/${PREFIX}.exons.txt"
ns=$(wc -l < "$OUT/${PREFIX}.splice_sites.txt")
ne=$(wc -l < "$OUT/${PREFIX}.exons.txt")
echo "splice_sites=$ns exons=$ne"
if [ "$ns" -eq 0 ] || [ "$ne" -eq 0 ]; then
    echo "[ERROR] 注释为空 —— 检查 GFF3 是否真的转成了带 gene_id/transcript_id 的 GTF" >&2
    exit 1
fi

# 3) 建索引（本地 /tmp 建，再拷贝，规避 NFS 满导致的静默损坏）
TMPDIR_IX=/tmp/hisat2_build_$$
mkdir -p "$TMPDIR_IX"
trap 'rm -rf "$TMPDIR_IX"' EXIT
hisat2-build -p "$THREADS" \
    --ss "$OUT/${PREFIX}.splice_sites.txt" \
    --exon "$OUT/${PREFIX}.exons.txt" \
    "$FASTA" "$TMPDIR_IX/${PREFIX}"

# 4) 拷贝 + 校验
cp "$TMPDIR_IX"/${PREFIX}.*.ht2 "$OUT/"
ls -la "$OUT"/${PREFIX}.*.ht2
# .3.ht2 / .7.ht2 / .8.ht2 天然很小（12-152 字节），只要非空即可
for i in 1 2 3 4 5 6 7 8; do
    f=$OUT/${PREFIX}.${i}.ht2
    [ -s "$f" ] || { echo "[ERROR] 缺 $f" >&2; exit 1; }
done
echo "[OK] index -> $OUT/${PREFIX}"
