#!/bin/bash
# 批量 HISAT2 比对
# 用法: bash run_hisat2.sh <sample_sheet.tsv>
#   sample_sheet.tsv 每行: <sample_id> <R1.fastq.gz> <R2.fastq.gz>
#
# MODE=sam_nfs  (默认)  SAM 写 NFS 的 tmp_sam/ 再排序 —— NFS 有空间、/tmp 小的场景
# MODE=lowio    管道直排 sort，不落 SAM —— NFS 紧张、追求最低 IO 的场景
#
# ⚠️ 开工前必须先看两个盘:
#     df -h /media/nfs2   # 输出盘
#     df -h /tmp          # 本地临时盘（大基因组每样本 SAM 10-21 GB，并行 N 个会瞬间吃满 /tmp）
#     exit 134 有两个原因（NFS 满 / /tmp 满），修错方向会二次失败。
set -euo pipefail

PROJ=${PROJ:?set PROJ=/path/to/project}
MAPPING=$PROJ/3.Mapping
CLEAN=$PROJ/2.clean
BAM=$MAPPING/bam
SUMMARY=$MAPPING/summary
LOGS=$MAPPING/logs
REF=${REF:-$MAPPING/ref/genome}          # ← index 前缀，不是 fasta 路径！
CLEAN_SUFFIX=${CLEAN_SUFFIX:-_clean}
MODE=${MODE:-sam_nfs}
THREADS=${THREADS:-16}
SORT_THREADS=${SORT_THREADS:-4}
SORT_MEM=${SORT_MEM:-4G}
PARALLEL=${PARALLEL:-2}

mkdir -p "$BAM" "$SUMMARY" "$LOGS"
[ -d "$MAPPING/tmp_sam" ] || mkdir -p "$MAPPING/tmp_sam"

SHEET=${1:?usage: run_hisat2.sh <sample_sheet.tsv>}

map_one() {
    local sid=$1 r1=$2 r2=$3
    local bam=$BAM/${sid}.sorted.bam
    local log=$LOGS/${sid}.hisat2.log

    # 幂等：已完成直接跳过
    if [ -s "$bam" ] && [ -s "${bam}.csi" ]; then echo "[skip] $sid"; return 0; fi

    if [ "$MODE" = "lowio" ]; then
        hisat2 -p "$THREADS" --no-unal -x "$REF" \
            -1 "$r1" -2 "$r2" 2> "$log" \
          | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "/tmp/whgsort_${sid}" -o "$bam" -
    else
        local sam=$MAPPING/tmp_sam/${sid}.sam
        hisat2 -p "$THREADS" -x "$REF" -1 "$r1" -2 "$r2" -S "$sam" 2> "$log"
        samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -o "$MAPPING/tmp_sam/${sid}.sorted.bam" "$sam"
        rm -f "$sam"
        mv "$MAPPING/tmp_sam/${sid}.sorted.bam" "$bam"
    fi

    # ⚠️ post-mapping guard：样本名带括号等会导致 0-read BAM 且 exit 0，必须显式拦截
    local pairs
    pairs=$(grep -oP '^\d+(?= reads; of these:)' "$log" || true)
    if [ -z "${pairs:-}" ] || [ "${pairs:-0}" -eq 0 ] || [ "$(stat -c%s "$bam")" -lt 1000000 ]; then
        echo "[FAIL_EMPTY] $sid (pairs=${pairs:-0}) —— 检查样本名是否含括号/空格"
        rm -f "$bam" "${bam}.csi"
        return 1
    fi

    samtools index -c "$bam"
    samtools flagstat "$bam" > "$SUMMARY/${sid}.flagstat"
    echo "[OK] $sid pairs=$pairs"
}
export -f map_one
export BAM SUMMARY LOGS REF MAPPING CLEAN MODE THREADS SORT_THREADS SORT_MEM

grep -v '^#' "$SHEET" | awk 'NF>=3' \
  | xargs -P "$PARALLEL" -I{} bash -c 'map_one {}'

echo "done -> $BAM"
