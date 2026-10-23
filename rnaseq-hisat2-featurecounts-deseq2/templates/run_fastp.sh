#!/bin/bash
# PE RNA-seq 批量质控（fastp）
# 用法: bash run_fastp.sh <sample_sheet.tsv>
#   sample_sheet.tsv 每行: <sample_id> <R1.fastq.gz> <R2.fastq.gz>
# 输出: <CLEAN>/<sid>_clean_R1.fastq.gz / _R2.fastq.gz + logs/<sid>.fastp.{html,json}
set -euo pipefail

PROJ=${PROJ:?set PROJ=/path/to/project}
CLEAN=$PROJ/2.clean
LOGS=$CLEAN/logs
THREADS=${THREADS:-8}
LEN_REQ=${LEN_REQ:-30}
Q_PHRED=${Q_PHRED:-15}
mkdir -p "$CLEAN" "$LOGS"

SHEET=${1:?usage: run_fastp.sh <sample_sheet.tsv>}

while read -r sid r1 r2; do
    [ -z "${sid:-}" ] && continue
    case "$sid" in \#*) continue;; esac
    out1=$CLEAN/${sid}_clean_R1.fastq.gz
    out2=$CLEAN/${sid}_clean_R2.fastq.gz
    # 幂等: 已完成则跳过
    if [ -s "$out1" ] && [ -s "$out2" ] && [ -s "$LOGS/${sid}.fastp.json" ]; then
        echo "[skip] $sid"; continue
    fi
    echo "[$(date '+%F %T')] fastp $sid"
    fastp \
        -i "$r1" -I "$r2" \
        -o "$out1" -O "$out2" \
        --detect_adapter_for_pe \
        --length_required "$LEN_REQ" \
        --qualified_quality_phred "$Q_PHRED" \
        --thread "$THREADS" \
        -h "$LOGS/${sid}.fastp.html" \
        -j "$LOGS/${sid}.fastp.json" \
        2> "$LOGS/${sid}.fastp.log" \
      || { echo "[FAIL] $sid"; rm -f "$out1" "$out2"; }
done < "$SHEET"

echo "done -> $CLEAN"
