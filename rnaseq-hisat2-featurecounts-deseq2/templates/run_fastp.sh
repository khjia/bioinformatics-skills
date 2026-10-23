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
# duzhong 标准配置（2026-10 实测统一）
LEN_REQ=${LEN_REQ:-36}          # PE150 的 1/4；取样长的 1/3~1/2 即可
Q_PHRED=${Q_PHRED:-20}          # fastp 默认 15 几乎不过滤
U_LIMIT=${U_LIMIT:-40}          # 允许不合格碱基比例；Q30 < 90% 的项目设 20
CUT_MEAN_Q=${CUT_MEAN_Q:-20}
CUT_WINDOW=${CUT_WINDOW:-4}
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
        --cut_front --cut_tail \
        --cut_window_size "$CUT_WINDOW" \
        --cut_mean_quality "$CUT_MEAN_Q" \
        --qualified_quality_phred "$Q_PHRED" \
        --unqualified_percent_limit "$U_LIMIT" \
        --n_base_limit 5 \
        --average_qual 20 \
        --length_required "$LEN_REQ" \
        --thread "$THREADS" \
        -h "$LOGS/${sid}.fastp.html" \
        -j "$LOGS/${sid}.fastp.json" \
        2> "$LOGS/${sid}.fastp.log" \
      || { echo "[FAIL] $sid"; rm -f "$out1" "$out2"; }
done < "$SHEET"

echo "done -> $CLEAN"
