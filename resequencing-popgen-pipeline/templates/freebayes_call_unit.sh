#!/bin/bash
# freebayes 联合 calling：跑一个「工作单元」（= 一批 100 kb 区间，单元内顺序处理）。
# 一个单元 = 一个 freebayes 进程；并发数 = SGE 数组任务数（配 -pe smp k 控槽位）。
# 用法: bash freebayes_call_unit.sh <unit_number|unit_file>     # 数组里靠 SGE_TASK_ID
#   BAMLIST=<list>  OUT=<dir>  可覆盖；已存在的块 VCF 直接 SKIP（可断点续跑 / 重提安全）
#
# 前置（按顺序）:
#   python3 make_regions.py 100000        # -> regions/chunks/<chr>.regions
#   python3 make_work_units.py 40         # -> lists/work_units/unit_NN.txt + lists/workunits.txt
set -uo pipefail
export PATH=<CONDA_ENV_WITH_FREEBAYES>/bin:/media/nfs1/hermes/miniforge3/bin:$PATH
P=<PROJECT_ROOT>                  # 例: /media/nfs1/hermes/project/dadou
V=$P/req/05.variants              # 改成本项目对应的编号目录
OUT=${OUT:-$V/vcf/parts_raw}
LOG=$V/logs
BAMLIST=${BAMLIST:-$V/lists/bamlist.txt}
REF=$P/ref/genome.fasta

if [ $# -ge 1 ] && [ -f "$1" ]; then UF=$1
else N=${1:-${SGE_TASK_ID:-}}; [ -z "$N" ] && { echo "usage: freebayes_call_unit.sh <unit|file>"; exit 2; }
     UF=$(sed -n "${N}p" "$V/lists/workunits.txt"); fi
[ -s "$UF" ] || { echo "no unit file ($UF)"; exit 2; }
[ -s "$BAMLIST" ] || { echo "no bamlist $BAMLIST"; exit 2; }
mkdir -p "$OUT" "$LOG"

TAG=$(basename "$UF" .txt); n=$(wc -l < "$UF"); ns=$(wc -l < "$BAMLIST")
echo "[$(date '+%F %T')] START $TAG chunks=$n samples=$ns host=$(hostname)"
fail=0; done_n=0; t_start=$(date +%s)
while read -r r; do
    [ -z "$r" ] && continue
    out="$OUT/${r//:/-}.vcf"
    if [ -s "$out" ]; then echo "SKIP $r"; done_n=$((done_n+1)); continue; fi
    t0=$(date +%s)
    if freebayes -f "$REF" -r "$r" \
        -q 20 -m 20 --use-mapping-quality --genotype-qualities --ploidy 2 \
        --min-alternate-count 2 --min-alternate-fraction 0.05 \
        --bam-list "$BAMLIST" -v "$out.tmp" >/dev/null 2>>"$LOG/fb_err_$TAG.log"; then
        mv "$out.tmp" "$out"; done_n=$((done_n+1))
        echo "OK $r $(( $(date +%s) - t0 ))s $(grep -vc '^#' "$out") variants"
    else
        fail=$((fail+1)); rm -f "$out.tmp"; echo "FAIL $r"
    fi
done < "$UF"
echo "[$(date '+%F %T')] DONE $TAG ok=$done_n/$n fail=$fail wall=$(( $(date +%s) - t_start ))s"
[ "$fail" -eq 0 ]
