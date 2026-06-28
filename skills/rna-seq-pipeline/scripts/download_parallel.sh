#!/bin/bash
# ENA parallel FASTQ downloader via sing-box proxy.
# Usage:
#   META=/path/to/metadata WORKLIST=/path/to/worklist.txt \
#     LOG=/path/to/download.log nohup bash download_parallel.sh > /dev/null 2>&1 &
#
# WORKLIST format (one record per line):
#   URL|DSTDIR|EXPECTED_MD5
#
# Env vars (with defaults):
#   PROXY=http://127.0.0.1:10809
#   CONCURRENCY=4
#   META=$PWD
#   LOG=$META/download.log
#   WORKLIST=$META/worklist.txt
set -u
META=${META:-$PWD}
LOG=${LOG:-$META/download.log}
WORKLIST=${WORKLIST:-$META/worklist.txt}
PROXY=${PROXY:-http://127.0.0.1:10809}
CONCURRENCY=${CONCURRENCY:-4}

dl_one() {
  local line=$1
  local url="${line%%|*}"
  local rest="${line#*|}"
  local dst="${rest%|*}"
  local md5="${rest##*|}"
  local fname; fname=$(basename "$url")

  mkdir -p "$dst"

  # Skip if already correct
  if [ -e "$dst/$fname" ]; then
    if [ "$(md5sum "$dst/$fname" | awk '{print $1}')" = "$md5" ]; then
      echo "[$(date +%H:%M:%S)] SKIP: $fname" >> "$LOG"
      return 0
    fi
    rm -f "$dst/$fname"
  fi

  for attempt in 1 2 3; do
    if curl -sS --max-time 7200 --retry 0 -x "$PROXY" -o "$dst/$fname" "$url"; then
      if [ "$(md5sum "$dst/$fname" | awk '{print $1}')" = "$md5" ]; then
        echo "[$(date +%H:%M:%S)] OK: $fname ($(du -h "$dst/$fname" | cut -f1)) attempt=$attempt" >> "$LOG"
        return 0
      else
        echo "[$(date +%H:%M:%S)] MD5_FAIL attempt=$attempt: $fname" >> "$LOG"
        rm -f "$dst/$fname"
      fi
    else
      echo "[$(date +%H:%M:%S)] CURL_FAIL attempt=$attempt: $fname" >> "$LOG"
    fi
    sleep 10
  done
  echo "[$(date +%H:%M:%S)] !!GIVE_UP: $fname" >> "$LOG"
  return 1
}
export -f dl_one
export PROXY LOG

[ -s "$WORKLIST" ] || { echo "WORKLIST empty: $WORKLIST" >&2; exit 1; }

echo "[$(date +%F\ %T)] START $(wc -l < "$WORKLIST") files concurrency=$CONCURRENCY" >> "$LOG"

# CRITICAL: pass {} as $1 via `_ {}`; do NOT use IFS=$'\t' read <<<"{}"
xargs -a "$WORKLIST" -P "$CONCURRENCY" -I {} bash -c 'dl_one "$@"' _ {}

echo "[$(date +%F\ %T)] DONE" >> "$LOG"

# Summary
ok=$(grep -c '^.*OK: ' "$LOG" || true)
giveup=$(grep -c '!!GIVE_UP' "$LOG" || true)
echo "[$(date +%F\ %T)] summary: OK=$ok GIVE_UP=$giveup" >> "$LOG"
