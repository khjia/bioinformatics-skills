#!/bin/bash
# Batch query ENA Portal API for library metadata of multiple runs.
# Usage: bash check_ena_strategy.sh <runs.txt> <output.tsv>
#
# Input: one run accession per line
# Output: TSV with columns: run strategy source layout platform

set -euo pipefail

INPUT=${1:?need input file (one run per line)}
OUT=${2:?need output TSV path}
PROXY=${PROXY:-http://127.0.0.1:10809}
CHUNK=${CHUNK:-40}  # ENA URL length limit ~2000 chars, ~40 runs safe

echo -e "run\tstrategy\tsource\tlayout\tplatform" > "$OUT"
total=$(wc -l < "$INPUT")
offset=0

while [ $offset -lt $total ]; do
  runs=$(tail -n +$((offset+1)) "$INPUT" | head -$CHUNK | tr '\n' ',' | sed 's/,$//')
  query="run_accession=$(echo "$runs" | sed 's/,/%20OR%20run_accession=/g')"

  curl -sS -x "$PROXY" \
    "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=${query}&fields=run_accession,library_strategy,library_source,library_layout,instrument_platform&format=tsv&limit=$CHUNK" \
    | tail -n +2 >> "$OUT"

  offset=$((offset + CHUNK))
  echo "Progress: $offset / $total" >&2
done

echo "Done. Lines: $(wc -l < "$OUT")" >&2
