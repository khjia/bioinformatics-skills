#!/bin/bash
# Speed sampling — read on-disk file growth (don't trust curl's stats).
#
# Usage: ./speed_check.sh /path/to/01.rawdata [sample_seconds]
set -u
DIR=${1:-$PWD}
SECS=${2:-60}
GLOB=${GLOB:-'*.fastq.gz'}

t1=/tmp/_speed_sz1.$$
t2=/tmp/_speed_sz2.$$
trap "rm -f $t1 $t2" EXIT

ls -l $DIR/$GLOB 2>/dev/null | awk '{print $9, $5}' > $t1
echo "Sampling for ${SECS}s..."
sleep "$SECS"
ls -l $DIR/$GLOB 2>/dev/null | awk '{print $9, $5}' > $t2

join -j 1 $t1 $t2 | awk -v secs="$SECS" '
{
  diff = $3 - $2
  if (diff > 0) {
    mb     = diff / 1048576
    speed  = mb / secs
    cur_mb = $3 / 1048576
    name=$1; sub(/.*\//,"",name)
    printf "  %-50s  +%6.1f MB  %.2f MB/s  (cur %.0f MB)\n",
           name, mb, speed, cur_mb
    total += diff
  }
}
END {
  if (total > 0) printf "\nAggregate: %.2f MB/s over %ds\n", total/1048576/secs, secs
  else           print  "\nNo growth observed."
}
'

echo
echo "Active curl processes: $(pgrep -c '^curl$' || echo 0)"
