#!/usr/bin/env python3
"""Select the best DIAMOND hit per query after strict identity/E-value filters.

Expected input: tabular DIAMOND output with columns in this order:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send
evalue bitscore qlen slen stitle
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "stitle",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Raw DIAMOND TSV")
    parser.add_argument("--output", required=True, type=Path, help="Best-hit TSV")
    parser.add_argument("--database", required=True, help="Database label for output")
    parser.add_argument("--min-identity", type=float, default=30.0,
                        help="Strict lower identity bound (default: >30 percent)")
    parser.add_argument("--max-evalue", type=float, default=1e-5,
                        help="Strict upper E-value bound (default: <1e-5)")
    args = parser.parse_args()

    best: dict[str, tuple[tuple[float, float, float], list[str]]] = {}
    rows = 0
    with args.input.open(encoding="utf-8") as source:
        for line_no, line in enumerate(source, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t", len(FIELDS) - 1)
            if len(fields) != len(FIELDS):
                raise SystemExit(
                    f"{args.input}:{line_no}: expected {len(FIELDS)} tab-separated fields, "
                    f"found {len(fields)}"
                )
            try:
                identity = float(fields[2])
                evalue = float(fields[10])
                bitscore = float(fields[11])
            except ValueError as exc:
                raise SystemExit(f"{args.input}:{line_no}: invalid numeric field: {exc}") from exc
            rows += 1
            if identity <= args.min_identity or evalue >= args.max_evalue:
                continue
            # Highest bit score wins; ties go to lower E-value, then higher identity.
            rank = (-bitscore, evalue, -identity)
            query = fields[0]
            if query not in best or rank < best[query][0]:
                best[query] = (rank, fields)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as output:
        output.write("database\t" + "\t".join(FIELDS) + "\n")
        for _, fields in best.values():
            output.write(args.database + "\t" + "\t".join(fields) + "\n")

    print(
        f"database={args.database} candidate_alignments={rows} "
        f"queries_with_passing_best_hit={len(best)} "
        f"identity>{args.min_identity} evalue<{args.max_evalue:g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
