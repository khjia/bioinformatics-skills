#!/usr/bin/env python3
"""
Extract upstream promoter regions (default 2kb) for a list of genes from GFF.
Output: BED file with strand-aware upstream coordinates.

Usage:
  python extract_promoters.py --gff <gff> --gene-list <ids.txt> \
         --upstream 2000 --out promoters.bed
Then:
  bedtools getfasta -fi genome.fa -bed promoters.bed -s -name+ -fo /dev/stdout \
    | sed 's/::.*$//' > promoters.fa
"""
import argparse
import re

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--gene-list', required=True, help='one gene ID per line')
    ap.add_argument('--upstream', type=int, default=2000)
    ap.add_argument('--out', required=True, help='output BED file')
    ap.add_argument('--feature', default='gene', help='GFF feature type (default: gene)')
    ap.add_argument('--id-attr', default='ID', help='attribute key for gene ID (default: ID)')
    args = ap.parse_args()

    with open(args.gene_list) as f:
        genes = {ln.strip() for ln in f if ln.strip()}

    rx = re.compile(rf'{args.id_attr}=([^;]+)')
    n = 0
    with open(args.gff) as f, open(args.out, 'w') as out:
        for ln in f:
            if ln.startswith('#') or not ln.strip():
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != args.feature:
                continue
            m = rx.search(cols[8])
            if not m:
                continue
            gid = m.group(1).split(':')[-1]  # strip gene: prefix if any
            if gid not in genes:
                continue
            chrom, start, end, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]
            if strand == '+':
                p_end = start - 1
                p_start = max(0, p_end - args.upstream)
            else:
                p_start = end
                p_end = end + args.upstream
            if p_end <= p_start:
                continue
            out.write(f"{chrom}\t{p_start}\t{p_end}\t{gid}\t.\t{strand}\n")
            n += 1
    print(f"[extract_promoters] {n} promoter intervals written to {args.out}")

if __name__ == '__main__':
    main()
