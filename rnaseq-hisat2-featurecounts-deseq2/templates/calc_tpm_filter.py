#!/usr/bin/env python3
"""Calculate TPM from featureCounts output and filter by tissue-level mean TPM >= 1.

Usage:
    python3 calc_tpm_filter.py <counts_raw.txt> [--tissue-regex '^(?P<tissue>.+?)_rep\\d+$'] [--out-dir .]

The tissue regex must capture a group named `tissue`. Default handles sample names
like `Gli_root1_rep1` -> tissue `root1`.

Outputs:
    counts_matrix.tsv / tpm_matrix.tsv
    counts_matrix_filtered.tsv / tpm_matrix_filtered.tsv

Filter rule: a gene is kept if ANY tissue has mean TPM >= 1.
Using the global mean across all samples removes tissue-specific genes -- do not do that.
"""
import argparse
import os
import re


def get_tissue(sample_name, pattern):
    m = re.match(pattern, sample_name)
    if not m:
        return sample_name
    if "tissue" in (m.groupdict() or {}):
        return m.group("tissue")
    if m.groups():
        return m.group(1)
    return sample_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts_file")
    ap.add_argument("--tissue-regex", default=r"^(?P<tissue>.+?)_rep\d+$")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--min-tpm", type=float, default=1.0)
    args = ap.parse_args()

    with open(args.counts_file) as f:
        lines = [l for l in f if not l.startswith("#")]

    header = lines[0].strip().split("\t")
    # featureCounts columns: Chr Start End Strand Length <bam...>
    sample_names = [os.path.basename(s).replace(".sorted.bam", "")
                    for s in header[6:]]

    gene_data, gene_lengths = {}, {}
    for line in lines[1:]:
        cols = line.strip().split("\t")
        gene_id = cols[0]
        gene_lengths[gene_id] = int(cols[5])
        gene_data[gene_id] = [int(c) for c in cols[6:]]

    # TPM = (count / length_kb) * 1e6 / sum(count_j / length_kb_j)
    tpm = {}
    for gene_id, counts in gene_data.items():
        tpm[gene_id] = [c / (gene_lengths[gene_id] / 1000.0) for c in counts]

    for j in range(len(sample_names)):
        total_rpk = sum(tpm[g][j] for g in tpm)
        scale = 1e6 / total_rpk if total_rpk > 0 else 0
        for gene_id in tpm:
            tpm[gene_id][j] *= scale

    tissue_groups = {}
    for s in sample_names:
        tissue_groups.setdefault(get_tissue(s, args.tissue_regex), []).append(s)

    filtered = [g for g in sorted(tpm)
                if any(sum(tpm[g][sample_names.index(s)] for s in ts) / len(ts) >= args.min_tpm
                       for ts in tissue_groups.values())]

    os.makedirs(args.out_dir, exist_ok=True)
    def fmt_int(v):
        return str(v)

    def fmt_float(v):
        return f"{v:.4f}"

    jobs = [(gene_data, "counts_matrix", fmt_int, sorted(gene_data)),
            (tpm, "tpm_matrix", fmt_float, sorted(tpm)),
            (gene_data, "counts_matrix_filtered", fmt_int, filtered),
            (tpm, "tpm_matrix_filtered", fmt_float, filtered)]
    for src, suffix, formatter, genes in jobs:
        out = os.path.join(args.out_dir, f"{suffix}.tsv")
        with open(out, "w") as f:
            f.write("GeneID\t" + "\t".join(sample_names) + "\n")
            for gene_id in genes:
                f.write(gene_id + "\t" + "\t".join(
                    formatter(src[gene_id][j]) for j in range(len(sample_names))) + "\n")
        print(f"Written: {out}")

    print(f"Total: {len(gene_data)} -> Filtered: {len(filtered)} "
          f"(removed {len(gene_data) - len(filtered)})")
    print(f"Tissues: {', '.join(sorted(tissue_groups))}")


if __name__ == "__main__":
    main()
