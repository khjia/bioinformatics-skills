#!/usr/bin/env python3
"""
Parse FIMO output + JASPAR motif metadata → per-gene TF family evidence JSON.

For each target gene (sequence_name in fimo.tsv), record the set of TF families
whose motifs hit its promoter. Family is inferred from motif alt_id keywords.

Usage:
  python build_tfbs_evidence.py --fimo fimo_output/fimo.tsv \
         --pvalue 1e-4 --out fimo_tf_binding_evidence.json
"""
import argparse
import json
from collections import defaultdict

# TF family → motif name keyword patterns (JASPAR plant core)
DEFAULT_FAMILY_PATTERNS = {
    'AP2/ERF':  ['ERF', 'RAP2', 'DREB', 'CBF'],
    'WRKY':     ['WRKY'],
    'MYB':      ['MYB', 'GAMYB'],
    'MYB-related': ['MYB-related', 'RVE', 'CCA1', 'LHY'],
    'NAC':      ['NAC', 'ANAC', 'CUC', 'NAM', 'ATAF', 'VND', 'SND'],
    'bZIP':     ['bZIP', 'TGA', 'ABF', 'AREB', 'GBF', 'HY5', 'ABI5', 'EmBP'],
    'bHLH':     ['bHLH', 'MYC', 'PIF', 'BEE', 'ICE', 'FBH', 'EIL'],
    'C2C2-Dof': ['DOF', 'Dof', 'PBF', 'OBP', 'CDF'],
    'C2H2':     ['ZAT', 'AZF', 'STZ'],
    'HSF':      ['HSF', 'HSFA', 'HSFB'],
    'TCP':      ['TCP'],
    'HD-ZIP':   ['ATHB', 'HAT', 'PHV', 'REV'],
    'ARF':      ['ARF'],
    'MADS':     ['MADS', 'AGL', 'SOC', 'SVP', 'AP1', 'AP3', 'PI', 'AG', 'SEP', 'squamosa'],
    'SBP':      ['SPL', 'SBP'],
    'GATA':     ['GATA'],
    'B3':       ['ABI3', 'FUS3', 'LEC2', 'VRN1'],
    'GRAS':     ['GRAS', 'SCR', 'SHR'],
    'Trihelix': ['GT-1', 'GT-2', 'GT-3'],
    'LBD':      ['LBD', 'ASL'],
    'CAMTA':    ['CAMTA'],
}

def classify_family(alt_id, patterns):
    if not alt_id:
        return None
    s = str(alt_id)
    for fam, keys in patterns.items():
        for k in keys:
            if k.lower() in s.lower():
                return fam
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fimo', required=True, help='FIMO fimo.tsv output')
    ap.add_argument('--pvalue', type=float, default=1e-4)
    ap.add_argument('--out', required=True, help='output JSON: gene -> [families]')
    ap.add_argument('--patterns-json', help='optional override of family patterns (JSON)')
    args = ap.parse_args()

    patterns = DEFAULT_FAMILY_PATTERNS
    if args.patterns_json:
        with open(args.patterns_json) as f:
            patterns = json.load(f)

    gene_fams = defaultdict(set)
    n_hits = n_kept = n_unknown = 0
    with open(args.fimo) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for ln in f:
            if not ln.strip() or ln.startswith('#'):
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < len(header):
                continue
            n_hits += 1
            try:
                p = float(cols[idx['p-value']])
            except (ValueError, KeyError):
                continue
            if p > args.pvalue:
                continue
            gene = cols[idx['sequence_name']]
            alt = cols[idx.get('motif_alt_id', idx.get('motif_id'))]
            fam = classify_family(alt, patterns)
            if fam is None:
                n_unknown += 1
                continue
            gene_fams[gene].add(fam)
            n_kept += 1

    out = {g: sorted(list(s)) for g, s in gene_fams.items()}
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=2)

    print(f"[build_tfbs] total hits: {n_hits}")
    print(f"[build_tfbs] kept (p<{args.pvalue}): {n_kept}")
    print(f"[build_tfbs] unknown family (skipped): {n_unknown}")
    print(f"[build_tfbs] genes with evidence: {len(out)}")
    print(f"[build_tfbs] written to {args.out}")

if __name__ == '__main__':
    main()
