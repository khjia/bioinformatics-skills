#!/usr/bin/env python3
"""
iTAK TF classification from hmmscan domtblout results.
Usage: Edit paths at bottom, then run: python3 itak_classify.py

Input:
  - hmmscan domtblout vs Tfam_domain.hmm
  - hmmscan domtblout vs TF_selfbuild.hmm
  - TF_Rule.txt from iTAK database
  - HMM files (for NAME->ACC mapping)

Output:
  - TSV: gene, gene_id, TF_family, type, full_family
"""
from collections import defaultdict, Counter
import re


def build_name_to_acc(hmm_file):
    """Build NAME -> ACC (Pfam ID) mapping from HMM file."""
    mapping = {}
    name = None
    with open(hmm_file) as f:
        for line in f:
            if line.startswith('NAME'):
                name = line.split()[1]
            elif line.startswith('ACC') and name:
                acc = line.split()[1].split('.')[0]
                mapping[name] = acc
                name = None
            elif line.startswith('//'):
                if name:
                    mapping[name] = name
                    name = None
    if name and name not in mapping:
        mapping[name] = name
    return mapping


def parse_domtbl(fn, name_to_acc):
    """Return {gene: {acc: count}}"""
    hits = defaultdict(lambda: defaultdict(int))
    with open(fn) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            domain_name = parts[0]
            gene = parts[3]
            acc = name_to_acc.get(domain_name, domain_name)
            hits[gene][acc] += 1
    return hits


def parse_rules(fn):
    rules = []
    current = {}
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if line == '//':
                if current:
                    rules.append(current)
                current = {}
                continue
            if ':' in line:
                key, val = line.split(':', 1)
                current[key.strip()] = val.strip()
    if current:
        rules.append(current)
    return rules


def parse_domain_requirement(req_str):
    """Parse 'PF00847#2' or 'PF00847#1--PF02362#1' or 'G2-like#1:PF00046#1'
    '--' = AND, ':' = OR. Returns list of AND groups, each = list of (id, min_count) OR alternatives.
    """
    if not req_str or req_str == 'NA':
        return []
    and_parts = req_str.split('--')
    result = []
    for part in and_parts:
        or_parts = part.split(':')
        or_group = []
        for op in or_parts:
            op = op.strip().strip('()')
            for item in op.split(','):
                item = item.strip()
                m = re.match(r'(.+?)#(\d+)$', item)
                if m:
                    or_group.append((m.group(1), int(m.group(2))))
        if or_group:
            result.append(or_group)
    return result


def check_requirement(gene_domains, req_groups):
    """Check if gene_domains satisfies requirement groups (AND of OR)"""
    if not req_groups:
        return True
    for or_group in req_groups:
        satisfied = False
        for pfam, min_count in or_group:
            if gene_domains.get(pfam, 0) >= min_count:
                satisfied = True
                break
        if not satisfied:
            return False
    return True


def classify_tf(tfam_domtbl, selfbuild_domtbl, tfam_hmm, selfbuild_hmm, rule_file, output_tsv):
    """Main classification pipeline."""
    print("[1] Building NAME->ACC mapping...")
    n2a = {}
    n2a.update(build_name_to_acc(tfam_hmm))
    n2a.update(build_name_to_acc(selfbuild_hmm))
    print(f"  Mapping: {len(n2a)} entries")

    print("[2] Parsing hmmscan results...")
    tfam = parse_domtbl(tfam_domtbl, n2a)
    selfbuild = parse_domtbl(selfbuild_domtbl, n2a)

    all_hits = defaultdict(lambda: defaultdict(int))
    for gene, doms in tfam.items():
        for d, c in doms.items():
            all_hits[gene][d] = max(all_hits[gene][d], c)
    for gene, doms in selfbuild.items():
        for d, c in doms.items():
            all_hits[gene][d] = max(all_hits[gene][d], c)
    print(f"  Genes with domain hits: {len(all_hits)}")

    print("[3] Parsing TF rules...")
    rules = parse_rules(rule_file)
    print(f"  Rules loaded: {len(rules)}")

    print("[4] Classifying...")
    tf_results = {}
    for gene, domains in all_hits.items():
        for rule in rules:
            req_str = rule.get('Required', 'NA')
            forb_str = rule.get('Forbidden', 'NA')
            req_groups = parse_domain_requirement(req_str)

            if check_requirement(domains, req_groups):
                has_forbidden = False
                if forb_str and forb_str != 'NA':
                    for forb_part in forb_str.replace('--', ':').split(':'):
                        forb_part = forb_part.strip().strip('()')
                        for item in forb_part.split(','):
                            item = item.strip()
                            m2 = re.match(r'(.+)#(\d+)$', item)
                            if m2 and domains.get(m2.group(1), 0) >= int(m2.group(2)):
                                has_forbidden = True
                                break
                        if has_forbidden:
                            break
                if has_forbidden:
                    continue
                tf_results[gene] = {
                    'family': rule.get('Name', ''),
                    'type': rule.get('Type', 'TF'),
                    'full_family': rule.get('Family', ''),
                }
                break

    print(f"  TF/TR classified: {len(tf_results)}")

    with open(output_tsv, 'w') as f:
        f.write('gene\tgene_id\tTF_family\ttype\tfull_family\n')
        for gene in sorted(tf_results.keys()):
            info = tf_results[gene]
            gene_id = gene.rsplit('.', 1)[0]
            f.write(f"{gene}\t{gene_id}\t{info['family']}\t{info['type']}\t{info['full_family']}\n")

    print("\nTop TF families:")
    fam_counts = Counter(v['family'] for v in tf_results.values())
    for fam, cnt in fam_counts.most_common(15):
        print(f"  {fam}: {cnt}")

    return tf_results


if __name__ == '__main__':
    # ---- EDIT THESE PATHS ----
    ITAK_DB = '/path/to/iTAK/database/iTAK-db-v1/database'
    DATA_DIR = '/path/to/project/data'

    classify_tf(
        tfam_domtbl=f'{DATA_DIR}/module_tfam.domtbl',
        selfbuild_domtbl=f'{DATA_DIR}/module_selfbuild.domtbl',
        tfam_hmm=f'{ITAK_DB}/Tfam_domain.hmm',
        selfbuild_hmm=f'{ITAK_DB}/TF_selfbuild.hmm',
        rule_file=f'{ITAK_DB}/TF_Rule.txt',
        output_tsv=f'{DATA_DIR}/module_itak_tf.tsv',
    )
