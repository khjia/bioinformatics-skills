#!/usr/bin/env python3
"""Create a deterministic synthetic 11-pair dataset; no biological sequences."""
import argparse
import json
from pathlib import Path
import random


def write(path, lines):
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def generate(folder):
    folder.mkdir(parents=True, exist_ok=False)
    rng = random.Random(17)
    pairs = [(f'chr{i:02}', f'A_chr{i:02}', f'B_chr{i:02}') for i in range(1, 12)]
    lengths = {chrom: 12_050_000 + rng.randrange(17_000_000) for _, a, b in pairs for chrom in (a, b)}
    config = {'assemblies': [], 'pairs': 'pairs.tsv', 'paf': 'B_to_A.paf', 'window_bp': 1_000_000}
    write(folder / 'pairs.tsv', ['label\tfirst\tsecond'] + ['\t'.join(pair) for pair in pairs])
    for index, side in enumerate(('A', 'B'), 1):
        files = {'label': f'Assembly {side}', 'fai': f'{side}.fa.fai', 'genes': f'{side}.genes.gff3',
                 'repeats': f'{side}.repeats.gff3', 'telomeres': f'{side}.ends.tsv', 'rdna': f'{side}.rdna.tsv'}
        config['assemblies'].append(files)
        fai, genes, repeats = [], ['##gff-version 3'], ['##gff-version 3']
        ends = ['chromosome\tend\tstatus']
        rdna = ['chromosome\tstart\tend\ttype']
        for n, pair in enumerate(pairs):
            chrom = pair[index]
            length = lengths[chrom]
            fai.append(f'{chrom}\t{length}\t0\t80\t81')
            ends.append(f'{chrom}\tleft\t{"detected" if n % 5 else "not_detected"}')
            if n % 3:
                ends.append(f'{chrom}\tright\tdetected')
            if n in (1, 7):
                rdna.append(f'{chrom}\t3050000\t3180000\t5S')
            if n in (3, 7):
                rdna.append(f'{chrom}\t4100000\t4230000\t18S')
            for i, start in enumerate(range(1, length - 490_000, 350_000)):
                for gene in range(rng.randrange(7)):
                    pos = start + gene * 12_000
                    genes.append(f'{chrom}\tdemo\tgene\t{pos}\t{pos+2780}\t.\t+\t.\tID=g{n}_{i}_{gene}')
                for j, repeat in enumerate(('LTR/Gypsy', 'LTR/Copia', 'DNA/hAT', 'LINE/L1')):
                    pos = start + j * 72_000
                    repeats.append(f'{chrom}\tdemo\trepeat_region\t{pos}\t{pos+rng.randrange(4_000, 140_000)}\t.\t.\t.\tClassification={repeat}')
        for key, rows in (('fai', fai), ('genes', genes), ('repeats', repeats), ('telomeres', ends), ('rdna', rdna)):
            write(folder / files[key], rows)
    paf = []
    for i, (_, a, b) in enumerate(pairs):
        for j in range(4):
            s1, s2, span = j * 2_800_000 + 100_000, j * 2_600_000 + 170_000, 2_100_000
            strand = '-' if i % 4 == 0 and j == 1 else '+'
            paf.append(f'{b}\t{lengths[b]}\t{s2}\t{s2+span}\t{strand}\t{a}\t{lengths[a]}\t{s1}\t{s1+span}\t{span-750}\t{span}\t60\ttp:A:P')
    write(folder / 'B_to_A.paf', paf)
    write(folder / 'manifest.json', [json.dumps(config, indent=2)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    generate(parser.parse_args().directory)
