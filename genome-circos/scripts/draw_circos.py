#!/usr/bin/env python3
"""Draw an interleaved, two-assembly genome overview and its English legend."""

import argparse
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re
import sys
from xml.sax.saxutils import escape

TRACKS = (
    ('Genes', '#2166AC'), ('LTR/Gypsy', '#862D2D'), ('LTR/Copia', '#17694B'),
    ('DNA transposons', '#915419'), ('LINE and unclassified repeats', '#613B91'),
)
# Chromosome accents are secondary; position and labels identify each chromosome.
ACCENTS = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300',
           '#6250d6', '#e34948', '#087f9e', '#98510a', '#7c449f')
TEL = {'detected': '#19744C', 'not_detected': '#BD382D', 'not_assessed': '#646E79'}
RDNA = {'5S': '#2563B8', '18S': '#CA6312'}
SIZE, CX, CY = 1080, 540, 440
INK = '#172333'


def input_text(path):
    return (gzip.open if path.suffix == '.gz' else open)(path, 'rt', encoding='utf-8')


def records(path):
    with input_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.strip() and not line.lstrip().startswith('#'):
                yield number, line.rstrip('\r\n').split('\t')


def table(path, header):
    with input_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        if next(reader, None) != header.split():
            raise ValueError(f'{path}: expected TSV header: {header}')
        for row in reader:
            if row:
                if len(row) != len(header.split()):
                    raise ValueError(f'{path}:{reader.line_num}: wrong column count')
                yield row


def interval(start, end, length, location):
    if not 0 <= start < end <= length:
        raise ValueError(f'{location}: interval [{start}, {end}) outside 0–{length} or empty')
    return start, end


def fai(path):
    lengths = {}
    for number, fields in records(path):
        if len(fields) < 2 or fields[0] in lengths or int(fields[1]) <= 0:
            raise ValueError(f'{path}:{number}: invalid FAI record')
        lengths[fields[0]] = int(fields[1])
    return lengths


def repeat_track(classification):
    value = classification.lower()
    if value.startswith('ltr/gypsy'):
        return 1
    if value.startswith('ltr/copia'):
        return 2
    if value.startswith(('dna/', 'mite/', 'rc/helitron')):
        return 3
    if value.startswith('line/') or value in ('line', 'unknown', 'ltr/unknown', 'sine/unknown'):
        return 4
    return None


def features(path, lengths, kind):
    with input_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            location = f'{path}:{number}'
            if len(fields) != 9:
                raise ValueError(f'{location}: expected nine GFF3 columns')
            chrom = fields[0]
            if chrom not in lengths:
                continue
            if kind == 'genes':
                if fields[2] != 'gene':
                    continue
                track = 0
            else:
                attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
                track = repeat_track(attrs.get('Classification', attrs.get('classification', '')))
                if track is None:
                    continue
            start, end = interval(int(fields[3]) - 1, int(fields[4]), lengths[chrom], location)
            yield chrom, track, start, end


def union_intervals(ranges):
    merged = []
    for start, end in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def densities(assembly, lengths, window):
    result = {c: [[0] * ((n + window - 1) // window) for _ in TRACKS] for c, n in lengths.items()}
    extents = defaultdict(list)
    for chrom, _, start, end in features(assembly['genes'], lengths, 'genes'):
        result[chrom][0][((start + end - 1) // 2) // window] += 1
    if not any(any(tracks[0]) for tracks in result.values()):
        raise ValueError(f'{assembly["genes"]}: no gene records match the selected chromosomes')
    for chrom, track, start, end in features(assembly['repeats'], lengths, 'repeats'):
        extents[chrom, track].append((start, end))
    if not extents:
        raise ValueError(f'{assembly["repeats"]}: no classified repeat records match the selected chromosomes')
    for (chrom, track), ranges in extents.items():
        for start, end in union_intervals(ranges):
            for index in range(start // window, (end - 1) // window + 1):
                result[chrom][track][index] += min(end, (index + 1) * window) - max(start, index * window)
    return result


def markers(assembly, lengths):
    ends, sites = {}, []
    if assembly.get('telomeres'):
        for chrom, end, status in table(assembly['telomeres'], 'chromosome end status'):
            if chrom not in lengths or end not in ('left', 'right') or status not in TEL or (chrom, end) in ends:
                raise ValueError(f'{assembly["telomeres"]}: invalid or duplicate telomere row for {chrom}/{end}')
            ends[chrom, end] = status
    if assembly.get('rdna'):
        for chrom, start, end, kind in table(assembly['rdna'], 'chromosome start end type'):
            if chrom not in lengths or kind not in RDNA:
                raise ValueError(f'{assembly["rdna"]}: invalid marker for {chrom}/{kind}')
            start, end = interval(int(start), int(end), lengths[chrom], assembly['rdna'])
            sites.append((chrom, start, end, kind))
    return ends, sites


def paf_links(path, lengths, pairs, minimum):
    chosen = {(a, b): index for index, (_, a, b) in enumerate(pairs)}
    seen = set()
    for number, row in records(path):
        if len(row) < 12 or row[4] not in ('+', '-'):
            raise ValueError(f'{path}:{number}: invalid PAF record')
        b, a = row[0], row[5]  # Query is the second assembly; target is the first.
        if (a, b) not in chosen or any(tag.startswith('tp:A:') and tag != 'tp:A:P' for tag in row[12:]):
            continue
        if (int(row[6]), int(row[1])) != (lengths[0][a], lengths[1][b]):
            raise ValueError(f'{path}:{number}: PAF lengths do not match both FAIs')
        s1, e1 = interval(int(row[7]), int(row[8]), lengths[0][a], f'{path}:{number}')
        s2, e2 = interval(int(row[2]), int(row[3]), lengths[1][b], f'{path}:{number}')
        if min(e1 - s1, e2 - s2) >= minimum:
            item = (chosen[a, b], a, s1, e1, b, s2, e2, row[4])
            if item not in seen:
                seen.add(item)
                yield item


def color(index):
    return ACCENTS[index] if index < len(ACCENTS) else '#526678'


def lighten(hex_color, factor=0.43):
    return '#' + ''.join(f'{round(int(hex_color[i:i+2], 16) * (1 - factor) + 255 * factor):02x}' for i in (1, 3, 5))


def point(r, a):
    radians = math.radians(a)
    return f'{CX + r * math.sin(radians):.3f},{CY - r * math.cos(radians):.3f}'


def arc(r, start, end, sweep=1):
    return f'A{r},{r} 0 {int(abs(end - start) > 180)} {sweep} {point(r, end)}'


def ring(r1, r2, start, end, fill, **kwargs):
    if r1 >= r2 or start >= end:
        return ''
    return svg_path(f'M{point(r2, start)} {arc(r2, start, end)} L{point(r1, end)} {arc(r1, end, start, 0)} Z', fill, **kwargs)


def svg_path(d, fill, **kwargs):
    attrs = ' '.join(f'{key.replace("_", "-")}="{val}"' for key, val in kwargs.items())
    return f'<path d="{d}" fill="{fill}" {attrs}/>'


def text(x, y, value, size=16, anchor='start', style=''):
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="sans-serif" '
            f'font-size="{size}" fill="{INK}" {style}>{escape(str(value))}</text>')


def ribbon(a, b, c, d, shade, reversed_strand):
    radius = 168
    enter, leave = (c, d) if reversed_strand else (d, c)
    path = (f'M{point(radius, a)} {arc(radius, a, b)} '
            f'Q{CX},{CY} {point(radius, enter)} {arc(radius, enter, leave, int(reversed_strand))} '
            f'Q{CX},{CY} {point(radius, a)} Z')
    return svg_path(path, shade, fill_opacity='0.22')


def draw(names, pairs, lengths, values, annotations, links, window):
    size = sum(n for sizes in lengths for n in sizes.values())
    gap = min(1.2, 36 / (2 * len(pairs)))
    scale = (343 - gap * (2 * len(pairs) - 1)) / size
    cursor = 8.5
    extents = {}
    for index, (_, a, b) in enumerate(pairs):
        for assembly, chrom in enumerate((a, b)):
            end = cursor + lengths[assembly][chrom] * scale
            extents[assembly, chrom] = (cursor, end, color(index))
            cursor = end + gap

    def angle(assembly, chrom, bp):
        return extents[assembly, chrom][0] + bp * scale

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{SIZE}" height="{SIZE}" viewBox="0 0 {SIZE} {SIZE}">',
           '<title>Paired genome Circos overview</title>', '<rect width="1080" height="1080" fill="white"/>']
    for index, a, s1, e1, b, s2, e2, strand in links:
        svg.append(ribbon(angle(0, a, s1), angle(0, a, e1), angle(1, b, s2), angle(1, b, e2), color(index), strand == '-'))
    for index, (label, a, b) in enumerate(pairs):
        pair_start, pair_end = extents[0, a][0], extents[1, b][1]
        for assembly, chrom in enumerate((a, b)):
            start, end, shade = extents[assembly, chrom]
            if assembly == 0:
                svg.append(ring(364, 391, start, end, shade, stroke=INK, stroke_width='0.85'))
            else:
                svg.append(ring(364, 391, start, end, lighten(shade)))
                # Emit the dashes as geometry so SVG-to-PDF converters preserve them.
                dashes = []
                for radius in (364, 391):
                    steps = max(1, round((end - start) / 0.9))
                    step = (end - start) / steps
                    for j in range(steps):
                        lo = start + (j + .1) * step
                        hi = start + (j + .7) * step
                        dashes.append(f'M{point(radius, lo)} {arc(radius, lo, hi)}')
                for edge in (start, end):
                    for r in (373, 383):
                        dashes.append(f'M{point(r, edge)} L{point(r + 5, edge)}')
                svg.append(svg_path(' '.join(dashes), 'none', stroke=INK, stroke_width='1.1'))
            ends, sites = annotations[assembly]
            # A narrow, separate strip prevents the terminal status from obscuring the sector shade.
            cap = min((end - start) / 5, 1.55)
            for side, lo, hi in (('left', start, start + cap), ('right', end - cap, end)):
                status = ends.get((chrom, side), 'not_assessed')
                svg.append(ring(364, 372, lo, hi, TEL[status]))
                if status == 'not_detected':
                    svg.append(svg_path(f'M{point(364, lo)} L{point(372, hi)}', 'none', stroke='white', stroke_width='1.1'))
            for track, (_, fill) in enumerate(TRACKS):
                radius = 331 - track * 32
                svg.append(ring(radius, radius + 25, start, end, '#F3F5F7'))
                bins = values[assembly][chrom][track]
                maximum = max(bins)
                if not maximum:
                    continue
                for i, count in enumerate(bins):
                    if count:
                        lo = angle(assembly, chrom, i * window)
                        hi = angle(assembly, chrom, min((i + 1) * window, lengths[assembly][chrom]))
                        svg.append(ring(radius, radius + 25 * count / maximum, lo, lo + 0.92 * (hi - lo), fill))
            svg.append(ring(177, 198, start, end, '#F3F5F7'))
            for seq, lo, hi, kind in sites:
                if seq == chrom:
                    # Minimum tick width is schematic and stays within the sector.
                    x, y = angle(assembly, chrom, lo), angle(assembly, chrom, hi)
                    mid = (x + y) / 2
                    low_radius = 177 if kind == '5S' else 188
                    svg.append(ring(low_radius, low_radius + 10, max(start, min(x, mid - .32)), min(end, max(y, mid + .32)), RDNA[kind]))

        middle = (pair_start + pair_end) / 2
        x, y = map(float, point(414, middle).split(','))
        rotation = middle if middle < 90 else middle - 180 if middle < 270 else middle - 360
        svg.append(text(x, y + 4, label, 15, 'middle', f'transform="rotate({rotation:.2f} {x} {y})"'))
    for letter, radius in zip('abcdefg', (378, 343, 311, 279, 247, 215, 187)):
        svg.append(text(CX, CY - radius + 5, letter, 14, 'middle', 'font-weight="bold"'))
    svg.append(svg_path('M132,886 H948', 'none', stroke='#CFD6DE', stroke_width='1'))
    for x, assembly, name in ((160, 0, names[0]), (550, 1, names[1])):
        shade = lighten('#406687') if assembly else '#406687'
        border = 'none' if assembly else INK
        svg.append(f'<rect x="{x}" y="906" width="25" height="15" fill="{shade}" stroke="{border}"/>')
        if assembly:
            dashes = [f'M{x+j},{y} h3' for y in (906, 921) for j in range(0, 25, 5)]
            dashes += [f'M{xx},{y} v3' for xx in (x, x + 25) for y in (907, 912, 917)]
            svg.append(svg_path(' '.join(dashes), 'none', stroke=INK, stroke_width='1.1'))
        svg.append(text(x + 36, 919, ('First: ' if assembly == 0 else 'Second: ') + name, 16))
    for index, (label, fill) in enumerate(TRACKS):
        column, row = index % 3, index // 3
        x = (160, 407, 665)[column]
        y = 947 + 25 * row
        svg.append(f'<rect x="{x}" y="{y-12}" width="14" height="12" fill="{fill}"/>')
        svg.append(text(x + 22, y, f'{chr(98 + index)}  {label}', 14))
    svg.append(text(687, 972, 'g  5S rDNA (blue); 18S rDNA (orange)', 14))
    svg.append(text(160, 1006, 'a  End caps: green, detected; red/slash, not detected; gray, not assessed.', 14))
    if links:
        svg.append(text(160, 1034, 'Center: alignments between paired chromosomes.', 14))
    svg.append('</svg>')
    return '\n'.join(svg) + '\n'


def caption(names, count, window, links, minimum):
    first, second = names
    window_label = f'{window // 1_000_000:,} Mb' if window % 1_000_000 == 0 else f'{window:,} bp'
    return (
        f'**Circular genome comparison of {first} and {second}.** Starting at the top gap and moving '
        f'clockwise, each of the {count} chromosome pairs shows {first} first (darker fill, solid outline) '
        f'and {second} second (lighter fill, dashed outline). From outside to inside, (a) shows chromosome '
        'ideograms and telomeric-repeat detection (green, detected; red with a white slash, not detected; '
        'gray, not assessed). Terminal caps are schematic. Tracks (b–f) show gene counts, LTR/Gypsy, '
        'LTR/Copia, DNA transposons, and LINE and unclassified repeats, respectively. Track (g) marks '
        f'provided 5S (blue) and 18S (orange) rDNA intervals; unmarked locations do not establish absence. '
        f'Gene models are counted by midpoint, and repeat tracks show non-overlapping annotated bases '
        f'within each repeat group per {window_label} window. Bar heights are normalized independently '
        'to the maximum for each track within each chromosome sector, and are not absolute comparisons '
        'between assemblies.' +
        (f' Central ribbons represent input PAF alignments between paired chromosomes with spans of at '
         f'least {minimum:,} bp on both assemblies; they do not represent structural-variant calls.' if links else '') + '\n'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path, help='JSON input manifest; paths relative to the manifest')
    parser.add_argument('output', type=Path, help='new output directory; existing directories are not overwritten')
    parser.add_argument('--export', action='store_true', help='also produce vector PDF and 300-dpi PNG (requires PyMuPDF)')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError(f'{args.output} already exists; use a new output directory')
    config = json.loads(args.manifest.read_text(encoding='utf-8'))
    source = config['assemblies']
    names = [str(item['label']).strip() for item in source]
    if len(names) != 2 or len(set(names)) != 2 or any(not name or len(name) > 32 for name in names):
        raise ValueError('provide two distinct assembly labels, each 1–32 characters long')
    window = int(config.get('window_bp', 1_000_000))
    minimum = int(config.get('min_link_bp', 10_000))
    if min(window, minimum) < 1:
        raise ValueError('window_bp and min_link_bp must be positive')
    root = args.manifest.resolve().parent
    assemblies = [{key: root / item[key] for key in ('fai', 'genes', 'repeats', 'telomeres', 'rdna') if item.get(key)} for item in source]
    pairs = list(table(root / config['pairs'], 'label first second'))
    if not pairs or any(len({row[i] for row in pairs}) != len(pairs) for i in range(3)):
        raise ValueError('chromosome pairs must be non-empty, with unique labels and unique names on each side')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]{1,20}', row[0]) for row in pairs):
        raise ValueError('pair labels must be 1–20 characters: ASCII letters, digits, underscore, dot or hyphen')
    sizes = [fai(item['fai']) for item in assemblies]
    for index, (_, a, b) in enumerate(pairs, 1):
        if a not in sizes[0] or b not in sizes[1]:
            raise ValueError(f'pair {index}: {a}/{b} not found in their respective FAIs')
    lengths = [{row[i + 1]: sizes[i][row[i + 1]] for row in pairs} for i in range(2)]
    values = [densities(item, lens, window) for item, lens in zip(assemblies, lengths)]
    marks = [markers(item, lens) for item, lens in zip(assemblies, lengths)]
    links = list(paf_links(root / config['paf'], lengths, pairs, minimum)) if config.get('paf') else []
    if config.get('paf') and not links:
        print('Warning: no PAF links passed the pair and length filters.', file=sys.stderr)
    svg = draw(names, pairs, lengths, values, marks, links, window)
    if args.export:
        try:
            import pymupdf
        except ImportError as exc:
            raise ValueError('--export requires PyMuPDF; omit it for standard-library SVG output') from exc
        with pymupdf.open(stream=svg.encode('utf-8'), filetype='svg') as document:
            pdf = document.convert_to_pdf()
            png = document[0].get_pixmap(dpi=300).tobytes('png')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'circos.svg').write_text(svg, encoding='utf-8')
    (args.output / 'legend.en.md').write_text(caption(names, len(pairs), window, links, minimum), encoding='utf-8')
    if args.export:
        (args.output / 'circos.pdf').write_bytes(pdf)
        (args.output / 'circos.png').write_bytes(png)
    with (args.output / 'window_values.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['assembly', 'chromosome', 'start', 'end', 'gene_count', 'gypsy_bp', 'copia_bp', 'dna_te_bp', 'line_unclassified_bp'])
        for name, lens, rows in zip(names, lengths, values):
            for chrom, tracks in rows.items():
                for i, counts in enumerate(zip(*tracks)):
                    writer.writerow([name, chrom, i * window, min((i + 1) * window, lens[chrom]), *counts])
    print(f'{len(pairs)} chromosome pairs, {len(links)} alignment ribbons: {args.output}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as error:
        sys.exit(f'Error: {error}')
