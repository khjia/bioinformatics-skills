#!/usr/bin/env python3
"""Validate a standardized per-library RNA-seq QC TSV, without modifying files."""

import argparse
import csv
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path
from statistics import mean, median
import sys

COUNT = ('raw_reads', 'raw_bases', 'clean_reads', 'clean_bases', 'mapped_reads', 'mapping_input_reads')
RATE = ('raw_q20_pct', 'raw_q30_pct', 'raw_gc_pct', 'clean_q20_pct', 'clean_q30_pct', 'clean_gc_pct', 'mapping_pct')
FIELDS = {'sample_id', *COUNT, *RATE}
MISSING = {'', 'na', 'n/a'}


def validate(path, expected):
    ids, failures, warnings, absent, valid = set(), [], [], {}, []
    with path.open(encoding='utf-8-sig', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t', strict=True)
        keys = reader.fieldnames or []
        if len(keys) != len(set(keys)):
            raise ValueError('Duplicate TSV column names')
        if FIELDS - set(keys):
            raise ValueError('Missing TSV columns: ' + ', '.join(sorted(FIELDS - set(keys))))
        for row in reader:
            where = f'line {reader.line_num}'
            n_failures = len(failures)
            if None in row or any(row[k] is None for k in keys):
                failures.append(f'{where}: row length differs from header')
                continue
            sid = row['sample_id'].strip()
            if not sid or sid.casefold() in MISSING:
                failures.append(f'{where}: missing sample_id')
                continue
            if sid in ids:
                failures.append(f'{where}: duplicate sample_id {sid!r}')
            ids.add(sid)
            data = {}
            for key in (*COUNT, *RATE):
                value = row[key].strip()
                if value.casefold() in MISSING:
                    absent[key] = absent.get(key, 0) + 1
                    data[key] = None
                    continue
                try:
                    n = Decimal(value)
                    if not n.is_finite() or n < 0:
                        raise ValueError()
                    if key in COUNT:
                        if n != n.to_integral_value():
                            raise ValueError()
                        n = int(n)
                    elif n > 100:
                        raise ValueError()
                    data[key] = n
                except (ValueError, InvalidOperation, OverflowError):
                    failures.append(f'{where}: invalid value for {key}: {value!r}')
            if len(failures) != n_failures:
                continue

            def check_le(left, right):
                if data[left] is not None and data[right] is not None and data[left] > data[right]:
                    failures.append(f'{where}: {left} exceeds {right}; check units and provenance')

            for field in ('reads', 'bases'):
                check_le(f'clean_{field}', f'raw_{field}')
            for stage in ('raw', 'clean'):
                check_le(f'{stage}_q30_pct', f'{stage}_q20_pct')
            check_le('mapped_reads', 'mapping_input_reads')
            num, den, rate = (data[k] for k in ('mapped_reads', 'mapping_input_reads', 'mapping_pct'))
            clean = data['clean_reads']
            if den is not None and clean is not None and den != clean:
                warnings.append(f'{where}: mapper input differs from clean reads; confirm the denominator')
            if den == 0 and rate is not None:
                failures.append(f'{where}: mapping rate needs a nonzero denominator; use NA when undefined')
            if den and num is not None and rate is not None:
                # Match the precision as shown in the TSV, e.g. 96.47 has a half-step of 0.005.
                decimals = min(12, max(0, -rate.as_tuple().exponent))
                unit = Decimal(1).scaleb(-decimals)
                computed = Decimal(100) * num / den
                # Allow either rounding direction on an exact halfway point.
                if abs(computed - rate) > unit / 2 + Decimal('1e-20'):
                    correct = computed.quantize(unit, rounding=ROUND_HALF_UP)
                    failures.append(f'{where}: mapping_pct disagrees with integer counts (computed {correct})')
            valid.append(data)
    if not ids:
        failures.append('No sample rows')
    if expected is not None:
        for diff, prefix in [(expected - ids, 'Missing expected samples'), (ids - expected, 'Unexpected samples')]:
            if diff:
                failures.append(prefix + ': ' + ', '.join(sorted(diff)))
    return valid, failures, warnings, absent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('table', type=Path, help='UTF-8 TSV; NA for missing values')
    parser.add_argument('--expected-samples', type=Path, help='Optional UTF-8 list of sample IDs, one per line without a header')
    args = parser.parse_args()
    try:
        expected = None
        if args.expected_samples:
            listed = [s.strip() for s in args.expected_samples.read_text(encoding='utf-8-sig').splitlines() if s.strip()]
            expected = set(listed)
            if not expected or len(listed) != len(expected):
                raise ValueError('Expected-samples list is empty or contains duplicate IDs')
        rows, errors, warnings, missing = validate(args.table, expected)
    except (OSError, UnicodeError, ValueError, csv.Error) as err:
        print(f'ERROR: {err}', file=sys.stderr)
        return 2
    for warning in warnings:
        print(f'WARN: {warning}', file=sys.stderr)
    for err in errors:
        print(f'ERROR: {err}', file=sys.stderr)
    if errors:
        return 2
    print(f'PASS: {len(rows)} sample rows; checked TSV structure and internal consistency.')
    if missing:
        print('Missing values by column: ' + ', '.join(f'{k}={v}' for k, v in sorted(missing.items())))
    rates = [r['mapping_pct'] for r in rows if r['mapping_pct'] is not None]
    if rates:
        print(f'Per-sample mapping (%; {len(rates)}/{len(rows)}): mean {mean(rates):.2f}, median {median(rates):.2f}, range {min(rates):.2f}–{max(rates):.2f}.')
    complete = [r for r in rows if r['mapped_reads'] is not None and r['mapping_input_reads'] is not None]
    if len(complete) == len(rows):
        num = sum(r['mapped_reads'] for r in complete)
        den = sum(r['mapping_input_reads'] for r in complete)
        if den:
            print(f'All-sample read-weighted mapping: {Decimal(100) * num / den:.2f}% ({num:,}/{den:,}).')
        else:
            print('All-sample read-weighted mapping not reported: total mapper input is zero.')
    else:
        print(f'All-sample read-weighted mapping not reported: exact counts for {len(complete)}/{len(rows)} samples.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
