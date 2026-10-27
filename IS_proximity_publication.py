#!/usr/bin/env python3
"""Reproducible IS–determinant proximity analysis (AMR, VF, METRES).

Requires Python >=3.9, pandas, numpy, openpyxl (for XLSX input).

Input conventions (1-based, inclusive coordinates):
  AMR XLSX: 'Contig id', 'Gene symbol', 'Start', 'Stop'.
  VF XLSX: 'ID', 'Unnamed: 7', 'Unnamed: 8', 'Unnamed: 14';
      gene label extracted from the first parenthesized name in description.
  METRES CSV/XLSX: headerless BLAST-like table; columns 0,1,7,8,14;
      gene label is text before first underscore in column 1.
  IS CSV: headerless table; columns 0,1,7,8; IS name is before first '|'.
  Lengths CSV (optional): columns 'Contig id', 'Length_bp'.

Distances are base-pair coordinate gaps between inclusive intervals (adjacent
intervals have distance 1); overlapping intervals have distance 0. Circular
mode requires lengths and uses the shortest distance around the circular origin.

IS–target–IS modules are *candidate* architectures, NOT validated transposons.
They comprise target genes fully enclosed by two consecutive, nonoverlapping
IS annotations, with no intervening IS and within a user-defined span.
Other unannotated genes may lie inside. Circular-origin-spanning modules are
not called, even in circular proximity mode.

Permutation null: each target interval is independently relocated uniformly
within its own plasmid, preserving its length and the observed IS coordinates.
This is a positional null model, not a phylogenetically independent test.
Use length metadata for all analyzed plasmids and interpret results alongside
plasmid size, IS density, gene density and host/Inc-group composition.
"""
from __future__ import annotations

import argparse
import json
import logging
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd

LOG = logging.getLogger('is_proximity')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--analysis', choices=['AMR', 'VF', 'METRES'], required=True)
    p.add_argument('--target', type=Path, required=True)
    p.add_argument('--is-file', type=Path, required=True)
    p.add_argument('--replicon', required=True, help='Output label (e.g. IncFIB)')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--lengths', type=Path, help='CSV: Contig id,Length_bp')
    p.add_argument('--circular', action='store_true', help='Circular distances; requires --lengths')
    p.add_argument('--thresholds', type=int, nargs='+', default=[1000, 5000, 10000])
    p.add_argument('--window', type=int, default=10000)
    p.add_argument('--module-max-span', type=int, default=30000)
    p.add_argument('--permutations', type=int, default=0, help='Within-plasmid positional null; requires --lengths')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--is-name-pattern', default=r'^IS', help='Regex matching classical IS labels')
    p.add_argument('--overwrite', action='store_true')
    a = p.parse_args()
    if any(t < 0 for t in a.thresholds) or a.window < 0 or a.module_max_span < 1 or a.permutations < 0:
        p.error('Distances/spans and permutation count must be nonnegative; module span must be positive')
    if (a.circular or a.permutations) and a.lengths is None:
        p.error('--lengths is required for --circular and --permutations')
    if len(set(a.thresholds)) != len(a.thresholds):
        p.error('Thresholds must be unique')
    return a


def load_table(path, *, header=0):
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() in ('.xlsx', '.xls'):
        return pd.read_excel(path, header=header)
    return pd.read_csv(path, header=header, low_memory=False)


def require_columns(df, columns, source):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f'{source}: missing columns {missing}. Available: {list(df.columns)}')


def load_targets(path, analysis):
    if analysis == 'AMR':
        raw = load_table(path)
        cols = ['Contig id', 'Gene symbol', 'Start', 'Stop']
        require_columns(raw, cols, str(path))
        df = raw[cols].copy()
        df.columns = ['Contig id', 'Gene', 'coord1', 'coord2']
    elif analysis == 'VF':
        raw = load_table(path)
        cols = ['ID', 'Unnamed: 7', 'Unnamed: 8', 'Unnamed: 14']
        require_columns(raw, cols, str(path))
        df = raw[cols].copy()
        df.columns = ['Contig id', 'coord1', 'coord2', 'Description']
        df['Gene'] = df['Description'].astype('string').str.extract(r'^\S+\s+\(([^)]+)\)', expand=False)
    else:
        raw = load_table(path, header=None)
        cols = [0, 1, 7, 8, 14]
        require_columns(raw, cols, str(path))
        df = raw[cols].copy()
        df.columns = ['Contig id', 'Hit', 'coord1', 'coord2', 'Description']
        df['Gene'] = df['Hit'].astype('string').str.split('_', n=1).str[0]
    return clean_intervals(df, 'Gene', 'targets')


def load_is(path, pattern):
    raw = load_table(path, header=None)
    require_columns(raw, [0, 1, 7, 8], str(path))
    df = raw[[0, 1, 7, 8]].copy()
    df.columns = ['Contig id', 'Full_IS_name', 'coord1', 'coord2']
    df['IS'] = df['Full_IS_name'].astype('string').str.split('|', n=1, regex=False).str[0].str.strip()
    before = len(df)
    df = df.loc[df['IS'].str.contains(pattern, regex=True, na=False)].copy()
    LOG.info('IS name filter retained %s / %s records', len(df), before)
    return clean_intervals(df, 'IS', 'IS')


def clean_intervals(df, label, source):
    before = len(df)
    for c in ('coord1', 'coord2'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    valid = (df['coord1'].notna() & df['coord2'].notna()
             & np.isfinite(df['coord1']) & np.isfinite(df['coord2'])
             & (df['coord1'] >= 1) & (df['coord2'] >= 1)
             & (df['coord1'] % 1 == 0) & (df['coord2'] % 1 == 0)
             & df['Contig id'].notna() & df[label].notna())
    df = df.loc[valid].copy()
    df['Contig id'] = df['Contig id'].astype(str).str.strip()
    df[label] = df[label].astype(str).str.strip()
    df = df.loc[(df['Contig id'] != '') & (df[label] != '')].copy()
    df['start'] = df[['coord1', 'coord2']].min(axis=1).astype('int64')
    df['end'] = df[['coord1', 'coord2']].max(axis=1).astype('int64')
    df = df.drop_duplicates(subset=['Contig id', 'start', 'end', label])
    df = df.sort_values(['Contig id', 'start', 'end', label]).reset_index(drop=True)
    df.insert(0, 'Feature_ID', [f'{source}_{i+1:08d}' for i in range(len(df))])
    LOG.info('%s: retained %d unique valid records; %d excluded/duplicates', source, len(df), before-len(df))
    return df


def read_lengths(path):
    df = load_table(path)
    require_columns(df, ['Contig id', 'Length_bp'], str(path))
    if df['Contig id'].isna().any():
        raise ValueError('Missing contig ID in length file')
    df['Contig id'] = df['Contig id'].astype(str).str.strip()
    df['Length_bp'] = pd.to_numeric(df['Length_bp'], errors='raise')
    if df['Contig id'].duplicated().any() or (df['Length_bp'] < 1).any() or (df['Length_bp'] % 1 != 0).any():
        raise ValueError('Lengths must be unique positive integers per contig')
    return dict(zip(df['Contig id'], df['Length_bp'].astype(int)))


def validate_lengths(target, is_df, lengths, required):
    if lengths is None:
        return
    for name, df in [('target', target), ('IS', is_df)]:
        missing = sorted(set(df['Contig id']) - set(lengths))
        if missing and required:
            raise ValueError(f'Missing length metadata for {name} contigs: {missing[:10]}')
        check = df.loc[df['Contig id'].isin(lengths)].copy()
        if not check.empty and (check['end'] > check['Contig id'].map(lengths)).any():
            raise ValueError(f'{name} coordinates exceed sequence length; check coordinate conventions')


def interval_distances(a, b, c, d, length=None):
    """Broadcast interval distances; inputs are start/end arrays, 1-based inclusive."""
    linear = np.maximum(np.maximum(c - b, a - d), 0)
    if length is None:
        return linear
    # Intervals are non-wrapping; shortest distance can cross origin.
    via_origin_1 = a + length - d
    via_origin_2 = c + length - b
    return np.minimum(linear, np.minimum(via_origin_1, via_origin_2))


def analyze_proximity(target, is_df, lengths, circular, thresholds, window):
    is_groups = {k: g.reset_index(drop=True) for k, g in is_df.groupby('Contig id', sort=False)}
    nearest_rows, window_rows = [], []
    for contig, tg in target.groupby('Contig id', sort=False):
        ig = is_groups.get(contig)
        for t in tg.itertuples(index=False):
            base = {'Contig id': contig, 'Target_ID': t.Feature_ID, 'Gene': t.Gene,
                    'Gene_start': t.start, 'Gene_end': t.end}
            if ig is None or ig.empty:
                nearest_rows.append({**base, 'Nearest_IS_ID': pd.NA, 'Nearest_IS': pd.NA,
                                     'Distance_bp': np.nan})
                continue
            dist = interval_distances(t.start, t.end, ig['start'].to_numpy(),
                                      ig['end'].to_numpy(), lengths[contig] if circular else None)
            j = int(np.argmin(dist))
            nearest_rows.append({**base, 'Nearest_IS_ID': ig.iloc[j]['Feature_ID'],
                                 'Nearest_IS': ig.iloc[j]['IS'], 'Distance_bp': int(dist[j])})
            for k in np.flatnonzero(dist <= window):
                s = ig.iloc[int(k)]
                window_rows.append({**base, 'IS_ID': s['Feature_ID'], 'IS': s['IS'],
                                    'IS_start': int(s['start']), 'IS_end': int(s['end']),
                                    'Distance_bp': int(dist[k])})
    nearest_cols = ['Contig id', 'Target_ID', 'Gene', 'Gene_start', 'Gene_end',
                    'Nearest_IS_ID', 'Nearest_IS', 'Distance_bp']
    window_cols = ['Contig id', 'Target_ID', 'Gene', 'Gene_start', 'Gene_end',
                   'IS_ID', 'IS', 'IS_start', 'IS_end', 'Distance_bp']
    nearest = pd.DataFrame(nearest_rows, columns=nearest_cols)
    window_df = pd.DataFrame(window_rows, columns=window_cols)
    for threshold in thresholds:
        nearest[f'Within_{threshold}bp'] = nearest['Distance_bp'].le(threshold).fillna(False)
    return nearest, window_df


def canonical(left, genes, right):
    fwd = (left, tuple(genes), right)
    rev = (right, tuple(reversed(genes)), left)
    return min(fwd, rev)


def call_modules(target, is_df, max_span):
    records = []
    tg = {k: g for k, g in target.groupby('Contig id', sort=False)}
    for contig, ig in is_df.groupby('Contig id', sort=False):
        ig = ig.sort_values(['start', 'end', 'IS']).reset_index(drop=True)
        if contig not in tg or len(ig) < 2:
            continue
        genes = tg[contig]
        for i in range(len(ig)-1):
            left, right = ig.iloc[i], ig.iloc[i+1]
            # Avoid overlapping/nested IS calls and ambiguously ordered features.
            if left['end'] >= right['start']:
                continue
            span = int(right['end'] - left['start'] + 1)
            if span > max_span:
                continue
            contained = genes.loc[(genes['start'] > left['end']) & (genes['end'] < right['start'])]
            if contained.empty:
                continue
            contained = contained.sort_values(['start', 'end', 'Gene'])
            names = tuple(contained['Gene'].astype(str))
            a, b, c = canonical(str(left['IS']), names, str(right['IS']))
            records.append({'Contig id': contig, 'Left_IS_ID': left['Feature_ID'],
                            'Right_IS_ID': right['Feature_ID'], 'Span_bp': span,
                            'Left_IS': left['IS'], 'Target_gene_order': ' | '.join(names),
                            'Right_IS': right['IS'], 'Target_IDs': ' | '.join(contained['Feature_ID']),
                            'Canonical_Left_IS': a, 'Canonical_gene_order': ' | '.join(b),
                            'Canonical_Right_IS': c,
                            'Canonical_pattern': ' — '.join((a, *b, c)),
                            'Same_IS_label': a == c})
    cols = ['Contig id', 'Left_IS_ID', 'Right_IS_ID', 'Span_bp', 'Left_IS',
            'Target_gene_order', 'Right_IS', 'Target_IDs', 'Canonical_Left_IS',
            'Canonical_gene_order', 'Canonical_Right_IS', 'Canonical_pattern', 'Same_IS_label']
    modules = pd.DataFrame(records, columns=cols)
    keys = ['Canonical_Left_IS', 'Canonical_gene_order', 'Canonical_Right_IS', 'Canonical_pattern']
    summary_cols = keys + ['n_plasmids', 'occurrences', 'median_span_bp', 'min_span_bp', 'max_span_bp']
    if modules.empty:
        return modules, pd.DataFrame(columns=summary_cols)
    summary = (modules.groupby(keys, dropna=False).agg(
        n_plasmids=('Contig id', 'nunique'), occurrences=('Contig id', 'size'),
        median_span_bp=('Span_bp', 'median'), min_span_bp=('Span_bp', 'min'),
        max_span_bp=('Span_bp', 'max')).reset_index()
        .sort_values(['n_plasmids', 'occurrences'], ascending=False))
    return modules, summary


def permutation_test(target, is_df, lengths, circular, thresholds, n_perm, seed):
    """Empirical one-sided enrichment p-values for proximity fractions."""
    if not n_perm or target.empty:
        return pd.DataFrame(columns=['Threshold_bp', 'Observed_fraction_all',
                                     'Null_mean', 'Null_sd', 'Empirical_p_enrichment',
                                     'Permutations', 'Seed'])
    rng = np.random.default_rng(seed)
    groups = []
    for contig, tg in target.groupby('Contig id', sort=False):
        length = lengths[contig]
        sizes = (tg['end'] - tg['start'] + 1).to_numpy(dtype=int)
        ig = is_df.loc[is_df['Contig id'] == contig]
        groups.append((length, sizes, ig['start'].to_numpy(dtype=int), ig['end'].to_numpy(dtype=int)))
    obs = np.zeros(len(thresholds), dtype=int)
    for contig, tg in target.groupby('Contig id', sort=False):
        ig = is_df.loc[is_df['Contig id'] == contig]
        if ig.empty:
            continue
        for t in tg.itertuples(index=False):
            ds = interval_distances(t.start, t.end, ig['start'].to_numpy(),
                                    ig['end'].to_numpy(), lengths[contig] if circular else None)
            obs += np.array([ds.min() <= th for th in thresholds], dtype=int)
    observed = obs / len(target)
    null = np.zeros((n_perm, len(thresholds)), dtype=float)
    for p in range(n_perm):
        count = np.zeros(len(thresholds), dtype=int)
        for length, sizes, starts, ends in groups:
            if len(starts) == 0:
                continue
            # Uniform starts; gene length preserved. No wrapped genes are generated.
            random_starts = rng.integers(1, length - sizes + 2)
            for start, size in zip(random_starts, sizes):
                ds = interval_distances(start, start + size - 1, starts, ends,
                                        length if circular else None)
                count += np.array([ds.min() <= th for th in thresholds], dtype=int)
        null[p] = count / len(target)
    return pd.DataFrame({
        'Threshold_bp': thresholds, 'Observed_fraction_all': observed,
        'Null_mean': null.mean(axis=0), 'Null_sd': null.std(axis=0, ddof=1) if n_perm > 1 else np.nan,
        'Empirical_p_enrichment': (1 + (null >= observed).sum(axis=0)) / (n_perm + 1),
        'Permutations': n_perm, 'Seed': seed})


def summarize(target, is_df, nearest, modules, families, thresholds):
    with_is = nearest.loc[nearest['Distance_bp'].notna()]
    records = [
        ('Determinant occurrences (all)', len(target)),
        ('Determinant-positive plasmids', target['Contig id'].nunique()),
        ('Classical IS occurrences', len(is_df)),
        ('IS-positive plasmids', is_df['Contig id'].nunique()),
        ('Determinants on IS-positive plasmids', len(with_is)),
        ('Determinants on IS-negative plasmids', len(nearest)-len(with_is)),
        ('Median nearest-IS distance (bp; IS-positive only)', with_is['Distance_bp'].median()),
        ('Candidate module occurrences', len(modules)),
        ('Candidate structural families', len(families)),
        ('Plasmids with candidate modules', modules['Contig id'].nunique()),
    ]
    for th in thresholds:
        n = int(nearest['Distance_bp'].le(th).sum())
        records.extend([
            (f'Determinants within {th} bp of IS', n),
            (f'Percent within {th} bp (all determinants)', 100*n/len(nearest) if len(nearest) else np.nan),
            (f'Percent within {th} bp (IS-positive plasmids only)', 100*n/len(with_is) if len(with_is) else np.nan),
        ])
    return pd.DataFrame(records, columns=['Metric', 'Value'])


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    if args.output.exists() and any(args.output.iterdir()) and not args.overwrite:
        raise FileExistsError(f'Output directory is not empty: {args.output}; use --overwrite')
    args.output.mkdir(parents=True, exist_ok=True)
    target = load_targets(args.target, args.analysis)
    is_df = load_is(args.is_file, args.is_name_pattern)
    lengths = read_lengths(args.lengths) if args.lengths else None
    validate_lengths(target, is_df, lengths, args.circular or bool(args.permutations))
    nearest, within = analyze_proximity(target, is_df, lengths, args.circular, args.thresholds, args.window)
    with_is = nearest.dropna(subset=['Distance_bp'])
    if with_is.empty:
        pairs = pd.DataFrame(columns=['Gene', 'Nearest_IS', 'occurrences', 'plasmids', 'median_distance_bp'])
        exact = pd.DataFrame(columns=['Gene', 'Nearest_IS', 'Distance_bp', 'occurrences', 'plasmids'])
    else:
        pairs = (with_is.groupby(['Gene', 'Nearest_IS']).agg(
            occurrences=('Target_ID', 'size'), plasmids=('Contig id', 'nunique'),
            median_distance_bp=('Distance_bp', 'median')).reset_index()
            .sort_values(['plasmids', 'occurrences'], ascending=False))
        exact = (with_is.groupby(['Gene', 'Nearest_IS', 'Distance_bp']).agg(
            occurrences=('Target_ID', 'size'), plasmids=('Contig id', 'nunique')).reset_index()
            .sort_values(['plasmids', 'occurrences'], ascending=False))
    modules, families = call_modules(target, is_df, args.module_max_span)
    same = families.loc[families['Canonical_Left_IS'] == families['Canonical_Right_IS']].copy()
    summary = summarize(target, is_df, nearest, modules, families, args.thresholds)
    perm = permutation_test(target, is_df, lengths, args.circular, args.thresholds,
                            args.permutations, args.seed)
    prefix = re.sub(r'[^A-Za-z0-9_.-]+', '_', f'{args.replicon}_{args.analysis}')
    outputs = {
        'global_summary': summary, 'nearest_IS': nearest, 'IS_within_window': within,
        'IS_pairs': pairs, 'IS_recurrent_distances': exact,
        'candidate_modules_all': modules, 'candidate_module_families': families,
        'candidate_same_IS_families': same, 'permutation_null': perm,
    }
    for suffix, df in outputs.items():
        df.to_csv(args.output / f'{prefix}_{suffix}.csv', index=False)
    metadata = {
        'analysis': args.analysis, 'replicon': args.replicon,
        'target_file': str(args.target.resolve()), 'IS_file': str(args.is_file.resolve()),
        'lengths_file': str(args.lengths.resolve()) if args.lengths else None,
        'circular_distance': args.circular, 'circular_module_calls': False,
        'thresholds_bp': args.thresholds, 'IS_window_bp': args.window,
        'module_max_span_bp': args.module_max_span, 'permutations': args.permutations,
        'random_seed': args.seed, 'IS_name_regex': args.is_name_pattern,
        'coordinate_system': '1-based inclusive; interval gap; overlap=0',
        'python_version': platform.python_version(), 'pandas_version': pd.__version__,
        'numpy_version': np.__version__,
        'limitations': [
            'IS–gene proximity does not prove mobilization or causality.',
            'Candidate IS–gene–IS modules are not validated transposons.',
            'Structural patterns omit unannotated intervening features and orientation strands.',
            'Circular-origin-spanning modules are not detected.',
            'Annotation quality and duplicated/clonal plasmids affect counts and inference.',
            'Positional permutations do not control for gene deserts, functional regions, or phylogeny.',
        ],
    }
    (args.output / f'{prefix}_run_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    LOG.info('Finished %s: %s determinants, %s IS, %s candidate modules',
             prefix, len(target), len(is_df), len(modules))
    print(summary.to_string(index=False))
    print(f'\nResults: {args.output.resolve()}')


if __name__ == '__main__':
    main()
