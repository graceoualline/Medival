# Process MEDIVAL output: compress overlapping coordinates, apply size filter,
# then cluster by gap distance.
#
# Each final interval records all reference metadata from every contributing
# input row.  Within a row the original comma separator is preserved (e.g.
# two refs per pair: "refA,refB").  Rows merged by clustering are joined with
# "|" so every "|"-delimited token corresponds to one original input pair.
#
# Usage:
#   python3 size_cluster_filter.py <input.tsv> <output.tsv>
#                                  [--size-filter N] [--cluster N] [--test-set]

import sys
import pandas as pd
from build_database_index import *

# Columns present in first-divergence-filter output but absent from overlap_div output.
# These are blat alignment details irrelevant to region filtering and are dropped on load.
_BLAT_ONLY_COLS = {
    'match', 'mismatch', 'rep. match', "N's",
    'Q gap count', 'Q gap bases', 'T gap count', 'T gap bases',
    'strand', 'block count', 'blockSizes', 'qStarts', 'tStarts',
}

# Reference-side columns aggregated across merged rows; order defines output column order
_META_COLS = [
    'T name',
    'T size',
    'T start',
    'T end',
    'Percent Identity',
    'Reference Species',
    'Divergence Time',
    'ANI<95(if div=unk)',
    'Div bt Ref Species',
    'ANI<95 bt Ref Seqs',
]


# ---------------------------------------------------------------------------
# Interval helpers — each interval carries the list of original rows that
# contributed to it so metadata can be aggregated at output time.
# ---------------------------------------------------------------------------

def _copy(iv):
    return {'start': iv['start'], 'end': iv['end'], 'rows': list(iv['rows'])}


def region_compress(ivs):
    """Merge overlapping or directly adjacent intervals, collecting contributing rows."""
    if not ivs:
        return []
    ivs = sorted(ivs, key=lambda x: x['start'])
    merged = [_copy(ivs[0])]
    for iv in ivs[1:]:
        last = merged[-1]
        if iv['start'] <= last['end'] + 1:
            last['end'] = max(last['end'], iv['end'])
            last['rows'].extend(iv['rows'])
        else:
            merged.append(_copy(iv))
    return merged


def size_filter(ivs, min_size):
    """Drop intervals shorter than min_size bp."""
    return [iv for iv in ivs if (iv['end'] - iv['start']) >= min_size]


def cluster_gap(ivs, gap):
    """Merge intervals whose gap is <= gap bp, collecting contributing rows."""
    if not ivs or gap <= 0:
        return ivs
    ivs = sorted(ivs, key=lambda x: x['start'])
    merged = [_copy(ivs[0])]
    for iv in ivs[1:]:
        last = merged[-1]
        if iv['start'] - last['end'] <= gap:
            last['end'] = max(last['end'], iv['end'])
            last['rows'].extend(iv['rows'])
        else:
            merged.append(_copy(iv))
    return merged


def process(ivs, min_size, gap):
    """Full pipeline: compress → size_filter → cluster_gap."""
    return cluster_gap(size_filter(region_compress(ivs), min_size), gap)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_data(file_path):
    """
    Returns:
        groups:          {q_name: [interval, ...]}
        plasmid_regions: {q_name: (mge_start, mge_end)}  when Q name is 4-part
    """
    df = pd.read_csv(file_path, sep='\t', on_bad_lines='skip', dtype=str, comment='#')
    df.columns = df.columns.str.strip()
    df.drop(columns=[c for c in _BLAT_ONLY_COLS if c in df.columns], inplace=True)

    required = {'Q name', 'Q start', 'Q end'}
    missing = required - set(df.columns)
    if missing:
        print(f"Error: Missing columns: {missing}")
        sys.exit(1)

    groups = {}           # q_name → list of raw single-row intervals
    plasmid_regions = {}  # q_name → (mge_start, mge_end)

    for _, row in df.iterrows():
        q_name  = str(row['Q name'])
        q_start = int(row['Q start'])
        q_end   = int(row['Q end'])

        # Q name format: plasmid_id,plasmid_start,plasmid_end,host_id
        parts = q_name.split(',')
        if len(parts) == 4:
            plasmid_regions[q_name] = (int(parts[1]), int(parts[2]))

        groups.setdefault(q_name, []).append(
            {'start': q_start, 'end': q_end, 'rows': [row.to_dict()]}
        )

    return groups, plasmid_regions


# ---------------------------------------------------------------------------
# Format one output row from a final interval
# ---------------------------------------------------------------------------

def _join(rows, col):
    """|'-join field values across all contributing rows (one token per pair)."""
    return '|'.join(str(r.get(col, '')) for r in rows)


def _summary_path(output_path):
    """Derive summary file path by inserting '_summary' before the extension."""
    if '.' in output_path.split('/')[-1]:
        base, ext = output_path.rsplit('.', 1)
        return f"{base}_summary.{ext}"
    return output_path + '_summary'


_SUMMARY_COLS = ['Q name', 'Q size', 'Q start', 'Q end', 'Query Species',
                 'Num Regions', 'Num Unique Species', 'Avg Divergence Time']


def format_summary_row(q_name, iv, index=None):
    first = iv['rows'][0]

    leaves = set()
    for r in iv['rows']:
        for ref_id in str(r.get('T name', '')).split(','):
            ref_id = ref_id.strip()
            if not ref_id:
                continue
            leaf = lookup_tree_leaf_name(index, ref_id) if index is not None else None
            leaves.add(leaf if leaf is not None else ref_id)

    div_vals = []
    for r in iv['rows']:
        for v in str(r.get('Divergence Time', '')).split(','):
            try:
                div_vals.append(float(v.strip()))
            except ValueError:
                pass

    return {
        'Q name':              q_name,
        'Q size':              first.get('Q size', ''),
        'Q start':             iv['start'],
        'Q end':               iv['end'],
        'Query Species':       first.get('Query Species', ''),
        'Num Regions':         len(iv['rows']),
        'Num Unique Species':  len(leaves),
        'Avg Divergence Time': round(sum(div_vals) / len(div_vals), 4) if div_vals else '',
    }


def format_row(q_name, iv, available_cols, test_set, plasmid_regions):
    first = iv['rows'][0]
    out = {
        'Q name':        q_name,
        'Q size':        first.get('Q size', ''),
        'Q start':       iv['start'],
        'Q end':         iv['end'],
        'Query Species': first.get('Query Species', ''),
    }
    for col in _META_COLS:
        if col in available_cols:
            out[col] = _join(iv['rows'], col)
    if test_set:
        if q_name in plasmid_regions:
            mge_s, mge_e = plasmid_regions[q_name]
            bp_mge = max(0, min(iv['end'], mge_e) - max(iv['start'], mge_s))
            out['bp_in_mge']  = bp_mge
            out['bp_in_host'] = (iv['end'] - iv['start']) - bp_mge
        else:
            out['bp_in_mge']  = None
            out['bp_in_host'] = None
    return out


def size_filter_cluster(input_path, output_path, summary_path, min_size, gap,
                        test_set=False, index_path=None, config_header=''):
    print(f"Loading: {input_path}")
    index = load_hash_table(index_path)
    groups, plasmid_regions = load_data(input_path)
    print(f"  {len(groups)} sequence(s) | size >= {min_size} bp | gap <= {gap} bp")
    if test_set:
        print("  --test-set: bp_in_mge / bp_in_host columns added")

    available_cols = set(
        pd.read_csv(input_path, sep='\t', nrows=0, comment='#').columns.str.strip()
    )
    meta_present = [c for c in _META_COLS if c in available_cols]
    fixed_cols = ['Q name', 'Q size', 'Q start', 'Q end', 'Query Species']
    test_cols  = ['bp_in_mge', 'bp_in_host'] if test_set else []
    out_cols   = fixed_cols + meta_present + test_cols

    rows_out = []
    summary_rows = []
    for q_name, ivs in groups.items():
        for iv in process(ivs, min_size, gap):
            summary_row = format_summary_row(q_name, iv, index)
            if summary_row['Num Unique Species'] <= 1:
                continue
            rows_out.append(format_row(q_name, iv, available_cols, test_set, plasmid_regions))
            summary_rows.append(summary_row)

    if config_header:
        with open(output_path, 'w') as f:
            f.write(config_header)
        pd.DataFrame(rows_out, columns=out_cols).to_csv(output_path, sep='\t', index=False, mode='a')
    else:
        pd.DataFrame(rows_out, columns=out_cols).to_csv(output_path, sep='\t', index=False)
    print(f"  {len(rows_out)} interval(s) written to: {output_path}")

    if config_header:
        with open(summary_path, 'w') as f:
            f.write(config_header)
        pd.DataFrame(summary_rows, columns=_SUMMARY_COLS).to_csv(summary_path, sep='\t', index=False, mode='a')
    else:
        pd.DataFrame(summary_rows, columns=_SUMMARY_COLS).to_csv(summary_path, sep='\t', index=False)
    print(f"  summary written to: {summary_path}")


def load_data_from_chunks(chunk_dir, chunk_size):
    """
    Load all *_overlap_div.tsv files from chunk_dir, applying the chunk coordinate
    offset (chunk_index * chunk_size) to Q start / Q end in-place.
    Returns the same (groups, plasmid_regions) structure as load_data().
    """
    import glob
    groups = {}
    plasmid_regions = {}

    chunk_files = sorted(glob.glob(os.path.join(chunk_dir, "*_overlap_div.tsv")))
    print(f"  Found {len(chunk_files)} chunk overlap_div files")

    for filepath in chunk_files:
        fname = os.path.basename(filepath)
        try:
            idx = int(fname.split("_")[1])
        except (IndexError, ValueError):
            print(f"  Warning: could not parse chunk index from {fname}, skipping")
            continue
        offset = idx * chunk_size

        try:
            df = pd.read_csv(filepath, sep='\t', on_bad_lines='skip', dtype=str, comment='#')
        except Exception as e:
            print(f"  Warning: could not read {filepath}: {e}")
            continue
        df.columns = df.columns.str.strip()
        df.drop(columns=[c for c in _BLAT_ONLY_COLS if c in df.columns], inplace=True)
        if df.empty:
            continue

        if offset > 0:
            df['Q start'] = pd.to_numeric(df['Q start'], errors='coerce').fillna(0).astype(int) + offset
            df['Q end']   = pd.to_numeric(df['Q end'],   errors='coerce').fillna(0).astype(int) + offset
        else:
            df['Q start'] = pd.to_numeric(df['Q start'], errors='coerce').fillna(0).astype(int)
            df['Q end']   = pd.to_numeric(df['Q end'],   errors='coerce').fillna(0).astype(int)

        for row in df.to_dict('records'):
            q_name  = str(row['Q name'])
            q_start = int(row['Q start'])
            q_end   = int(row['Q end'])

            parts = q_name.split(',')
            if len(parts) == 4:
                plasmid_regions[q_name] = (int(parts[1]), int(parts[2]))

            groups.setdefault(q_name, []).append(
                {'start': q_start, 'end': q_end, 'rows': [row]}
            )

    return groups, plasmid_regions


def size_filter_cluster_from_chunks(chunk_dir, chunk_size, output_path, summary_path,
                                     min_size, gap, test_set=False, index_path=None,
                                     config_header=''):
    """
    Like size_filter_cluster() but reads directly from per-chunk overlap_div files
    instead of one large combined file, avoiding the 900 MB+ memory spike.
    """
    print(f"Loading chunks from: {chunk_dir}")
    index = load_hash_table(index_path)
    groups, plasmid_regions = load_data_from_chunks(chunk_dir, chunk_size)
    print(f"  {len(groups)} sequence(s) | size >= {min_size} bp | gap <= {gap} bp")

    # Infer available columns from the first chunk file
    import glob
    available_cols = set()
    for fp in sorted(glob.glob(os.path.join(chunk_dir, "*_overlap_div.tsv"))):
        try:
            available_cols = set(
                pd.read_csv(fp, sep='\t', nrows=0, comment='#').columns.str.strip()
            )
            break
        except Exception:
            continue

    meta_present = [c for c in _META_COLS if c in available_cols]
    fixed_cols = ['Q name', 'Q size', 'Q start', 'Q end', 'Query Species']
    test_cols  = ['bp_in_mge', 'bp_in_host'] if test_set else []
    out_cols   = fixed_cols + meta_present + test_cols

    rows_out = []
    summary_rows = []
    for q_name, ivs in groups.items():
        for iv in process(ivs, min_size, gap):
            summary_row = format_summary_row(q_name, iv, index)
            if summary_row['Num Unique Species'] <= 1:
                continue
            rows_out.append(format_row(q_name, iv, available_cols, test_set, plasmid_regions))
            summary_rows.append(summary_row)

    if config_header:
        with open(output_path, 'w') as f:
            f.write(config_header)
        pd.DataFrame(rows_out, columns=out_cols).to_csv(output_path, sep='\t', index=False, mode='a')
    else:
        pd.DataFrame(rows_out, columns=out_cols).to_csv(output_path, sep='\t', index=False)
    print(f"  {len(rows_out)} interval(s) written to: {output_path}")

    if config_header:
        with open(summary_path, 'w') as f:
            f.write(config_header)
        pd.DataFrame(summary_rows, columns=_SUMMARY_COLS).to_csv(summary_path, sep='\t', index=False, mode='a')
    else:
        pd.DataFrame(summary_rows, columns=_SUMMARY_COLS).to_csv(summary_path, sep='\t', index=False)
    print(f"  summary written to: {summary_path}")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 size_cluster_filter.py <input.tsv> <output.tsv> "
              "[--size-filter N] [--cluster N] [--test-set] [--db <medival_db>]")
        sys.exit(1)
    # testset is for internal benchmarking of known MGE regions
    input_path  = sys.argv[1]
    output_path = sys.argv[2]
    test_set    = '--test-set' in sys.argv

    min_size = 0
    if '--size-filter' in sys.argv:
        i = sys.argv.index('--size-filter')
        min_size = int(sys.argv[i + 1])

    gap = 0
    if '--cluster' in sys.argv:
        i = sys.argv.index('--cluster')
        gap = int(sys.argv[i + 1])

    index = None
    if '--db' in sys.argv:
        i = sys.argv.index('--db')
        index = f'{sys.argv[i + 1]}/medival_db_index.pkl'

    size_filter_cluster(input_path, output_path, _summary_path(output_path),
                        min_size, gap, test_set, index)
