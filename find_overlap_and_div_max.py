# Drop-in replacement for find_overlap_and_div that pairs overlapping intervals
# with the LARGEST mutual overlap first instead of leftmost-start-first.
#
# Why: greedy left-to-right pairing can consume a row for a small overlap and
# leave a larger, higher-quality overlap unformed.  Largest-first pairing
# ensures the most informative pairs are always preferred.
#
# Complexity: O(n * k * log(n*k)) per Q-name group, where k is the average
# number of concurrently overlapping rows (bounded by local read density,
# not n).  Groups are processed independently so rows from different query
# sequences never interact.
#
# Usage (standalone):
#   python3 find_overlap_and_div_max.py <first_div.tsv> <out.tsv> <tree_dir> <medival_db>
#
# Usage (pipeline): import find_overlap_and_div_max and call
#   find_overlap_and_div_max(rows, output_file, tree, triangle_dict, index)
#   with the same signature as find_overlap_and_div.

import heapq
import sys
import pickle
from collections import defaultdict
from filter_blat import *
from build_database_index import *
import csv

_qs    = 2   # index into care-tuple: Q start
_qe    = 3   # Q end
_rsp   = 10  # Reference Species
_tname = 4   # T name

_HEADER = (
    "Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\t"
    "Percent Identity\tQuery Species\tReference Species\tDivergence Time\t"
    "ANI<95(if div=unk)\tDiv bt Ref Species\tANI<95 bt Ref Seqs\n"
)

def crude_ani_overlap(s1, e1, s2, e2, len1, len2):
    overlap = overlaps(s1, e1, s2, e2)
    if overlap == None: return -1
    if 0 in [len1, len2]: return -1
    return (overlap/min(len1, len2))*100

def lookup_ani_triangle(id1, id2, triangle_dict, index, s1=None, e1=None, s2=None, e2=None, cutoff=500):
    """
    Returns True if ANI < 95% (keep), False if ANI >= 95% (discard).
    Presence in triangle_dict means ANI >= 95.
    Falls back to crude_ani_overlap when either sequence is < 500 bp.
    """
    len1 = lookup_length(index, id1)
    len2 = lookup_length(index, id2)
    if len1 <= cutoff or len2 <= cutoff:
        ani = crude_ani_overlap(s1, e1, s2, e2, len1, len2)
        return isinstance(ani, (int, float)) and ani < 95
        
    key = (min(id1, id2), max(id1, id2))
    return key not in triangle_dict  # True = ANI < 95, False = ANI >= 95

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None

#if it finds regions that overlap and map to the same species, it will combine them
def compress(input_file):
    #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7

    groups = defaultdict(set)

    with open(input_file, 'r') as f:
        last_pos = f.tell()
        line = f.readline()
        while line and line.startswith('#'):
            last_pos = f.tell()
            line = f.readline()

        # Go back to the start of the first non-comment line
        f.seek(last_pos)
        tsv = csv.reader(f, delimiter="\t")
        for line in tsv:
            if line[1].isnumeric():
                care = line[9:17] + line[21:] #q name to ref Name
                groups[care[0]].add(tuple(care))

    compressed = []
    for row_set in groups.values():
        rows = sorted(list(row_set), key=lambda x: (int(x[qs]), int(x[qe])))
        if len(rows) <= 1:
            continue

        #Initialize first merged interval
        i = 0

        s_cur, e_cur, species_cur, row_cur = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], list(rows[i])

        while species_cur == "unclassified" and i < len(rows):
            # add those whose species is unknown
            s_cur, e_cur, species_cur, row_cur = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], list(rows[i])
            compressed.append(tuple(row_cur))
            i += 1

        for row in rows[i:]:
            s2, e2, species2, row2 = int(row[qs]), int(row[qe]), row[rsp], list(row)
            # if this new species is unclassified, just add it and move on
            if species2 == "unclassified":
                compressed.append(tuple(row2))
                continue

            if species2 == species_cur and overlaps(s_cur, e_cur, s2, e2):
                #then update the end
                e_cur = max(e_cur, e2)
            #if it does not overlap or does not share the same species, start the next merge
            else:
                #update this row to be the new values of the new end
                row_cur[qe] = str(e_cur)
                compressed.append(tuple(row_cur))
                #start a new current fam
                s_cur, e_cur, species_cur, row_cur = int(row[qs]), int(row[qe]), row[rsp], list(row)

        # add the final interval
        row_cur[qe] = str(e_cur)
        compressed.append(tuple(row_cur))

    return compressed

def build_overlap_row(row1, row2, new_start, new_end, ani_value):
    new_row = []
    for h in range(len(row1)):
        if h in [0, 1, 9]:  # qname, rname, q_species
            new_row.append(row1[h])
        elif h == 2:
            new_row.append(str(new_start))
        elif h == 3:
            new_row.append(str(new_end))
        else:
            new_row.append(f"{row1[h]},{row2[h]}")
    new_row.append(str(ani_value))
    return "\t".join(new_row)

def find_overlap_and_div_max(rows, output_file, tree, triangle_dict, index):
    """
    Same contract as find_overlap_and_div but pairs intervals with the
    largest mutual overlap first.

    rows         - list of care-tuples from compress()
    output_file  - path to write the overlap-div TSV
    tree         - Divergence_Tree_Preprocessed instance
    triangle_dict - (min_id, max_id) -> True for all db pairs with ANI >= 95
    index        - hash table from load_hash_table()
    """
    div_cache = {}
    ani_cache = {}
    new_rows  = set()

    # Group by Q name so coordinates from different query sequences never mix.
    groups = defaultdict(list)
    for row in rows:
        groups[row[0]].append(row)

    for qrows in groups.values():
        _pair_group(qrows, new_rows, tree, triangle_dict, index,
                    div_cache, ani_cache)

    with open(output_file, "w") as out:
        out.write(_HEADER)
        for line in new_rows:
            out.write(line + "\n")


def _pair_group(rows, new_rows, tree, triangle_dict, index,
                div_cache, ani_cache):
    """Process all rows for a single Q name.

    Two optimisations over the naive O(n*k) approach:

    1. Lazy quality check — tree.divergence / lookup_ani_triangle are only
       called at *pop* time, not during the sweep.  Pairs whose rows get
       consumed by a larger-overlap pair are popped, skipped via `used`, and
       never pay for the expensive check at all.


    2. Per-row sweep cap — the inner loop is limited to
       max_candidates_per_row overlapping rows per row i.  Because rows are
       sorted by Q start the closest-starting (highest-overlap) candidates
       are always seen first, so the cap rarely drops a real best pair while
       bounding worst-case complexity to O(n * max_candidates_per_row).
    """
    rows = sorted(rows, key=lambda x: (int(x[_qs]), int(x[_qe])))
    n = len(rows)
    if n <= 1:
        return

    # ---------------------------------------------------------------------- #
    # Sweep: find all geometrically overlapping pairs.                        #
    # Only apply the cheap species-string equality check here.                #
    # Everything else (leaf lookup, div, ANI) is deferred to the pop phase.  #
    # ---------------------------------------------------------------------- #
    heap = []  # (-overlap_size, i, j) — no pair_meta needed

    for i, row_i in enumerate(rows):
        s1  = int(row_i[_qs])
        e1  = int(row_i[_qe])
        sp1 = row_i[_rsp]
        best_ov = 0  # largest overlap already added to heap for row i

        for j in range(i + 1, n):
            s2 = int(rows[j][_qs])
            if s2 > e1:
                break  # sorted by start; nothing further can overlap row i

            # Upper bound on overlap for this j and ALL future j' (s2 only increases):
            # overlap(i, j') ≤ e1 - s2'.  Once that ceiling ≤ best_ov, break.
            if e1 - s2 < best_ov:
                break

            ov = overlaps(s1, e1, s2, int(rows[j][_qe]))
            if not ov:
                continue

            # cheapest possible early exit — same non-unclassified species string
            if sp1 == rows[j][_rsp] and sp1 != "unclassified":
                continue

            heapq.heappush(heap, (-ov, i, j))
            if ov > best_ov:
                best_ov = ov
            if best_ov == e1 - s1:
                break  # full-coverage pair found — nothing better possible for row i

    # ---------------------------------------------------------------------- #
    # Pop in decreasing overlap order.                                        #
    # Expensive quality check runs only for pairs where both rows are still   #
    # available — pairs consumed by a larger overlap are skipped for free.    #
    # ---------------------------------------------------------------------- #
    used = set()
    while heap:
        if len(heap) % 10000000 == 0: print(len(heap))
        _neg_ov, i, j = heapq.heappop(heap)
        if i in used or j in used:
            continue  # already consumed — skip with zero quality-check cost

        row1, row2 = rows[i], rows[j]
        s1 = int(row1[_qs]); e1 = int(row1[_qe])
        s2 = int(row2[_qs]); e2 = int(row2[_qe])
        sp1, sp2 = row1[_rsp], row2[_rsp]

        sp1_leaf = lookup_tree_leaf_name(index, row1[_tname])
        sp2_leaf = lookup_tree_leaf_name(index, row2[_tname])

        if (sp1_leaf == sp2_leaf or sp1 == sp2) and \
           "unclassified" not in (sp1, sp2):
            continue

        div = check_cache(sp1_leaf, sp2_leaf, div_cache)
        if div is None:
            div = tree.divergence(sp1_leaf, sp2_leaf)
            div_cache[(sp1_leaf, sp2_leaf)] = div

        ani = "NA"
        if isinstance(div, str):
            id1, id2 = row1[_tname], row2[_tname]
            ani = check_cache(id1, id2, ani_cache)
            if ani is None:
                ani = lookup_ani_triangle(
                    id1, id2, triangle_dict, index, s1, e1, s2, e2
                )
                ani_cache[(id1, id2)] = ani

        if not ((not isinstance(div, str) and div >= 1) or ani is True):
            continue  # pair fails quality filter

        new_start = max(s1, s2)
        new_end   = min(e1, e2)
        out_row = []
        for h in range(len(row1)):
            if h in (0, 1, 9):
                out_row.append(row1[h])
            elif h == _qs:
                out_row.append(str(new_start))
            elif h == _qe:
                out_row.append(str(new_end))
            else:
                out_row.append(f"{row1[h]},{row2[h]}")
        out_row.extend([str(div), str(ani)])

        used.add(i)
        used.add(j)
        new_rows.add("\t".join(out_row))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage: python3 find_overlap_and_div_max.py "
            "<first_div.tsv> <output.tsv> <tree_dir> <medival_db>"
        )
        sys.exit(1)

    input_file  = sys.argv[1]
    output_file = sys.argv[2]
    tree        = Divergence_Tree_Preprocessed(sys.argv[3])
    medival_db  = sys.argv[4]
    index       = load_hash_table(f"{medival_db}/medival_db_index.pkl")
    with open(f"{medival_db}/skani_triangle_ani95.pkl", "rb") as f:
        triangle_dict = pickle.load(f)

    rows = compress(input_file)
    find_overlap_and_div_max(rows, output_file, tree, triangle_dict, index)
    print("Filtered file written to", output_file)
