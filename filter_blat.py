#this will intake a psl blat file and only output the results where there is divergence between the two species
# will also add a column that will show the divergence time between the query and reference species, and those two species

#assuming we figured out the query species in the beginning using Kraken

import sys
import pickle
from build_database_index import *

def crude_ani(len1, len2, matches):
    if 0 in [len1, len2]: return -1
    value = matches / min(len1, len2)
    return value*100

def check_cache(a, b, cache):
    if (a, b) in cache: return cache[(a, b)]
    elif (b, a) in cache: return cache[(b, a)]
    return None

def filter_blat(inf, outf, q_species, index, tree, q_seq, blat_db, minIdentity, skani_ani_dict):
    div_cache = dict() # will be (species1, species2): divergence time

    # get the name of the q species that is in the tree
    q_og_species = q_species
    if q_species != "unclassified": q_species = tree.get_path_identifier(tree.get_path(q_species))

    with open(inf, 'r') as infile, open(outf, 'w') as outfile:

        outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \tQ gap bases\tT gap count\tT gap bases\tstrand\tQ name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI<95(if div=unk)\n")
        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) != 21 or columns[0] == 'match':
                continue
            query_id = columns[9]
            ref_id = columns[13]
            match = int(columns[0])
            Qstart = int(columns[11])
            Qsize = int(columns[10])
            Qend = int(columns[12])
            perIdent = round((match / (Qend - Qstart)) * 100, 5)

            if perIdent < minIdentity:
                continue

            ref_species = lookup_species(index, ref_id)
            if ref_species is None:
                print("ERROR, your index mapping database genomes to their species and locations was not made correctly.")
                print(f"Unable to find the species for {ref_id} in {index}. Skipping.")
                continue

            if q_species in ["unclassified", "NA"] or ref_species in ["unclassified", "NA"]:
                div = "unk:unclassified_species"
            else:
                div = check_cache(q_species, ref_species, div_cache)
                if div is None:
                    ref_leaf_name = lookup_tree_leaf_name(index, ref_id)
                    div = tree.divergence(q_species, ref_leaf_name)
                    div_cache[(q_species, ref_species)] = div
            if div is None:
                div = "unk:unable_to_find_ref_species_in_tree"

            if isinstance(div, (int, float)) and div < 1:
                continue  # Skip lines with less than 1 MYA

            if isinstance(div, str):
                ref_len = lookup_length(index, ref_id)
                # Fall back to crude ANI when either sequence is too short for skani
                if Qsize < 500 or ref_len < 500:
                    ani = crude_ani(Qsize, ref_len, match)
                    ani_under_95 = isinstance(ani, (int, float)) and 0 <= ani < 95
                else:
                    # Use precomputed skani search results: presence means ANI >= 95
                    ani_under_95 = ref_id not in skani_ani_dict.get(query_id, [])

                if not ani_under_95:
                    continue  # ANI >= 95, discard
                new_line = "\t".join(columns)
                outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{ref_species}\t{div}\t{ani_under_95}\n")
            else:
                new_line = "\t".join(columns)
                outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{ref_species}\t{div}\tNA\n")

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python3 filter_blat.py <input psl file> <output .tsv file> <query species> <medival_db> <div tree> <input sequence> <minIdentity> <skani_ani_dict.pkl>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    q_species = sys.argv[3]
    medival_db = sys.argv[4]
    tree = Divergence_Tree_Preprocessed(sys.argv[5])
    q_seq = sys.argv[6]
    index = load_hash_table(f'{medival_db}/medival_db_index.pkl')
    blat_db = f'{medival_db}/blat_2bit_db'
    minIdentity = float(sys.argv[7])
    with open(sys.argv[8], 'rb') as f:
        skani_ani_dict = pickle.load(f)
    filter_blat(input_file, output_file, q_species, index, tree, q_seq, blat_db, minIdentity, skani_ani_dict)