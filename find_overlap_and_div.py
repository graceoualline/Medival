#same as find overlap, but we also filter for if the two lines we are compressing 
# map to two different speceis that are divergently different

#this will only return regions that have been repeated, 
# and are different species

# will be O(N^2)
import csv
import sys
from filter_blat import *
import os
from build_database_index import *
from find_overlap import *

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

#if it finds regions that overlap and map to the same species, it will combine them

def find_overlap_and_div(rows, output_file, tree, skani_db, index):
    div_cache = dict()
    ani_cache = dict() # (id1, id2) : ani
    #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7
    tname = 4
    new_rows = set()
    rows = sorted(list(rows), key=lambda x:int(x[qs]))
    used = set()

    i = 0

    for i, row1 in enumerate(rows):
        if i in used:
            continue

        s1, e1, species1 = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp]
        
        if i + 1 < len(rows) and e1 < int(rows[i+1][qs]):
            continue
        
        for j in range(i+1, len(rows)):
            if j in used:
                continue
            s2, e2, species2, row2 = int(rows[j][qs]), int(rows[j][qe]), rows[j][rsp], rows[j]
            # if we completely surpassed where we can overlap, skip
            if s2 > e1:
                break

            if overlaps(s1, e1, s2, e2) == None:
                continue

            # if they are both equal and known
            sp1_leaf = lookup_tree_leaf_name(index, row1[tname])
            sp2_leaf = lookup_tree_leaf_name(index, row2[tname])

            # if the leafs are the same and they are both classfied, then skip
            if (sp1_leaf == sp2_leaf or species1 == species2) and "unclassified" not in [species1, species2]:
                continue
            div = check_cache(sp1_leaf, sp2_leaf, div_cache)
            if div == None:
                # Get the divergence time between query species and reference species
                div = tree.divergence(sp1_leaf, sp2_leaf)
                div_cache[(sp1_leaf, sp2_leaf)] = div
            
            ani = "NA"

            # if div is not a number, then it is unknown
            if isinstance(div, str):
                #get the ani instead
                id1 = row1[tname]
                id2 = row2[tname]
                ani = check_cache(id1, id2, ani_cache)
                if ani == None:
                    ani = find_ani_overlap(id1, id2, skani_db, index)
                    if ani == "unk:too_short": ani = crude_ani_overlap(s1, e1, s2, e2, int(row1[5]), int(row2[5]))
                    ani_cache[(id1, id2)] = ani

                
            if (not isinstance(div, str) and div >= 1) or (isinstance(ani, int) and ani < 95):
                #print("merging rows")
                new_row = []
                new_start = max(s1, s2)
                new_end = min(e1, e2)
                for h in range(len(row1)):
                    if h in [0, 1, 9]: #0 and 1 are q name and length, and q species
                        new_row.append(row1[h])
                    elif h == 2: new_row.append(str(new_start)) #put new start or end
                    elif h==3: new_row.append(str(new_end)) # put in new end
                    else:
                        new_row.append(f"{row1[h]},{row2[h]}")
                new_row.extend([str(div), str(ani)])
                #print(i, j, "merged")
                used.add(rows[i])
                used.add(rows[j])
                new_rows.add("\t".join(new_row))
                break
            
    with open(output_file, "w") as out:
                 #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tDiv bt ref species\tAni bt ref species\n")
        #if we do have new rows, then write them
        if len(new_rows) > 0:
            for l in new_rows:
                out.write(l + "\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    #
    if len(sys.argv) != 5:
        print("Usage: python3 find_overlap_and_div.py <file of blatdiver output> <output file name> <tree> <medival_db>")
        sys.exit(1)
    input_file = sys.argv[1]
    output = sys.argv[2]
    tree = Divergence_Tree_Preprocessed(sys.argv[3])
    index = load_hash_table(f"{sys.argv[4]}/medival_db_index.pkl")
    skani_db = f"{sys.argv[4]}/skani_db"

    #print("Compressing")
    rows = compress(input_file)
    #print("Rows before:", len(rows))
    find_overlap_and_div(rows, output, tree, skani_db, index)
    print("Filtered file written to", output)