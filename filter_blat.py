#this will intake a psl blat file and only output the results where there is divergence between the two species
# will also add a column that will show the divergence time between the query and reference species, and those two species

#assuming we figured out the query species in the beginning using Kraken

import pandas as pd
import sys
from ete3 import Tree
import subprocess
#import glob
import os
import tempfile
#import bisect
from build_database_index import *
from Bio import SeqIO

def crude_ani(len1, len2, matches):
    if 0 in [len1, len2]: return -1
    return (matches / min(len1, len2))*100

def batch_extract_sequences(ref_ids, index, skani_db, cutoff = 500):
    """
    Make a file that has a list all of the relevant sequences
    And returns the sequences that are too short
    """
    too_short = []
    # a dictionary opf [ path to ref] = ref
    long_refs = []
    # make a temp list file:
    with tempfile.NamedTemporaryFile(prefix = f"reference_list", suffix=".txt", delete=False) as ref_fasta:
        ref_fasta_list = ref_fasta.name
    with open(ref_fasta_list, "w") as f:
        for ref in ref_ids:
            length = lookup_length(index, ref)
            if length < cutoff: too_short.append(ref)
            else:
                long_refs.append(ref)
                f.write(f"{skani_db}/{ref}.fasta\n")

    return ref_fasta_list, long_refs, too_short

def batch_calculate_ani(query_path, ref_list_path, ref_list):
    """
    Calculate ANI for all query-ref pairs using skani in batch mode.
    Returns dict of ref_id -> ANI value
    """
    if not ref_list_path:
        return {}
    
    ani_results = {}
    
    try:
        # Run skani with query against all references
        print(f"skani dist {query_path} --rl {ref_list_path}")
        result = subprocess.run(
            ["skani", "dist", query_path, "--rl", ref_list_path],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse output
        lines = result.stdout.strip().split("\n")
        for line in lines[1:]:  # Skip header
            parts = line.split()
            if len(parts) >= 3:
                ref_path = parts[0]
                ref = os.path.splitext(ref_path)[0]
                ani_value = float(parts[2])
                ani_results[ref] = ani_value
        # Return 0 for the ones that skani didnt calculate
        for r in ref_list:
            if r not in ani_results:
                ani_results[r] = 0
                
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"Error in batch ANI calculation: {e}")
        #with open("debugging_batch_skani_err.txt", "w") as f: f.write(str(e))
        #input("Press enter to continue")
        # Return -1 for all refs that failed
        for ref_id in ref_list_path.keys():
            if ref_id not in ani_results:
                ani_results[ref_id] = "unk:skani_err"
    finally:
        os.remove(ref_list_path)
    
    return ani_results

def calculate_distance(seq1, seq2):
    """
    Calculates the ANI distance between two sequences using skani.
    
    Parameters:
    - seq1 (str): Path to the query sequence (FASTA format).
    - seq2 (str): Path to the reference sequence (FASTA format).

    Returns:
    - float: ANI value (or 0 if not found).
    """
    #skani wont work if sequence under 500 bp
    try:
        # Run the skani command to calculate distance
        result = subprocess.run(
            ["skani", "dist", seq1, seq2],
            capture_output=True,
            text=True,
            check=True
        )
        # Capture the output and process it to extract the ANI value
        output = result.stdout
        lines = output.strip().split("\n")
        
        # Assume the second line contains the data, split by whitespace to extract the ANI
        if len(lines) > 1:
            ani_value = lines[1].split()[2]
            return float(ani_value)
        else:
            #that means nothing was found, so the ani is 0
            return 0

    except (IndexError, ValueError, subprocess.CalledProcessError) as e:
        print(f"Error running skani: {str(e)}")
        print(f"Common issue is the fasta files given to skani do not exist")
        print("Please make sure the skani_db portion of the medival database was create correctly.")
        return "unk:skani_err"  # Return 0 in case of errors

def check_cache(a, b, cache):
    if (a, b) in cache: return cache[(a, b)]
    elif (b, a) in cache: return cache[(b, a)]
    return None

def filter_blat(inf, outf, q_species, index, tree, q_seq, blat_db, minIdentity, skani_db):
    div_cache = dict() # will be (species1, species2): divergence time

    # get the name of the q species that is in the tree
    q_og_species = q_species
    if q_species != "unclassified": q_species = tree.get_path_identifier(tree.get_path(q_species))
    
    lines_to_process = []
    ref_ids_needing_ani = set()

    with open(inf, 'r') as infile, open(outf, 'w') as outfile:
        
        outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \tQ gap bases\tT gap count\tT gap bases\tstrand\tQ name\tQ size\tQ start\t Q end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\n")
        # Process each line in the BLAT file
        for line in infile:
            columns = line.strip().split('\t')
            #skip the header lines
            if len(columns) != 21 or columns[0] == 'match':
                continue
            query_id = columns[9]  # Assuming query sequence ID is in the 10th column
            ref_id = columns[13]  # Assuming reference sequence ID is in the 14th column
            match = int(columns[0])
            Qstart = int(columns[11])
            Qend = int(columns[12])
            perIdent = round((match / (Qend-Qstart))*100, 5)

            #skip this line if it doesnt meet the threshold
            if perIdent < minIdentity:
                continue
            # Get the species of the reference sequence using Kraken and grep
            ref_species = lookup_species(index, ref_id)
    
            if ref_species is None:
                print("ERROR, your index mapping database genomes to their species and locations was not made correctly.")
                print(f"Unable to find the species for {ref_id} in {index}. Skipping.")
                continue  # Skip this line if no species found for reference sequence
            if q_species in ["unclassified", "NA"] or ref_species in ["unclassified", "NA"]: # or path1 == None:
                div = "unk:unclassified_species"
            else:
                div = check_cache(q_species, ref_species, div_cache)
                if div == None:
                    ref_leaf_name = lookup_tree_leaf_name(index, ref_id)
                    # Get the divergence time between query species and reference species
                    div = tree.divergence(q_species, ref_leaf_name)
                    div_cache[(q_species, ref_species)] = div
            if div == None:
                div = "unk:unable_to_find_ref_species_in_tree"

            ani = "NA" #ani is NA for now unless defined later
            if isinstance(div, (int, float)) and div < 1:
                continue  # Skip lines with less than 1 MYA
            elif isinstance(div, str): #if div is unk, then find ani
                lines_to_process.append({
                'columns': columns,
                'ref_id': ref_id,
                'perIdent': perIdent,
                'ref_species': ref_species,
                'div': div
                })
                ref_ids_needing_ani.add(ref_id)
                continue

            # write everyone that made it through the filters
            new_line = "\t".join(columns)
            outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{ref_species}\t{div}\t{ani}\n")
        
        # now go through everyone that needs ani taken
        if ref_ids_needing_ani:
            print(f"Extracting {len(ref_ids_needing_ani)} reference sequences...")
            ref_fasta_list, long_refs, too_short = batch_extract_sequences(ref_ids_needing_ani, index, skani_db)
            
            print(f"Calculating ANI for {len(long_refs)} sequences...")
            ani_results = batch_calculate_ani(q_seq, ref_fasta_list, long_refs)

            # then write them to file if they pass
            for line_data in lines_to_process:
                columns = line_data['columns']
                ref_id = line_data['ref_id']
                perIdent = line_data['perIdent']
                div = line_data['div']
                
                ani = "NA"
                ani = ani_results.get(ref_id, f"unk:too_short")
                # if one of the sequences is too short, just do a crude ani
                if ani == "unk:too_short": ani = crude_ani(int(columns[10]), int(columns[14]), int(columns[0]))
                if isinstance(ani, (int, float)):
                    if ani < 0 or ani >= 95:
                        continue
                # Write the line with the added divergence time, query species, and reference species
                new_line = "\t".join(columns)
                outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{lookup_species(index, ref_id)}\t{div}\t{ani}\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 8:
        print("Usage: python3 filter_blat.py <input psl file> <output .tsv file> <query species (make sure _ instead of space)> <medival_db> <div tree> <input sequence> <minIdentity>")
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
    skani_db = f'{medival_db}/skani_db'
    filter_blat(input_file, output_file, q_species, index, tree, q_seq, blat_db, minIdentity, skani_db)