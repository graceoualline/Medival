# this should intake the folder that contains a bunch of split of 2bit files
# and it will return a .pkl file which is a hash table that essential works like
# index[sequence_id] = (sequence_species, sequence_location_in_the_database, name_in_tree)
from ete3 import Tree
import os
import subprocess
import tempfile
import pickle
from extract_species_from_kraken import *
import pandas as pd
import numpy as np

def get_one(path_iter):
    paths = list(path_iter)
    if len(paths) != 1:
        raise ValueError(f"Expected 1 file, found {len(paths)}: {paths}. You do not have the proper files needed for the divergence tree.")
    return str(paths[0])

def load_np(file):
    return np.load(file, allow_pickle=True)

class Divergence_Tree_Preprocessed():
    """
    Class for finding divergence along a tree that has been preprocessed
    """
    def __init__(self, tree_dir):
        tree_path = Path(tree_dir)
        nwk_file   = get_one(tree_path.glob("*.nwk"))
        plk_file   = get_one(tree_path.glob("*.pkl"))
        index_file = get_one(tree_path.glob("*.index.npy"))
        mins_file  = get_one(tree_path.glob("*.mins.npy"))
        tour_file  = get_one(tree_path.glob("*.tour.npy"))

        print("Loading Divergence Tree")
        self.tree = Tree(nwk_file)
        self.species_dist = load_hash_table(plk_file)
        self.index = load_np(index_file)
        self.mins = load_np(mins_file)
        self.tour = load_np(tour_file)
        print("Divergence Tree Loaded!")
        
    def get_min_index(self, i,j):
        """Range min query on L[i:j]"""
        m = int(np.log2(j-i))
        minimum = min(self.mins[m, i], self.mins[m, j-2**m])
        min_index = self.index[m, i] if self.mins[m, i] == minimum else self.index[m, j-2**m]
        return min_index

    def divergence(self, a, b):
        """Finds divergence between nodes a and b"""
        if a == "NA" or b == "NA":
            return "unk:unable_to_find_ref_species_in_tree"
        a_tour_step = self.species_dist[a][1]
        b_tour_step = self.species_dist[b][1]
        a_dist = self.species_dist[a][0]
        b_dist = self.species_dist[b][0]
        lca_step = self.get_min_index(min(a_tour_step, b_tour_step), max(a_tour_step, b_tour_step)+1)
        lca = self.tour[lca_step]
        lca_dist = self.species_dist[lca][0]
        return (a_dist + b_dist - 2*lca_dist)/2

    def get_path(self, sp):
        if sp == "unclassified":
            return "NA"

        sp = "_".join(sp.split(" "))
        sp = "'" + sp + "'"

        try:
            path = self.tree & sp
            return path
        except:
            first = sp.split("_")
            try:
                new = first[0] + "_" + first[1] + "'"
                path = self.tree & new
                return path
            except:
                try: 
                    new = first[0] + "'"
                    path = self.tree & new
                    return path
                except:
                    #if all else fails, find the first instance of the name within another name
                    for node in self.tree.traverse("preorder"):
                        if first[0] in node.name:
                            return node
                    return "NA"

    # return the actual node name that is within the tree
    def get_path_identifier(self, node):
        if node == "NA":
            return "NA"
        path_names = []
        current = node
        while current is not None:
            path_names.append(current.name if current.name else "unnamed")
            current = current.up
        return path_names[0]

def get_seq_id_length(seq_lengths):
    df = pd.read_csv(str(seq_lengths), header = None, sep = "\t")
    return list(df.itertuples(index=False, name=None))

def combine_seq_spec_len(species_file, seq_id_len, tree):
    # Step 1: Load species data into a dictionary
    seq_to_species = {}
    all_species = set()
    if species_file and os.path.exists(species_file):
        print(f"Loading existing species data from {species_file}...")
        with open(species_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    seq_id, species = parts[0], parts[1]
                    all_species.add(species)
                    seq_to_species[seq_id] = species
        print(f"Loaded {len(seq_to_species)} existing species annotations")
    else:
        raise Exception("No species file provided or file doesn't exist")
    
    #Step 2: Create a dictionary of all species, and where they are in the tree
    print("Obtaining species path in the tree")
    print(f"There are {len(all_species)} species")
    species_leaf_name = dict()
    for i, sp in enumerate(all_species):
        if i % 50 == 0:
            print(f"{i}/{len(all_species)} paths found.")
        path = tree.get_path(sp)
        leaf_name = tree.get_path_identifier(path)
        species_leaf_name[sp] = leaf_name
    # print("Species leaf name dictionary", species_leaf_name)
    print("Paths obtained. Creating hash table.")
    #Step 3: Process each sequence
    result_dict = {}
    for i, seq_len in enumerate(seq_id_len, 1):
        #make sure there are two entries
        if len(seq_len) >= 2:
            seq_id, length = seq_len[0], seq_len[1]

            # Progress indicator
            if i % 100 == 0 or i == len(seq_id_len):
                print(f"Progress: {i}/{len(seq_id_len)} sequences processed")

            # if the sequence does not have a species, go find it with kraken
            species = seq_to_species.get(seq_id, "unclassified")
            result_dict[seq_id] = (species, length, species_leaf_name[species])
    return result_dict

def load_hash_table(index_file):
    """
    Loads the precomputed hash table from a binary file.
    :param index_file: Path to the stored hash table
    :return: Dictionary {seq_id: species}
    """
    with open(index_file, "rb") as f:
        return pickle.load(f)

def lookup_species(hash_table, seq_id):
    """
    Looks up the species name for a given seq_id using the hash table.
    Returns "unclassified" if not found.
    """
    tuple = hash_table.get(seq_id, "unclassified")
    return tuple[0]

def lookup_length(hash_table, seq_id):
    tuple = hash_table.get(seq_id, 0)
    return int(tuple[1])

def lookup_tree_leaf_name(hash_table, seq_id):
    """
    Looks up the species name for a given seq_id using the hash table.
    Returns "unclassified" if not found.
    """
    tuple = hash_table.get(seq_id, "NA")
    return tuple[2]

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 build_database_index.py <medival_built_db> <species_file> <tree>")
        print("  the medival_db that was built by make_medival_db.py")
        print("  species_file: Existing species annotation file")
        print("  tree: Time Tree of Life .nwk file")
        sys.exit(1)
    medival_db = sys.argv[1]
    species_file = sys.argv[2]

    seq_lengths = Path(f"{sys.argv[1]}/seq_lengths.tsv")
    if not seq_lengths.exists():
        raise FileNotFoundError(f"File does not exist: {seq_lengths}, please make sure make_medival_db.py was run correctly, or manually make a file with the first column being seq_id, and second column seq_length.")
    final_output = f"{sys.argv[1]}/medival_db_index.pkl"
    tree = Divergence_Tree_Preprocessed(sys.argv[3])

    print("\n" + "="*60)
    print("STEP 1: Extracting sequence IDs and their lengths")
    print("="*60)
    seq_id_length = get_seq_id_length(seq_lengths)
    #seq_id_loc = [('JAHAAN010000031.1', 'f1'), ('QHWJ01000176.1', 'f2'), ('WQSI01000087.1', 'f3'), ('NZ_QJUG01000292.1', 'f4'), ('JAJXTT010000076.1', 'f5'), ('JAADEE010000229.1', 'f6')]
    # next, we use a file that connects seq_id and species determined by kraken
    # and add that to the index # used to be combine_species_location.py
    print(f"{seq_lengths} loaded.")
    print("\n" + "="*60) 
    print("STEP 2: Combining with species data and Kraken classification")
    print("="*60)
    seq_species_len = combine_seq_spec_len(species_file, seq_id_length, tree)

    # # now make a loadable index from this new table of id, species, len, leaf name
    print("\n" + "="*60)
    print("STEP 3: Saving database index")
    print("="*60)
    with open(final_output, "wb") as f:
        pickle.dump(seq_species_len, f)

    print(f"Hash index saved to {final_output}")
    print("="*60)
    print("DATABASE INDEX CREATION COMPLETED!")
    print(f"Once you've ensured your build database and index is functional, you can delete {medival_db}/blat_fasta_db and {seq_lengths} to save space.")
    print("="*60)
