#this will intake a kraken output file
# and return a .txt list of all sequence id's and species that were found
import tempfile
import os
import subprocess
import sys
from pathlib import Path

def process_all_sequences_kraken(kraken_file):
    # Use grep to find the sequence ID in the Kraken file
    all_species = set()
    with open(kraken_file, "r") as infile:
        # Extract the species name from the grep output
        # The species name is in the 3rd column (assuming Kraken output format is consistent)

        for line in infile:
            id = line.split('\t')[1]
            species = line.split('\t')[2].split('(')[0].strip()
            #print(id, species)
            all_species.add((id, species))
    return all_species

def get_sp_from_kraken(kraken_output_file):
    # Use grep to find the sequence ID in the Kraken file
    try:
        with open(kraken_output_file, "r") as infile:
            # Extract the species name from the grep output
            # The species name is in the 3rd column (assuming Kraken output format is consistent)
            lines = infile.readlines()
            if not lines:
                print("Kraken output is empty.")
                return "unclassified"

            last_line = lines[-1]
            last_line = last_line.split()
            species = last_line[5:]
            species = "_".join(species)
            return species
    except Exception as e:
            print(f"Error parsing Kraken output: {e}")
            return "unclassified"

def get_q_species(fasta_file, kraken_db):
    with tempfile.NamedTemporaryFile(prefix="kraken_output", delete=False, mode='w+') as tmp:
        tmp_filename = tmp.name
    command = [
        "kraken2", "--db", kraken_db, "--threads", "1",
        "--report", tmp_filename, fasta_file
    ]
    try:
        #run the command
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        #print("Done running kraken")
        species = get_sp_from_kraken(tmp_filename)
    
    except subprocess.CalledProcessError as e:
        print(f"Error running Kraken on {fasta_file}:\n{e.stderr.decode().strip()}")
        species = "unclassified"
    finally:
        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)

    return species

def make_seq_id_species_txt(seq_file, output, kraken_db):
    kraken_output = output + "_kraken_output.txt"

    command = [
        "kraken2", "--db", kraken_db, "--output", kraken_output, "--use-names", seq_file
    ]
    #print(command)
    try:
        #run the command
        print("Running Kraken on", seq_file)
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Kraken completed. Kraken output written to {kraken_output}. \nExtracting sequence id's and species...")
        species = process_all_sequences_kraken(kraken_output)
        
        #write all of the species to .txt file
        species_file = output + "_seq_species.txt"
        print("Writing sequence id's and species to", species_file)
        with open(species_file, 'w') as f:
            for seq, spec in species:
                f.write(f"{seq}\t{spec}\n")
        print(f"Data successfully written to {species_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Kraken on {seq_file}:\n{e.stderr.strip()}")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    # if q_species = unk, then run the get query species
    if len(sys.argv) != 3:
        print("Usage: python3 extract_species_from_kraken.py <input fasta file> <kraken db> ")
        sys.exit(1)
    q_seq = Path(sys.argv[1])
    output = q_seq.stem
    kraken_db = sys.argv[2]    
    make_seq_id_species_txt(q_seq, output, kraken_db)