import os
import pickle
import tempfile
import subprocess
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from datetime import datetime
from functools import partial
from extract_species_from_kraken import *
from tqdm import tqdm
from combine_medival_output import *
from find_overlap_and_div_max import *
from parse_args import *
from build_database_index import *
import time
from concurrent.futures import ThreadPoolExecutor
import threading
from size_cluster_filter import *

_worker_tree = None
_worker_index = None
_worker_triangle_dict = None

# initialize worker so dont need to reload tree and hash table every time
# _worker_triangle_dict is inherited via fork — do not pass via initargs
def _init_worker(tree_dir, index_path):
    global _worker_tree, _worker_index
    devnull = open(os.devnull, 'w')
    sys.stdout = devnull
    _worker_tree = Divergence_Tree_Preprocessed(tree_dir)
    _worker_index = load_hash_table(index_path)
    sys.stdout = sys.__stdout__
    devnull.close()

# surpressed what workers are printing out
def _suppress_stdout(func, args):
    old_stdout = sys.stdout
    devnull = open(os.devnull, 'w')
    sys.stdout = devnull
    try:
        return func(args)
    finally:
        sys.stdout = old_stdout
        devnull.close()

def combine_and_cleanup_psl_files(args):
    header_lines = 5
    output_chunk_dir, combined_path, remove, config_header = args

    if os.path.exists(combined_path):
        print(f"Skipping {combined_path}, already exists.")
    else:
        input_files = glob.glob(os.path.join(output_chunk_dir, "*.psl"))

        with open(combined_path, "w") as outfile:
            if config_header:
                outfile.write(config_header)
            outfile.write("match\tmismatch\trep. match\tN's\tQ gap count\tQ gap bases\tT gap count\tT gap bases\tstrand\tQ name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\n")
            for fname in input_files:
                with open(fname) as infile:
                    lines = infile.readlines()
                    # Now find where the alignment lines begin
                    outfile.writelines(lines[header_lines:])
                if remove: os.remove(fname)
        if remove: os.rmdir(output_chunk_dir)

def run_div_filter(job):
    all_blat_raw, output_file, species, index, tree, tmp_fasta_path, database, minIdentity, skani_ani_dict = job
    #index = load_hash_table(index)
    try:
        #run the first round of filtering on the raw blat data
        filter_blat(all_blat_raw, output_file, species, _worker_index, _worker_tree, tmp_fasta_path, database, minIdentity, skani_ani_dict)
    except subprocess.CalledProcessError as e:
        print(f"Error running filter_blat on chunk {output_file}:\n{e.stderr.decode().strip()}")
    return 1

def run_specified_filter(job):
    job_type, args = job
    try:
        if job_type == "overlap_div":
            to_filter_file, overlap_div_file = args
            rows = compress(to_filter_file)
            find_overlap_and_div_max(rows, overlap_div_file, _worker_tree, _worker_triangle_dict, _worker_index)
    except subprocess.CalledProcessError as e:
        print(f"[{job_type}] Error on {args[1] if len(args) > 1 else 'unknown'}:\n{e.stderr.decode().strip()}")
    return 1

def run_all_filters(chunk_jobs, c: Config):
    blat_combine_jobs = []
    jobs = []
    for args in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = args
        output_chunk_dir = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}")
        
        #make one big file that combines all of the blat outputs
        all_blat_raw = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if not(os.path.exists(all_blat_raw)): blat_combine_jobs.append((output_chunk_dir, all_blat_raw, c.remove, c.config_header))
        
        #then filter all blat results and add divergence or ani
        output_file = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}_first_div_output.tsv")
        jobs.append((chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c))

    output_name = os.path.basename(c.output_dir)
    blat_final_name = os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")
    overlap_div_final_name = os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")
    medival_final_name = os.path.join(c.output_dir, f"{output_name}_first_div_output.tsv")
    
    #first, combine all blat files
    if os.path.exists(blat_final_name): print(f"Skipping {blat_final_name}, already exists.")
    else:
        if len(blat_combine_jobs) > 0:
            with ThreadPoolExecutor(max_workers=c.max_threads) as executor:
                #Progress bar
                with tqdm(total=len(blat_combine_jobs), desc="Combining blat data") as pbar:
                    for _ in executor.map(lambda args: _suppress_stdout(combine_and_cleanup_psl_files, args), blat_combine_jobs):
                        pbar.update()
        print("Finished combining blat data.")
    
    #run divergence filter
    div_jobs = []
    if os.path.exists(medival_final_name): print(f"Skipping {medival_final_name}, already exists.")
    else:
        for j in jobs:
            chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = j
            if os.path.exists(output_file): print(f"Skipping {output_file}, already exists.")
            else: div_jobs.append((all_blat_raw, output_file, species, c.index, c.tree, tmp_fasta_path, c.blat_db, c.minIdentity, c.skani_ani_dict))
        if len(div_jobs) > 0:
            #print whatever command you are running
            with Pool(processes=c.max_threads, initializer=_init_worker, initargs=(c.tree, c.index)) as pool:
                #Progress bar
                with tqdm(total=len(div_jobs), desc="Running chunks through first divergence filter") as pbar:
                    for _ in pool.imap_unordered(run_div_filter, div_jobs):
                        pbar.update()
        print("Finished running divergence filter.")

    #run overlap_div filter
    make_overlap_div = c.overlap_div_filter
    if make_overlap_div and os.path.exists(overlap_div_final_name):
        print(f"Skipping {overlap_div_final_name}, already exists.")
        make_overlap_div = False
    filter_jobs = []
    for j in jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, medival_output_file, c = j

        # overlap div filter
        if make_overlap_div:
            overlap_div = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}_overlap_div.tsv")
            if os.path.exists(overlap_div): print(f"Skipping {overlap_div}, already exists.")
            else: filter_jobs.append(("overlap_div", (medival_output_file, overlap_div)))

    #run this first round of filtering
    if len(filter_jobs) > 0:
        global _worker_triangle_dict
        _worker_triangle_dict = c.skani_triangle_dict  # set before fork so workers inherit via copy-on-write
        with Pool(processes=c.max_threads, initializer=_init_worker, initargs=(c.tree, c.index)) as pool:
            #Progress bars
            with tqdm(total=len(filter_jobs), desc="Running chunks through overlap div filter") as pbar:
                for _ in pool.imap_unordered(run_specified_filter, filter_jobs):
                    pbar.update()
    print("All filters complete")
    
def run_blat(args):
    blat_dir, blat_file, ooc_file, query, output_path, minScore = args
    command = ["blat", str(blat_dir / blat_file), query,
           f"-ooc={blat_dir / ooc_file}", "-tileSize=11",
           f"-minScore={minScore}", output_path, "-q=dna", "-t=dna"]
    for attempt in range(3):
        try:
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running blat on {output_path} (attempt {attempt+1}/3):\n{e.stderr.decode().strip()}")
            continue
        if os.path.getsize(output_path) > 0:
            return 1
        print(f"Warning: blat produced empty output for {output_path} (attempt {attempt+1}/3), retrying...")
    print(f"Error: blat failed to produce output for {output_path} after 3 attempts")
    return 1

def run_blat_on_chunk(chunk_jobs, blat_files, ooc_files, c: Config):
    num_blat_db = len(blat_files)

    jobs = []
    output_name = os.path.basename(c.output_dir)
    final_blat_output = os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")
    # for each chunk, make all of the blat search jobs
    if os.path.exists(final_blat_output): 
        print(f"Skipping blat, {final_blat_output} already exists.")
        return
    for chunk in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = chunk
    
        # make a temporary file that isolates the fasta on its own
        #get the name and path for the output file
        all_blat_output = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if os.path.exists(all_blat_output): print(f"Skipping blat on chunk_{idx}_{original_id}, chunk_{idx}_{original_id}_blat_output.tsv already exists.")
        else:
            output_chunk = os.path.join(c.intermediate_dir, f"chunk_{idx}_{original_id}")
            if not os.path.exists(output_chunk): os.makedirs(output_chunk)
            
            for i in range(num_blat_db):
                output_path = f"chunk_{idx}_{original_id}_part_{i}.psl"
                output_path = os.path.join(output_chunk, output_path)
                if os.path.exists(output_path) and os.path.getsize(output_path) > 0: print(f"{output_path} skipped, already exists")
                else:
                    jobs.append((c.blat_db, blat_files[i], ooc_files[i], tmp_fasta_path, output_path, c.minScore))
    if len(jobs) > 0:
        from concurrent.futures import as_completed
        with ThreadPoolExecutor(max_workers=c.max_threads) as executor:
            with tqdm(total=len(jobs), desc="Running Sequence(s) through blat") as pbar:
                futures = [executor.submit(run_blat, job) for job in jobs]
                for _ in as_completed(futures):
                    pbar.update()

def run_combine_results(job):
    output_dir, chunk_size, to_combine, output_path, remove, config_header = job
    try:
        to_remove = adjust_and_merge_tsvs(output_dir, chunk_size, output_path, to_combine, config_header)
        if remove: 
            for f in to_remove: os.remove(f)
    except subprocess.CalledProcessError as e:
        print(f"Error combining results {output_path}: {e}")
    return 1

def combine_all_results(c: Config):
    #once all jobs are finished, combine all of the filters into one place:
    #parent_dir = os.path.dirname(c.output_dir)
    output_name = os.path.basename(c.output_dir)
    #python3 combine_medival_output.py <medival_output_directory> <chunk_size> <which files to combine ex. *first_div_output.tsv> <output_file>
    input_output = [("blat_output.tsv", os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")),
                    ("first_div_output.tsv", os.path.join(c.output_dir, f"{output_name}_first_div_output.tsv"))]
    if c.overlap_div_filter:
        input_output.append(("overlap_div.tsv", os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")))

    jobs = []
    for to_combine, output in input_output:
        if not os.path.exists(output): jobs.append((c.intermediate_dir, c.chunk_size, to_combine, output, c.remove, c.config_header))
    if len(jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(jobs), desc="Combining all results") as pbar:
                for _ in pool.imap_unordered(run_combine_results, jobs):
                    pbar.update()

    # then run the final size and clustering filter
    if c.overlap_div_filter:
        input_path = os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")
    else:
        input_path = os.path.join(c.output_dir, f"{output_name}_first_div_output.tsv")
    medival_final_name = os.path.join(c.output_dir, f"{output_name}_final_regions.tsv")
    medival_summary_final_name = os.path.join(c.output_dir, f"{output_name}_final_regions_summary.tsv")
    if os.path.exists(medival_final_name): print(f"Skipping {medival_final_name}, already exists.")
    if os.path.exists(medival_summary_final_name) and os.path.exists(medival_final_name): print(f"Skipping {medival_final_name}, {medival_summary_final_name}, already exists.")
    else:
        size_filter_cluster(input_path, medival_final_name, medival_summary_final_name, c.size_filter, c.cluster_size, False, c.index, c.config_header)


def write_species_to_file(sequence_id, species, species_file, mode='a'):
    """
    Write a sequence ID and its species to the auto-generated species file.
    
    Args:
        sequence_id: The sequence identifier
        species: The species name
        species_file: Path to the species file
        mode: 'a' for append, 'w' for write
    """
    try:
        with open(species_file, mode) as f:
            f.write(f"{sequence_id}\t{species}\n")
    except Exception as e:
        print(f"Error writing to species file {species_file}: {e}")

def detect_species(record, config):
    with tempfile.NamedTemporaryFile(prefix=f"{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
        SeqIO.write(record, tmp_fasta, "fasta")
        tmp_fasta_path = tmp_fasta.name
    species = get_q_species(tmp_fasta_path, config.kraken_db)
    os.remove(tmp_fasta_path)
    return species

def process_single_record(record, config):
    chunk_jobs = []
    tmp_fastas = []
    # Get the length of the sequence
    seq_len = len(record.seq)
    record.id = record.id.replace("'", "").replace('"', "")
    # Get the species of the sequence
    if config.species_mode == "auto":
        try:
            species = detect_species(record, config)
            #print(f"Species of {record.id}: {species}")
            write_species_to_file(record.id, species, config.speciesFile, mode='a')
        except Exception as e:
            print(f"Issue extracting species for {record.id}: {e}")
            species = "unclassified"
            write_species_to_file(record.id, species, config.speciesFile, mode='a')
    elif config.species_mode == 'file':
        if record.id in config.seq_species_dict:
            species = config.seq_species_dict[record.id]
        else:
            try:
                print(f"Warning: {record.id} not found in species file {config.speciesFile}. Auto-detecting using Kraken...")
                species = detect_species(record, config)
                #print(f"Species of {record.id}: {species}")
                write_species_to_file(record.id, species, config.speciesFile, mode='a')
            except Exception as e:
                print(f"Issue auto-detecting species for {record.id}: {e}")
                species = "unclassified"
                write_species_to_file(record.id, species, config.speciesFile, mode='a')
            
    else:
        species = config.species
    #print(f"Species of {record.id}: {species}")
    
    # If the sequence length is greater than the chunk size we want
    if seq_len > config.chunk_size:
        for i in range(0, seq_len, config.chunk_size):
            # Then we chunk it by size
            chunk_seq = record.seq[i:i+config.chunk_size]
            chunk_record = record[:0]  # copy header
            # Set the sequence to be this chunk
            chunk_record.seq = chunk_seq
            # Set the id to be the chunk number
            idx = i // config.chunk_size
            chunk_record.id = f"{record.id}"
            chunk_record.description = f"{record.description} chunk {idx}"
            
            # Get a tmp file of the sequence
            with tempfile.NamedTemporaryFile(prefix=f"chunk_{idx}_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                SeqIO.write(chunk_record, tmp_fasta, "fasta")
                tmp_fasta_path = tmp_fasta.name
                tmp_fastas.append(tmp_fasta_path)
            
            # Record the actual sequence and its record, its chunk number, and its id
            chunk_jobs.append((chunk_record, idx, record.id, species, tmp_fasta_path))
    else:
        # Get a tmp file of the sequence
        with tempfile.NamedTemporaryFile(prefix=f"chunk_0_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
            SeqIO.write(record, tmp_fasta, "fasta")
            tmp_fasta_path = tmp_fasta.name
            tmp_fastas.append(tmp_fasta_path)
        
        # If it is smaller than the chunk size, just add it
        chunk_jobs.append((record, 0, record.id, species, tmp_fasta_path))
    return chunk_jobs, tmp_fastas

def prepare_jobs_parallel(c: Config, n_processes=None):
    if n_processes is None:
        n_processes = min(c.max_threads, mp.cpu_count())
    
    # If auto mode, initialize the species file (clear it if it exists)
    if c.species_mode == 'auto':
        try:
            # Create parent directories if needed
            Path(c.speciesFile).parent.mkdir(parents=True, exist_ok=True)
            print(f"Auto-detecting species and saving to: {c.speciesFile}")
        except Exception as e:
            print(f"Error initializing species file {c.speciesFile}: {e}")
            return [], []
    elif c.species_mode == 'file':
        # Load the species dictionary from the file if in file mode
        try:
            c.seq_species_dict = species_seq_dict(c.speciesFile)
            print(f"Loaded species assignments from: {c.speciesFile}, Length: {len(c.seq_species_dict)}")
        except Exception as e:
            print(f"Error loading species file {c.speciesFile}: {e}")
            c.seq_species_dict = {}

    # Read all records first
    records = list(SeqIO.parse(c.input_fasta, "fasta"))
    print(f"Processing {len(records)} sequence(s) using {n_processes} threads...")
    
    # Create a partial function with the config
    process_func = partial(process_single_record, config=c)
    
    # Combine all results
    all_chunk_jobs = []
    all_tmp_fastas = []

    lock = threading.Lock()  # Protect writes to shared lists
    def worker_wrapper(record):
        return process_func(record)
    with ThreadPoolExecutor(max_workers=n_processes) as executor:
        with tqdm(total=len(records), desc="Processing sequence(s)") as pbar:
            for result in executor.map(worker_wrapper, records):
                chunk_jobs, tmp_fastas = result
                with lock:
                    all_chunk_jobs.extend(chunk_jobs)
                    all_tmp_fastas.extend(tmp_fastas)
                pbar.update()
 
        # Print species mode summary
    if c.species_mode == 'auto':
        print(f"Species assignments saved to: {c.speciesFile}")
    elif c.species_mode == 'single':
        print(f"Using single species for all sequences: {c.species}")
    else:
        print(f"Using species assignments from: {c.speciesFile} (updated with auto-detected entries)")
    
    # Run skani search to find all reference sequences with >= 95% ANI to each query
    skani_pkl = os.path.join(c.output_dir, f"skani_ani_dict_{c.input_fasta.split("/")[-1].split(".")[0]}.pkl")

    if os.path.exists(skani_pkl):
        print(f"Skipping skani search: {skani_pkl} already exists. Loading...")
        with open(skani_pkl, "rb") as f:
            skani_ani_dict = pickle.load(f)
    else:
        skani_tsv = os.path.join(c.output_dir, f"skani_search_{c.input_fasta.split("/")[-1].split(".")[0]}_95ani.tsv")

        if not os.path.exists(skani_tsv):
            print("\nRunning skani search to find sequences with >= 95% ANI matches...")
            try:
                result = subprocess.run(
                    [
                        "skani", "search",
                        "--qi", str(c.input_fasta),
                        "-t", str(n_processes),
                        "-s", "90",
                        "-d", str(c.skani_sketch_db),
                        "-o", skani_tsv
                    ],
                    check=True,
                    capture_output=True,
                    text=True
                )
                if result.stderr:
                    print(result.stderr.strip())
            except subprocess.CalledProcessError as e:
                print(f"ERROR running skani search: {e.stderr}")
                raise
            except FileNotFoundError:
                print("ERROR: skani not found in PATH.")
                raise
        else:
            print(f"Skipping skani search: {skani_tsv} already exists.")

        # Build dict: query_id -> [ref_ids with ANI >= 95]
        skani_ani_dict = {}
        with open(skani_tsv, "r") as f:
            header = f.readline().strip().split("\t")
            try:
                ani_col = header.index("ANI")
                query_col = header.index("Query_name")
                ref_col = header.index("Ref_name")
            except ValueError:
                ani_col, query_col, ref_col = 2, 6, 5

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) <= max(ani_col, query_col, ref_col):
                    continue
                try:
                    ani = float(parts[ani_col])
                except ValueError:
                    continue
                if ani < 95.0:
                    continue
                query_id = parts[query_col].split()[0]
                ref_id = parts[ref_col].split()[0]
                if query_id not in skani_ani_dict:
                    skani_ani_dict[query_id] = []
                skani_ani_dict[query_id].append(ref_id)

        with open(skani_pkl, "wb") as f:
            pickle.dump(skani_ani_dict, f)
        print(f"  Saved ANI dict to {skani_pkl}")

    print(f"  Found ANI >= 95% matches for {len(skani_ani_dict):,} query sequence(s).")
    c.skani_ani_dict = skani_ani_dict
    return all_chunk_jobs, all_tmp_fastas

def main():
    args = parse_args()

    config_data = {}
    if args.config:
        config_data = load_config_file(args.config)
        if config_data:  # Only validate if config was successfully loaded
            validate_config_keys(config_data)
            print(f"Loaded and validated configuration from: {args.config}")
        else:
            print(f"Warning: Config file {args.config} was empty or could not be loaded")
    
    # Merge config file and command line arguments
    merged_config = merge_config_and_args(args, config_data)

    # Validate required parameters
    if not validate_required_params(merged_config):
        return

    c = Config(
        input_fasta=merged_config['query'],
        chunk_size=merged_config['chunk'],
        output_dir=merged_config['output'],
        blat_db=Path(f"{merged_config['database']}/blat_2bit_db"),
        skani_sketch_db=Path(f"{merged_config['database']}/skani_sketch_db"),
        index=merged_config['index'],
        max_threads=merged_config['threads'],
        kraken_db=merged_config['kraken'],
        tree=merged_config['tree'],
        remove=merged_config['remove'],
        species=merged_config['species'],
        minScore=merged_config['minScore'],
        minIdentity=merged_config['minIdentity'],
        overlap_div_filter=merged_config['overlap_div_filter'],
        speciesFile = merged_config["speciesFile"],
        species_mode = merged_config["species_mode"],
        seq_species_dict = merged_config["seq_species_dict"],
        size_filter = merged_config["size_filter"],
        cluster_size = merged_config["cluster_size"],
        intermediate_dir=os.path.join(merged_config['output'], f"intermediate_{merged_config['output']}_files")
    )
    c.config_header = format_config_header(merged_config)

    triangle_pkl = Path(f"{merged_config['database']}/skani_triangle_ani95.pkl")
    if triangle_pkl.exists():
        with open(triangle_pkl, "rb") as f:
            c.skani_triangle_dict = pickle.load(f)
        print(f"Loaded skani triangle dict: {len(c.skani_triangle_dict):,} pairs")
    else:
        print(f"Warning: skani triangle pkl not found at {triangle_pkl}. Overlap-div ANI lookup will fall back to crude ANI.")

    print("Configuration:")
    print(f"  Input FASTA: {c.input_fasta}")
    print(f"  Output Directory: {c.output_dir}")
    print(f"  Blat Database: {c.blat_db}")
    print(f"  Skani Sketch Database: {c.skani_sketch_db}")
    print(f"  Skani Triangle PKL: {triangle_pkl}")
    print(f"  Index: {c.index}")
    print(f"  Threads: {c.max_threads}")
    print(f"  Chunk Size: {c.chunk_size}")
    print(f"  Species: {c.species or 'Auto-detect with Kraken'}")
    print(f"  Min Score: {c.minScore}")
    print(f"  Min Identity: {c.minIdentity}")
    print(f"  Overlap Div Filter: {c.overlap_div_filter}")
    print(f"  Size Filter: {c.size_filter} bp")
    print(f"  Cluster Size: {c.cluster_size} bp")
    print(f"  Clean up files: {c.remove}")
    print(f"  Species detect mode: {c.species_mode}")
    print(f"  Species file: {c.speciesFile}\n")

    if not os.path.exists(c.output_dir):
        os.makedirs(c.output_dir)
    
    if not os.path.exists(c.intermediate_dir):
        os.makedirs(c.intermediate_dir)

    print("Start time:", datetime.now())

    #check and ensure the database is sound
    blat_files = [Path(f) for f in os.listdir(c.blat_db) if f.endswith(".2bit")]
    blat_files.sort()
    ooc_files = [Path(f) for f in os.listdir(c.blat_db) if f.endswith(".ooc")]
    ooc_files.sort()
    if len(blat_files) != len(ooc_files):
        raise ValueError(
            f"Mismatch in file counts:\n"
            f"- Found {len(blat_files)} blat file(s)\n"
            f"- Found {len(ooc_files)} OOC file(s)\n"
            "Please ensure every blat file has a matching OOC file."
        )
    
    # go through and build jobs
    chunk_jobs, tmp_fastas = prepare_jobs_parallel(c)
    
    print(f"Prepared {len(chunk_jobs)} jobs(s) to process.")

    #this will parallelize it
    #run the function run_blat_on_chunk on all different threads
    print(f"Running on {c.max_threads} threads")
    run_blat_on_chunk(chunk_jobs, blat_files, ooc_files, c)
    print("All blat jobs finished.")

    print("Now filtering blat Output")
    run_all_filters(chunk_jobs, c)
    print("All filtering complete")

    print("Creating final result files")
    combine_all_results(c)

    #remove all of the tmp fastas that were created
    for tmp in tmp_fastas:
        os.remove(tmp)
        
    print("MEDIVAL FINISHED")

    print("Endtime time:", datetime.now())
    
main()