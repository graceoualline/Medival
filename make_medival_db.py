# this code will intake a gtdb multifasta file or fasta file of your choice
# and give you a fasta database that you can use for blat
from Bio import SeqIO
import sys
import os
import subprocess
import time
from pathlib import Path
import shutil
import pathlib

import pickle
import tempfile

def write_group(group_seqs, index, output_dir):
        # if somehow the list is empty, skip
        if not group_seqs:
            return
        out_path = os.path.join(output_dir, f"split_{index}.fasta")
        with open(out_path, "w") as handle:
            SeqIO.write(group_seqs, handle, "fasta")
        print(f"Wrote {out_path}")

def progress(processed_seqs, total_seqs, start_time, file_index):
    # Progress indicator
    if processed_seqs % 10000 == 0 or processed_seqs == total_seqs:
        elapsed = time.time() - start_time
        percent = (processed_seqs / total_seqs) * 100
        rate = processed_seqs / elapsed if elapsed > 0 else 0
        eta = (total_seqs - processed_seqs) / rate if rate > 0 else 0
        
        print(f"Progress: {processed_seqs:,}/{total_seqs:,} sequences ({percent:.1f}%) | "
                f"Rate: {rate:.0f} seq/s | ETA: {eta:.0f}s | Files created: {file_index-1}")

def count_sequences(input_fasta):
    print(f"Counting sequences in {input_fasta}...")
    count = 0
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
                if count % 10000 == 0: print("Current count:", count)
    return count

def count_sequences_in_directory(directory):
    """Count total sequences in all FASTA files in a directory"""
    if not os.path.exists(directory):
        return 0
    
    total_count = 0
    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
    
    for fasta_file in fasta_files:
        fa = os.path.join(directory, fasta_file)
        total_count += count_sequences(fa)
    return total_count

def check_step1_completion(input_fasta, output_dir, original_count):
    """Check if step 1 is complete by comparing sequence counts"""
    if not os.path.exists(output_dir):
        return False
    
    # Count sequences in split files
    split_count = count_sequences_in_directory(output_dir)
    
    is_complete = original_count == split_count and split_count > 0
    
    print(f"Step 1 status check:")
    print(f"  Original file sequences: {original_count:,}")
    print(f"  Split files sequences: {split_count:,}")
    print(f"  Step 1 complete: {'Yes' if is_complete else 'No'}")
    
    return is_complete

def split_fasta_by_bp(input_fasta, blat_fasta_db, max_bp, total_seqs):
    output_dir = str(blat_fasta_db)
    seq_lengths = dict()
    # Count total sequences for progress tracking

    group = []
    group_bp = 0
    file_index = 1
    processed_seqs = 0
    start_time = time.time()

    print("\n" + "="*60)
    print("STEP 1: Splitting FASTA file by base pairs")
    print("="*60)
    complete = check_step1_completion(input_fasta, output_dir, total_seqs)
    if complete: 
        print(f"Skipping splitting fastas, step already completed.")
        return None

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)

        processed_seqs += 1
        progress(processed_seqs, total_seqs, start_time, file_index)
        seq_lengths[record.id] = seq_len
        
        # If one sequence is longer than max_bp, write it by itself
        if seq_len > max_bp:
            print(f"Warning: {record.id} is {seq_len} bp, longer than max ({max_bp} bp). Saving in its own file.")
            write_group([record], file_index, output_dir)
            file_index += 1
            continue

        # If adding this record would exceed max_bp, flush and start new group
        if group_bp + seq_len > max_bp:
            write_group(group, file_index, output_dir)
            file_index += 1
            group = [record]
            group_bp = seq_len
        else:
            group.append(record)
            group_bp += seq_len

    # Write the final group
    if group:
        write_group(group, file_index, output_dir)

    total_time = time.time() - start_time
    print(f"\nStep 1 completed in {total_time:.1f}s")
    print(f"Created {file_index} FASTA files")
    return seq_lengths

def make_skani_sketch_db(input_fasta, sketch_dir, threads):
    """
    Create a pre-sketched skani database from a multi-fasta file.
    Each sequence is written to its own temporary FASTA file so that skani
    keeps every sequence as a separate sketch entry.  A file-list is passed
    to `skani sketch -l` to avoid hitting shell argument-length limits.

    This sketch database is used later by `skani search` during the main
    medival pipeline (filter_blat.py) to quickly compute ANI between a
    query chunk and all database sequences in a single pass.

    Args:
        input_fasta: Path to the (multi-)fasta file containing all DB sequences.
        sketch_dir:  Output directory where .sketch files will be written.
        threads:     Number of threads to pass to skani.
    """
    print("\n" + "="*60)
    print("STEP 4: Creating skani sketch database")
    print("="*60)

    sketch_dir = Path(sketch_dir)

    if sketch_dir.is_dir():
        print(f"Skipping: {sketch_dir} already exists, either database is already created, or skani cannot write over preexisting directories")
        print(f"If sketched directory was not created, please delete {sketch_dir} and try again.")
        return

    print(f"Splitting sequences from {input_fasta} into individual temp files...")
    print(f"  Using {threads} thread(s)")

    start_time = time.time()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        list_file = tmp_path / "seq_file_list.txt"

        seq_count = 0
        with open(list_file, "w") as lf:
            for i, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
                safe_id = record.id.replace("/", "_").replace("\\", "_")
                seq_file = tmp_path / f"{i}_{safe_id}.fa"
                with open(seq_file, "w") as sf:
                    SeqIO.write(record, sf, "fasta")
                lf.write(f"{seq_file}\n")
                seq_count += 1
                if seq_count % 1000 == 0:
                    elapsed = time.time() - start_time
                    print(f"  Written {seq_count:,} temp files ({elapsed:.1f}s elapsed)...")

        print(f"  Wrote {seq_count:,} individual sequence files")
        print(f"  Running skani sketch...")

        try:
            result = subprocess.run(
                [
                    "skani", "sketch",
                    "-t", str(threads),
                    "-l", str(list_file),
                    "-o", str(sketch_dir)
                ],
                check=True,
                capture_output=True,
                text=True
            )
            if result.stderr:
                print(result.stderr.strip())
        except subprocess.CalledProcessError as e:
            print(f"ERROR creating skani sketch database:")
            print(f"  stderr: {e.stderr}")
            raise
        except FileNotFoundError:
            print("ERROR: skani command not found. Please ensure it's installed and in your PATH.")
            raise

    elapsed = time.time() - start_time
    print(f"Done! Created {sketch_dir} ({elapsed:.1f}s)")

def make_skani_triangle(input_fasta, output_dir, threads):
    """
    Pre-compute pairwise ANI for all database sequences using skani triangle.
 
    This produces two output files inside output_dir:
      1.  skani_triangle_ani95.tsv   — human-readable filtered edge list
      2.  skani_triangle_ani95.pkl   — pickle dict for fast O(1) lookup
 
    The pickle dict uses normalised keys:  (min(id1, id2), max(id1, id2)) → ANI
    so that order doesn't matter.  If a pair is *not* in the dict it means
    ANI < 95 %  (or the pair was never compared because ANI < 94 %, skani's
    pre-filter threshold set by -s 94).
 
    This step can take a very long time on large databases
    (e.g. ~14 hours on the full GTDB with 11 M sequences).  The resulting
    lookup table is used by the overlap filters (find_overlap.py,
    find_overlap_and_div.py) to replace the old per-pair skani dist calls.
 
    Args:
        input_fasta: Path to the (multi-)fasta with all DB sequences.
        output_dir:  Directory where output files are written.
        threads:     Number of threads.
    """
    print("\n" + "="*60)
    print("STEP 5: Pre-computing all-vs-all ANI (skani triangle)")
    print("="*60)
 
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
 
    raw_output      = output_dir / "skani_triangle_raw.tsv"
    filtered_tsv    = output_dir / "skani_triangle_ani95.tsv"
    pkl_output      = output_dir / "skani_triangle_ani95.pkl"

     # ---- resume: if the pickle already exists we're done ----
    if pkl_output.exists():
        print(f"Skipping: {pkl_output} already exists")
        return
    # ---- 5a. run skani triangle ----
    if not raw_output.exists():
        print(f"Running skani triangle on {input_fasta}...")
        print(f"  Threads: {threads}")
        print(f"  Pre-filter threshold: 90 % ANI  (-s 90)")
        print("  WARNING: This could take a very long time.")
        print("           (≈14 hours for the full GTDB of 11 M sequences if -s=94" \
        "                  about 7 days if -s = 85))")
 
        start_time = time.time()
        try:
            # Redirect stdout to the raw output file instead of using shell '>'
            with open(raw_output, "w") as outfile:
                result = subprocess.run(
                    [
                        "skani", "triangle",
                        "-t", str(threads),
                        "-s", "90",
                        "-i", str(input_fasta),
                        "-E"
                    ],
                    stdout=outfile,
                    text=True,
                    check=True
                )
                if result.stderr:
                    # skani prints progress info to stderr; echo it
                    print(result.stderr.strip())
        except subprocess.CalledProcessError as e:
            print(f"ERROR running skani triangle:")
            print(f"  stderr: {e.stderr}")
            # Clean up partial output so resume doesn't skip it
            if raw_output.exists():
                raw_output.unlink()
            raise
        except FileNotFoundError:
            print("ERROR: skani command not found. Please ensure it's installed and in your PATH.")
            raise
 
        elapsed = time.time() - start_time
        print(f"  skani triangle finished in {elapsed:.1f}s")
    else:
        print(f"Skipping skani triangle run: {raw_output} already exists")
 
    # ---- 5b. filter to ANI >= 95 % and build lookup dict ----
    print("Processing raw triangle output...")
    print("  Filtering to ANI >= 95 percent and building lookup table...")
 
    ani_lookup = {}   # (min_id, max_id) → best ANI
    total_lines = 0
    kept_lines = 0
    start_time = time.time()
 
    with open(raw_output, "r") as infile, \
         open(filtered_tsv, "w") as outfile:
 
        outfile.write("Seq_ID_1\tSeq_ID_2\tANI\n")
 
        # Skip the header line
        header = infile.readline()
 
        for line in infile:
            total_lines += 1
 
            # Progress every 10 M lines
            if total_lines % 10_000_000 == 0:
                elapsed = time.time() - start_time
                rate = total_lines / elapsed if elapsed > 0 else 0
                print(f"  Processed {total_lines:,} lines  "
                      f"(kept {kept_lines:,})  "
                      f"| {rate:,.0f} lines/s")
 
            parts = line.split("\t")
            if len(parts) < 3:
                continue
 
            try:
                ani = float(parts[2])
            except ValueError:
                continue
 
            if ani < 95.0:
                continue

            id1 = parts[5].split()[0]
            id2 = parts[6].split()[0]
 
            # skip self-comparisons
            if id1 == id2:
                continue
 
            # Normalise key so lookup(a,b) == lookup(b,a)
            key = (min(id1, id2), max(id1, id2))
 
            # Keep the best ANI if a pair appears more than once
            if key not in ani_lookup or ani > ani_lookup[key]:
                ani_lookup[key] = ani
 
            outfile.write(f"{id1}\t{id2}\t{ani}\n")
            kept_lines += 1
 
    elapsed = time.time() - start_time
    print(f"  Finished processing in {elapsed:.1f}s")
    print(f"  Total lines read:   {total_lines:,}")
    print(f"  Pairs kept (>=95%): {kept_lines:,}")
    print(f"  Unique pairs in lookup: {len(ani_lookup):,}")
 
    # ---- 5c. save the lookup dict ----
    print(f"  Saving lookup table to {pkl_output}...")
    with open(pkl_output, "wb") as f:
        pickle.dump(ani_lookup, f)
 
    file_size_mb = pkl_output.stat().st_size / (1024 * 1024)
    print(f"  Pickle size: {file_size_mb:.1f} MB")
    print(f"Done! Triangle outputs in {output_dir}:")
    print(f"  Filtered TSV:  {filtered_tsv}")
    print(f"  Lookup pickle: {pkl_output}")

def lookup_triangle_ani(triangle_dict, id1, id2):
    """
    Look up the ANI between two reference sequences in the triangle dict.
 
    Returns:
        float  — the ANI value if the pair was found (>= 95 %).
        None   — if the pair is not in the table, meaning ANI < 95 %
                 (or < 95 %, skani's pre-filter).
    """
    key = (min(id1, id2), max(id1, id2))
    return triangle_dict.get(key, None)

def convert_2bit(input_directory, output_directory):
    # Get list of files to process
    fasta_files = [f for f in os.listdir(input_directory) 
                   if os.path.isfile(os.path.join(input_directory, f))]
    total_files = len(fasta_files)

    # check if this has already been done
    twobit_files = list(output_directory.glob("*.2bit"))
    
    if len(twobit_files) == total_files:
        print(f"Skipping blat convert fasta to 2bit: {len(twobit_files)} files already exist in {output_directory}.")
        return

    print("\n" + "="*60)
    print("STEP 2: Converting FASTA files to 2bit format")
    print("="*60)
    print(f"Total files to convert: {total_files}")

    file_type = "2bit"
    start_time = time.time()

    for i, input_file in enumerate(fasta_files, 1):
        file_name, file_extension = os.path.splitext(input_file)
        
        # Call your Python script with the input and output file names
        input_path = os.path.join(input_directory, input_file)
        output_path = os.path.join(output_directory, f"{file_name}.{file_type}")
        # if the output file already exists, we don't need to make it again
        if os.path.exists(output_path): continue
        # Progress indicator
        percent = (i / total_files) * 100
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_files - i) / rate if rate > 0 else 0
        
        print(f"Progress: {i}/{total_files} files ({percent:.1f}%) | "
              f"Rate: {rate:.1f} files/s | ETA: {eta:.0f}s")

        try:
            result = subprocess.run(["faToTwoBit", input_path, output_path], 
                                  capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting {input_file}: {e}")
            print(f"stderr: {e.stderr}")
        except FileNotFoundError:
            print("ERROR: faToTwoBit command not found. Please ensure it's installed and in your PATH.")
            break

    total_time = time.time() - start_time
    print(f"\nStep 2 completed in {total_time:.1f}s")
    print(f"Converted {total_files} files to 2bit format")

def make_ooc(input_directory):
    bit_files = [f for f in os.listdir(input_directory) 
                 if f.endswith('.2bit') and os.path.isfile(os.path.join(input_directory, f))]

    total_files = len(bit_files)

    # check if this has already been done
    ooc_files = list(input_directory.glob("*.ooc"))
    
    if len(ooc_files) == total_files:
        print(f"Skipping making ooc files: {len(ooc_files)} files already exist in {input_directory}.")
        return
    
    print("\n" + "="*60)
    print("STEP 3: Generating OOC files from 2bit files")
    print("="*60)
    print(f"Total 2bit files to process: {total_files}")

    if total_files == 0:
        print("ERROR: No 2bit files found in directory:", input_directory)
        return

    start_time = time.time()

    for i, bit_file in enumerate(bit_files, 1):
        base_name = bit_file.replace('.2bit', '')
        # Check if the file is a regular file
        input_path = os.path.join(input_directory, bit_file)
        output_file = f"{base_name}.ooc"
        output_path = os.path.join(input_directory, output_file)
        # if the output file already exists, we don't need to make it again
        if os.path.exists(output_path): continue

        percent = (i / total_files) * 100
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_files - i) / rate if rate > 0 else 0
        print(f"Progress: {i}/{total_files} files ({percent:.1f}%) | "
              f"Rate: {rate:.1f} files/s | ETA: {eta:.0f}s")

        try:
            cmd = ["blat", input_path, "/dev/null", "/dev/null", f"-makeOoc={output_path}", "-repMatch=1024"]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR generating OOC for {bit_file}: {e}")
            print(f"stderr: {e.stderr}")
        except FileNotFoundError:
            print("ERROR: blat command not found. Please ensure it's installed and in your PATH.")
            break
    total_time = time.time() - start_time
    print(f"\nStep 3 completed in {total_time:.1f}s")
    print(f"Generated OOC files for {total_files} 2bit files")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python make_medival_db.py <input.fasta> <output_dir> <max_bil_bp> <threads>")
        sys.exit(1)
    print("Started")
    input_fasta = sys.argv[1] #full_database = "/usr1/shared/gtdb_combined.fa"
    output_dir = Path(sys.argv[2])

    blat_fasta_split_dir = output_dir / "blat_fasta_db"
    blat_2bit_split_dir = output_dir / "blat_2bit_db"
    skani_dir = output_dir / "skani_sketch_db"

    for d in [output_dir, blat_fasta_split_dir, blat_2bit_split_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
    
    max_bp = int(float(sys.argv[3]) * 1000000000)
    threads = int(sys.argv[4])
    og_dir = os.getcwd()
    print("Loaded")
    seq_lengths = Path(f"{output_dir}/seq_lengths.tsv")
    print(seq_lengths)
    if os.path.exists(seq_lengths):
        print("Skipping obtaining sequence lengths, seq_lengths.tsv already exists")
    else:
        total_seqs = count_sequences(input_fasta)
        seqs_lengths = split_fasta_by_bp(input_fasta, blat_fasta_split_dir, max_bp, total_seqs)
        if seqs_lengths: 
            os.chdir(output_dir)
            with open("seq_lengths.tsv", "w") as f:
                for s in seqs_lengths:
                    f.write(f"{s}\t{seqs_lengths[s]}\n")
            os.chdir(og_dir)


    print("\nFinished splitting the query file. Converting to 2bit.")

    convert_2bit(blat_fasta_split_dir, blat_2bit_split_dir)

    print("\nFinished 2bit conversion. Generating OOC files.")

    make_ooc(blat_2bit_split_dir)

    print("\n Finished making ooc files, now making skani database")
    
    make_skani_sketch_db(input_fasta, Path(skani_dir), threads)

    print("\n Finished making skani database, now calculating skani triangle associations")
    
    make_skani_triangle(input_fasta, output_dir, threads)
        
    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)



# python3 /usr1/gouallin/blat/blat_pipeline/make_medival_db.py /usr1/shared/gtdb_combined.fa gtdb_smart_split_fa 2
