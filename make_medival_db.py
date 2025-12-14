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
    for _ in SeqIO.parse(input_fasta, "fasta"):
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

def make_skani_db(split_fasta_dir, skani_db_dir, total_seqs):
    """
    Split FASTA files into individual sequences for skani database.
    
    Args:
        split_fasta_dir: Directory containing FASTA files to process
        output_dir: Base output directory where medival_skani_db will be created
    """
    #split_fasta_dir = Path(split_fasta_dir)
    
    # Create the skani database directory
    #skani_db_dir = output_dir / "medival_skani_db"
    print("\n" + "="*60)
    print("STEP 4: Make Skani DB")
    print("="*60)
    # Get all FASTA files in the input directory

    fasta_files = [f.resolve() for f in split_fasta_dir.glob("*.fasta")] + \
              [f.resolve() for f in split_fasta_dir.glob("*.fa")] + \
              [f.resolve() for f in split_fasta_dir.glob("*.fna")]
    
    check = len(list(skani_db_dir.glob("*")))


    
    if check == total_seqs and total_seqs > 0:
        print(f"Skipping skani db build: {check} sequences already exist in {skani_db_dir}.")
        return
    
    total_files = len(fasta_files)
    print(f"Processing {total_files} FASTA files...")
    
    # AWK command to split FASTA into individual sequence files
    awk_script = r'''
    /^>/{
        split($0, a, " ");
        id = a[1];
        gsub("^>", "", id);
        filename = id ".fasta";
    }
    { print > filename }
    '''
    
    # Process each FASTA file
    start_time = time.time()
    og_dir = os.getcwd()
    os.chdir(skani_db_dir)
    for i, fasta_file in enumerate(fasta_files):
        # Change to the skani_db directory so files are created there
        percent = (i / total_files) * 100
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_files - i) / rate if rate > 0 else 0
        print(f"Progress: {i}/{total_files} files ({percent:.1f}%) | "
              f"Rate: {rate:.1f} files/s | ETA: {eta:.0f}s")

        try:
            # Run awk command on the FASTA file
            subprocess.run(
                ['awk', awk_script, str(fasta_file)],
                check=True,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error processing {fasta_file}: {e.stderr}")
    os.chdir(og_dir)
    print(f"Done! Individual sequences saved to {skani_db_dir}")

def convert_2bit(input_directory, output_directory):
    # input_directory = output_name + "_fa"
    # output_directory = output_name + "_2bit"
    # if not os.path.exists(output_directory):
    #         os.makedirs(output_directory)

    # Get list of files to process
    fasta_files = [f for f in os.listdir(input_directory) 
                   if os.path.isfile(os.path.join(input_directory, f))]
    total_files = len(fasta_files)

    # check if this has already been done
    twobit_files = list(output_directory.glob("*.2bit"))
    
    if len(twobit_files) == total_files:
        print(f"Skipping blat convert fasta to 2bit: {len(twobit_files)} sequences already exist in {output_dir}.")
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
    #input_directory = output_name + "_2bit"
    bit_files = [f for f in os.listdir(input_directory) 
                 if f.endswith('.2bit') and os.path.isfile(os.path.join(input_directory, f))]

    total_files = len(bit_files)

    # check if this has already been done
    ooc_files = list(input_directory.glob("*.ooc"))
    
    if len(ooc_files) == total_files:
        print(f"Skipping making ooc files: {len(ooc_files)} sequences already exist in {input_directory}.")
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
        # Generate output filename (replace .2bit with .ooc)
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
    if len(sys.argv) != 4:
        print("Usage: python make_medival_db.py <input.fasta> <output_dir> <max_bil_bp>")
        sys.exit(1)
    
    input_fasta = sys.argv[1] #full_database = "/usr1/shared/gtdb_combined.fa"
    output_dir = Path(sys.argv[2])

    blat_fasta_split_dir = output_dir / "blat_fasta_db"
    blat_2bit_split_dir = output_dir / "blat_2bit_db"
    skani_dir = output_dir / "skani_db"

    for d in [output_dir, blat_fasta_split_dir, skani_dir, blat_2bit_split_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
    
    max_bp = int(sys.argv[3]) * 1000000000

    total_seqs = count_sequences(input_fasta)
    og_dir = os.getcwd()

    
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
    
    make_skani_db(Path(blat_fasta_split_dir), Path(skani_dir), total_seqs)
        
    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)
    

# python3 /usr1/gouallin/blat/blat_pipeline/make_medival_db.py /usr1/shared/gtdb_combined.fa gtdb_smart_split_fa 2
