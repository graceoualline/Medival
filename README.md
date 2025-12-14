# medival

**M**obile (Genetic) **E**lement finder using **DIV**ergence and **AL**ignment. A parallelized sequence alignment tool for detecting Mobile Genetic Elements (MGEs) using BLAT and phylogenetic divergence analysis.
## Overview

medival is designed to find novel MGEs through optimized threading of BLAT (BLAST-Like Alignment Tool) operations, and advanced filtering. It splits large sequences into manageable chunks, runs BLAT searches in parallel across multiple database files, and applies divergence filtering and sequence analysis to identify potential MGEs. The tool is particularly effective for horizontal gene transfer detection.

## Installation
Install medival:
```bash
git clone https://github.com/graceoualline/medival.git
cd medival
```
Install supporting files and databases here: 
TODO: please upload all files that are in /usr1/shared/medival_github_large_files.
You should have the following files:
```
├── divergence_tree.tar.gz
├── kraken2_db.tar.gz
└── medival_gtdb_db
    ├── blat_2bit_db -> /usr1/shared/medival_gtdb_db
    ├── medival_db_index.tar.gz
    └── skani_db.tar.gz
# you can also download these if you want to build the db/index yourself:
gtdb_all_seqs.tar.gz
all_gtdb_seq_kraken_species.tar.gz
```
Decompress the files:
```
tar -xzvf divergence_tree.tar.gz
tar -xzvf kraken2_db.tar.gz
tar -xzvf medival_gtdb_db
tar -xzvf blat_2bit_db -> /usr1/shared/medival_gtdb_db
tar -xzvf medival_db_index.tar.gz
tar -xzvf skani_db.tar.gz
```
After decompressing the files, you should have the following:
```
divergence_tree/
├── TimeTree_v5_Final.nwk
├── TimeTree_v5.hashtable.pkl
├── TimeTree_v5.index.npy
├── TimeTree_v5.mins.npy
└── TimeTree_v5.tour.npy

medival_gtdb_db/
├── blat_2bit_db
│   ├── split_1_output.2bit
│   ├── split_1_output.ooc
...
│   ├── split_137_output.2bit
│   └── split_137_output.ooc
├── medival_db_index.pkl
└── skani_db
    └── (all sequences of the gtdb)

kraken2_db
├── hash.k2d
├── opts.k2d
├── taxo.k2d
└── unmapped.txt
```
### Required Python Packages
```bash
pip install biopython tqdm
```
### Prerequisites
Please ensure you have the following tools downloaded:
- Python 3.7+
- BLAT: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/
  ```bash
  conda install -c bioconda blat
  ```
- Kraken2: https://github.com/DerrickWood/kraken2/wiki/Manual
  ```bash
  conda install -c bioconda kraken2
  ```

## Usage

### Quick Start

```bash
# To see all input parameters
python3 medival.py -h

# command line with only required arguments
python3 medival.py -q input.fasta -o output_directory -d /path/to/medival_database -tr /path/to/divergence_tree -k /path/to/kraken_db --threads 20

# with config
python medival.py --config config_example.yaml

# with config and command line. Command-line arguments take priority over config file values.
python medival.py --config config_example.yaml -q input.fasta -o output_directory --threads 20
```
#### Ready-to-Run Example (only on FAUST) Will take ~25 minutes.
```bash
# with config
python3  /usr1/gouallin/blat/blat_pipeline/medival.py --config /usr1/gouallin/blat/blat_pipeline/config_example.yaml

# command line only
python3 /usr1/gouallin/blat/blat_pipeline/medival.py \
  -q /usr1/gouallin/blat/blat_pipeline/test/acrB.fasta \
  -o medival_test_results_acrB \
  -d /usr1/shared/medival_gtdb_db/ \
  -tr /usr1/gouallin/blat/divergence_tree \
  -k /usr1/shared/kraken2_custom_db/ \
  -t 20 \
  --no-remove \
  -minIdentity 90 \
  --overlap_filter \
  --overlap_div_filter
```
### Parameters

#### Required Arguments
| Parameter | Description |
|-----------|-------------|
| `-q, --query`| Path to the query FASTA file |
| `-o, --output`| Name of your output directory 
| `-d, --database`| Path to the medival database directory| |
| `-tr, --tree`| Path to the phylogenetic tree directory |
| `-k, --kraken`| Path to Kraken2 database |

#### Optional Arguments
| Parameter | Default | Description |
|-----------|---------|-------------|
|`-t, --threads`| 1 | Number of threads to use - **Highly recommend using more threads to speed up the program.**| 
| `-i, --index` | Taken from prebuilt medival database | Is a pkl formated file that describes all sequences in the Blat DB, its species, length, and leaf name in the phylogenetic tree. If you built your own index with separate species, you can define that index here.
| `-c, --chunk`| 100000| Chunk size for sequence splitting |
| `-s, --species`| auto-detect| Predefined species name (replace spaces with '_'). For multifasta files, applies to all sequences.|
| `-minScore`| 30| Minimum BLAT alignment score |
| `-minIdentity`|0| Minimum percent identity threshold|
| `--remove, --no-remove`|False| Clean up redundant files after running. Use --remove to enable, --no-remove to disable.|
| `--overlap-filter, --no-overlap-filter`|False| Enable overlap filter. Use --overlap-filter to enable, --no-overlap-filter to disable.|
| `--overlap-div-filter, --no-overlap-div-filter`|False| Enable overlap and divergence filter. Use --overlap-div-filter to enable, --no-overlap-div-filter to disable.|

#### Config File
Create a YAML configuration file for repeated analyses:
```yaml
# medival Configuration File
# Required parameters
query: sequences.fasta
output: medival_results
database: medival_gtdb_db/
tree: divergence_tree/
kraken: kraken2_custom_db/

# Optional parameters
threads: 20
chunk: 100000
minIdentity: 95
index: medival_gtdb_db/medival_db_index.pkl
overlap_filter: false
overlap_div_filter: false
remove: false
species: null  # Auto-detect with Kraken2
minScore: 30
minIdentity: 95 
```

### Example Commands

#### Basic Command with GTDB Database and 20 threads
```bash
python medival.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d medival_gtdb_db \
  -tr divergence_tree \
  -k kraken2_custom_db \
  -t 20 \
```

#### Run with additional filtering:
```bash
python medival.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d medival_gtdb_db \
  -tr divergence_tree \
  -k kraken2_custom_db \
  -t 20 \
  --overlap_filter \
  --overlap_div_filter \
  -minIdentity 95
```
## Filters
Filters are described in further detail, and their processes are illustrated in our paper (add cite).
### Divergence Filtering
- Divergence filtering is the main method for detecting MGEs and produces the file ```{output_name}_medival_output.tsv```.
- By examining the species of the query genome and the genome it aligned to, we use the TimeTree of Life to calculate the divergence time between the two species. If the species diverged over 1 million years ago (```divergence >= 1 MYA```), the alignment is retained.
- This filter is effective at identifying horizontal gene transfer events because MGEs transferred between distantly related species will show high sequence similarity despite ancient species divergence.
- For detailed information on how this filter detects MGEs, please refer to our paper: (citation tba).
### Additional filtering
These filters further refine the alignments that were identified as divergently distant in ```{output_name}_medival_output.tsv```:
#### Overlap Filtering
- Enabling ```--overlap_filter``` produces the file ```{output_name}_overlap.tsv```
- This filter processes the results from ```{output_name}_medival_output.tsv```
- It identifies two alignments that overlap spatially on the query sequence but originate from different species
- While effective at reducing false positives, it also reduces medival's sensitivity for detecting true MGEs (see paper for details)
#### Overlap Divergence Filtering
- Enabling ```--overlap_div_filter``` produces the file ```{output_name}_overlap_div.tsv```
- This filter processes the results from {output_name}_medival_output.tsv
-This filter identifies overlapping alignments where the source species are divergently distant (>= 1 MYA) 
- It is the most stringent filter that effectively reduces false positives but may decrease sensitivity for detecting true MGEs and requires longer processing time.
- Recommended for high-confidence MGE detection when processing time is not a constraint

### Standard Output (with remove enabled, and all filtering options on):
```
output_directory/
├── output_name_blat_results.tsv      # Raw BLAT alignments
├── output_name_medival_output.tsv    # Divergence-filtered results
├── output_name_overlap.tsv           # Overlap-filtered results (if enabled)
└── output_name_overlap_div.tsv       # Overlap+divergence filtered (if enabled)
```

#### Example output (with remove disabled, and all filtering options on):
```
medival_test_results_test/
├── chunk_0_test
│   ├── chunk_0_test_part_0.psl
....
│   └── chunk_0_test_part_136.psl
├── chunk_1_test
│   ├── chunk_1_test_part_0.psl
....
│   └── chunk_1_test_part_136.psl
.....
├── chunk_0_test_blat_output.tsv
├── chunk_0_test_medival_output.tsv
├── chunk_0_test_overlap_div.tsv
├── chunk_0_test_overlap.tsv
├── chunk_1_test_blat_output.tsv
├── chunk_1_test_medival_output.tsv
├── chunk_1_test_overlap_div.tsv
├── chunk_1_test_overlap.tsv
.....
├── medival_test_results_test_blat_results.tsv
├── medival_test_results_test_medival_output.tsv
├── medival_test_results_test_overlap_div.tsv
└── medival_test_results_test_overlap.tsv
```
chunk\_{i}\_{output_name} is the directory that will contain all of the raw BLAT output that is run on that chunk of the input sequence. The GTDB-BLAT database is split into 136 pieces to enable efficient parallelization of BLAT. The results are combined, and filtering is done separately for each chunk, generating chunk\_{i}_\{output_name}\_{output or filter type}.tsv for each chunk. Then, at the very end, all chunks are combined into the final output files described below.

### Resume Functionality
Important: The program is designed to resume from interruptions by checking for existing files. If a run is stopped prematurely, it will restart from where it left off. Avoid creating files with names that could overlap with medival's output to prevent conflicts.

## Database Setup

#### Kraken and Kraken Database
We provide the kraken database ```kraken2_custom_db```. If you want to recreate it locally, do the following:
```bash
mkdir -p kraken2_custom_db

kraken2-build --download-taxonomy --db kraken2_custom_db

kraken2-build --download-library archaea --db kraken2_custom_db
kraken2-build --download-library bacteria --db kraken2_custom_db
kraken2-build --download-library plasmid --db kraken2_custom_db
kraken2-build --download-library viral --db kraken2_custom_db
kraken2-build --download-library fungi --db kraken2_custom_db
kraken2-build --download-library protozoa --db kraken2_custom_db
kraken2-build --download-library nt --db kraken2_custom_db

kraken2-build --build --db kraken2_custom_db
kraken2-build --clean --db kraken2_custom_db
```

#### Pre-build GTDB-Blat Database
We provide a ready-to-use BLAT-compatible database created from the Genome Taxonomy Database (GTDB):
- **Database:** ```medival_gtdb_db/blat_2bit_db```
- **Index:** ```medival_gtdb_db/medival_gtdb_2bil_index.pkl```

### Creating a Custom Blat database
-  If you wish to create your own database and index, follow these steps:
#### Step 1: Build the Blat and Skani database
The database must be in 2bit format with a maximum of ~2 billion base pairs per file. We recommend a size of 2 billion base pairs.
```
python3 make_medival_db.py [fasta_file] [output_name] [size_in_bil_bp]
```
For example:
```
python3 make_medival_db.py gtdb.fa medival_gtdb_db 2
```
This script will:
- Generate a file (output_name/seq_lengths.tsv) that contains the lengths of all sequences in your input fasta file
- Split the multifasta file into chunks of specified size
- Convert split files to 2bit format
- Generate .ooc files for all 2bit files
- Ensure each .2bit file has a matching .ooc file
- Split every individual sequence into a separate file for quick use for calculating ani.
#### Step 2: Extract Species Information
Classify sequences using Kraken2:
```
python3 extract_species_from_kraken.py [fasta file] [kraken_db]
```
##### Outputs:
- ```{fasta_file}_kraken_output.txt```: Kraken's detailed classification output
- ```{fasta_file}_seq_species.txt```: Tab-separated file with sequence ID and species:
```
seq_id1  species1
seq_id2  species2
...  ...
```
You can manually edit species assignments in the ```{fasta_file}_seq_species.txt``` file if needed.
#### Step 3: Build the Database Index
Create a hash table index linking sequence IDs to their species, length, and leaf name in the phylogenetic tree:
```
python3 build_database_index.py [medival_database_from_make_medival_db.py] [{fasta_file}_seq_species.txt] [divergence_tree]
```
This creates an index where ```index[seq_id] = (species, length, leaf name in tree)```.
This index is crucial to the pipeline's speed.

## Required Files Summary
1. **BLAT Database:** Directory containing .2bit and .ooc files
2. **Index File:** .pkl file mapping sequence IDs to species, length, and leaf name in the phylogeny tree.
3. **TimeTree Directory:** Phylogenetic tree with relevant species, preprocessed,
4. **Kraken Database:** For taxonomic classification

## Phylogenetic Tree
We use the Time Tree of Life to calculate divergence times between species. From the original .nwk TimeTree file, we preprocessed the tree (as shown in this Google Colab notebook (link)
) to allow finding the closest common ancestor in O(1) time. If a new .nwk file from the Time Tree becomes available, this notebook can be used to generate updated indexes and preprocess the tree for efficient queries. (cite timetree)

## Performance Tips
1. Use Multiple Threads: Set -t to utilize available CPU cores (We use 46 in experimentation, although thread use will depend on one's resources.)
2. Optimal Chunk Size: Default 100kb works well; adjust based on sequence lengths. Increasing chunk size will significantly slow down the program.
3. Enable Filtering: Enabling overlap and overlap divergence filters will decrease false positives, but will also cause medival to miss more MGEs. Overlap divergence takes a very long time to run.
4. Monitor Resources: Large databases require substantial RAM
5. Resume Feature: Take advantage of the resume capability for long runs

## Workflow
**Input Processing:** Reads FASTA sequences and (if not specified by the user) determines species classification via Kraken2

**Chunking:** Splits sequences larger than chunk size into manageable pieces

**Parallel BLAT Search:** Runs BLAT searches against database files using multiple threads

**Filtering Pipeline:** Applies divergence and other specified filters


## Citation

If you use medival in your research, please cite:
[Add later]

## Support

For questions and support, please submit an issue ticket on the GitHub repository.
