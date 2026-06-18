import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import argparse
from datetime import datetime

@dataclass
class Config:
    input_fasta: str #query fasta
    chunk_size: int # how we want to split up the query to run in blat
    output_dir: str # output dir name
    blat_db: Path # blat db path
    skani_sketch_db: Path # skani sketch db path
    index: str #index that maps sequences in the database to their species and locations
    max_threads: int # number of threads the user wants to use
    kraken_db: Path # path to the kraken database
    tree: str # the folder that contains the tree and its preprocessing files
    species: str # user defined species of all of the sequences given in a query
    minScore: int # Is the minimum alignment score threshold, sets blat's parameter -minScore, default 30
    minIdentity: int # Is the threshold of the minimum percent identity a hit must have. Default: 90. Percent identity = ( match / Q_end - Q_start )*100
    speciesFile: str # is a file that defines the species of every sequence in the fasta file
    species_mode: str # code defined mode to help know how species are being defined
    seq_species_dict: str # if file given, the dictionary created relating dict[seq_id]: species
    intermediate_dir: str # the directory that will contain intermediate files that are used to build the final files
    size_filter: int # Filter out any regions smaller than this many bp
    cluster_size: int # Combine any regions within this many bp of each other
    skani_ani_dict: Dict[str, List[str]] = field(default_factory=dict) # query_id -> ref_ids with >= 95% ANI
    skani_triangle_dict: Dict = field(default_factory=dict) # (min_id, max_id) -> ANI for all db pairs >= 95%
    config_header: str = ''

def parse_args():
    parser = argparse.ArgumentParser(
        description="""Take a FASTA file and run it through medival.\n
        
        Config File Usage:\n
    You can create a YAML config file to set default values for all parameters.\n
    Command-line arguments will override config file values.\n
    
    An Example config file is available in config_example.yaml.\n
        
    Example commands:\n
    # Use config file with command-line overrides\n
    python medival.py --config my_config.yaml --threads 46 --chunk 100000\n
    
    # Use only command-line arguments (traditional way)\n
    python medival.py -q input.fasta -o output -d database -tr tree -i index -k kraken"""
    )
    parser.add_argument(
        "--config", type=str, help="Path to a YAML config file with default arguments.")
    parser.add_argument(
        "-q", "--query", help="Path to the query FASTA file")
    parser.add_argument(
        "-o", "--output", help="Name of your output directory")
    parser.add_argument(
        "-d", "--database", help="Path to the medival database")
    parser.add_argument(
        "-tr", "--tree", help="Directory containing The Time Tree of Life tree (TimeTree_v5_Final.nwk), and its preprocessed files")
    parser.add_argument(
        "-k", "--kraken", help="Path to Kraken2 database")
    parser.add_argument(
        "-s", "--species", type=str, help="(Optional) Species of all sequences in your input FASTA file (replace spaces with '_'). Cannot be used together with --speciesFile. If neither --species nor --speciesFile is provided, species will be auto-detected and saved to a file.")
    parser.add_argument(
        "--speciesFile", type=str, help="(Optional) Tab-separated file where the first column is sequence ID and the second column is its species. Cannot be used together with --species. If neither --species nor --speciesFile is provided, species will be auto-detected and saved to a file.")
    parser.add_argument(
        "-i", "--index", type=str, help="(Optional) Index of sequences in your db, of their species, length, and tree leaf names. Default is the one built in the medival database. Can specify a different index here if you want to manually define each sequence's species.")
    parser.add_argument(
        "-t", "--threads", type=int, help="(Optional) Number of threads to use. Highly recommend a large number of threads (default: 1)")
    parser.add_argument(
        "-c", "--chunk", type=int, help="(Optional) Size for splitting sequences to feed to blat (default: 100000)")
    parser.add_argument(
        "-minScore", type = int, help = "(Optional) Minimum alignment score threshold for blat (default: 30)")
    parser.add_argument(
        "-minIdentity", type = int, help = "(Optional) Minimum percent identity threshold (default: 90). Percent identity = ( match / Q_end - Q_start )*100")
    parser.add_argument(
        "--size-filter", dest="size_filter", type=int, help="(Optional) Filter out any regions smaller than this many bp (default: 150)")
    parser.add_argument(
        "--cluster-size", dest="cluster_size", type=int, help="(Optional) Combine any regions within this many bp of each other (default: 0)")

    return parser.parse_args()

def load_config_file(config_path):
    try:
        with open(config_path, 'r') as file:
            config_data = yaml.safe_load(file)
        return config_data
    except FileNotFoundError:
        print(f"Config file {config_path} not found.")
        return {}
    except yaml.YAMLError as e:
        print(f"Error parsing config file: {e}")
        return {}

def validate_config_keys(config_data):
    """Validate that all keys in config file are recognized parameters."""
    if not config_data:
        return
    
    # Define all valid configuration keys
    valid_keys = {
        'query', 'output', 'database', 'tree', 'index', 'kraken',
        'threads', 'chunk', 'species', 'minScore', 'minIdentity',
        'speciesFile', 'size_filter', 'cluster_size'
    }
    
    # Check for unrecognized keys
    config_keys = set(config_data.keys())
    unrecognized_keys = config_keys - valid_keys
    
    if unrecognized_keys:
        print("Error: Unrecognized parameters in config file:")
        for key in sorted(unrecognized_keys):
            print(f"  - '{key}'")
        print(f"\nValid parameters are:")
        for key in sorted(valid_keys):
            print(f"  - {key}")
        print(f"\nPlease check your config file for typos or unsupported parameters.")
        raise ValueError(f"Invalid config file parameters: {', '.join(sorted(unrecognized_keys))}")
    
    # Validate integer parameters
    int_keys = {'threads', 'chunk', 'minScore', 'minIdentity', 'size_filter', 'cluster_size'}
    for key in int_keys:
        if key in config_data:
            value = config_data[key]
            if not isinstance(value, int):
                try:
                    config_data[key] = int(value)
                except (ValueError, TypeError):
                    raise ValueError(f"Config parameter '{key}' must be an integer, got '{value}'")
            if config_data[key] < 0:
                raise ValueError(f"Config parameter '{key}' must be non-negative, got {config_data[key]}")

def species_seq_dict(speciesFile):
    """
    Load species mapping from a tab-separated file.
    Returns a dictionary mapping sequence IDs to species.
    """
    UNCLASSIFIED = "unclassified"
    NULL_VALUES = {"none", "null", "na", "n/a", "", "unknown"}

    species_map = {}
    try:
        with open(speciesFile, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    print(
                        f"Warning: Line {line_num} malformed "
                        f"(expected ≥2 columns). Marking as unclassified."
                    )
                    seq_id = parts[0]
                    species_map[seq_id] = UNCLASSIFIED
                    continue
                seq_id, species = parts[0], parts[1]

                species_clean = species.strip()

                # Catch undefined / null species values
                if species_clean.lower() in NULL_VALUES:
                    species_map[seq_id] = UNCLASSIFIED
                else:
                    species_map[seq_id] = species_clean
        return species_map
    
    except FileNotFoundError:
        print(f"Error: Species file {speciesFile} not found.")
        return {}
    except Exception as e:
        print(f"Error reading species file {speciesFile}: {e}")
        return {}

def merge_config_and_args(args, config_data=None):
    """Merge configuration file data with command line arguments.
    Command line arguments take priority over config file values."""
    
    if config_data is None:
        config_data = {}
    
    # Create a dictionary with all parameters
    # Start with config file values, then override with command line args
    merged_config = {}
    
    # Define parameter mappings and defaults
    param_mapping = {
        'query': ('query', None),
        'output': ('output', None),
        'database': ('database', None),
        'tree': ('tree', None),
        'index': ('index', 'default_placeholder'),
        'kraken': ('kraken', None),
        'threads': ('threads', 1),
        'chunk': ('chunk', 100000),
        'species': ('species', None),
        'speciesFile': ('speciesFile', None),
        'minScore': ('minScore', 30),
        'minIdentity': ('minIdentity', 90),
        'size_filter': ('size_filter', 150),
        'cluster_size': ('cluster_size', 0),
        'seq_species_dict': ('seq_species_dict', None),
    }
    
    # Process each parameter
    for config_key, (arg_key, default_value) in param_mapping.items():
        # Get value from config file first
        config_value = config_data.get(config_key, default_value)
        
        # Get command line argument value
        arg_value = getattr(args, arg_key, None)
        
        # Command line takes precedence if provided
        if arg_value is not None:
            merged_config[config_key] = arg_value
        else:
            merged_config[config_key] = config_value
    
    if merged_config['index'] == 'default_placeholder' or merged_config['index'] is None:
        if merged_config['database']:
            merged_config['index'] = f"{merged_config['database']}/medival_db_index.pkl"
        else:
            merged_config['index'] = None


    return merged_config

def validate_required_params(config_dict):
    required_params = ['query', 'output', 'database', 'tree', 'index', 'kraken']
    missing_params = []
    
    for param in required_params:
        if not config_dict.get(param):
            missing_params.append(param)
    
    if missing_params:
        print("Error: Missing required parameters:")
        for param in missing_params:
            print(f"  --{param}")
        print("\nEither provide them via command line or in your config file.")
        return False
        
    # Validate species configuration
    species = config_dict.get("species")
    speciesFile = config_dict.get("speciesFile")
    
    # Check that both species and speciesFile are not defined
    if species is not None and speciesFile is not None:
        print("Error: Both --species and --speciesFile are defined.")
        print("Please only define one or the other:")
        print("  --species assigns one species to every sequence in the input file")
        print("  --speciesFile is a tab-separated file that assigns a different species to each sequence")
        print("If you want us to automatically assign species, leave both undefined.")
        return False

    species_output_file = f"{config_dict['output']}/species_{Path(config_dict['query']).stem}.tsv"
    
    # Determine species handling mode
    if species is None and speciesFile is None:
        # Auto-detect mode: species will be written to and read from auto-generated file
        if Path(species_output_file).exists():
            # File already exists just use it
            config_dict['species_mode'] = 'file'
            config_dict['speciesFile'] = species_output_file
            config_dict['seq_species_dict'] = species_seq_dict(species_output_file)
        else:
            # File does not exist auto-detect and create it
            config_dict['species_mode'] = 'auto'
            config_dict['speciesFile'] = species_output_file
    elif species is not None:
        # User-defined single species for all sequences
        config_dict['species_mode'] = 'single'
    else:
        # User-provided species file
        config_dict['species_mode'] = 'file'
        config_dict['speciesFile'] = speciesFile
        config_dict['seq_species_dict'] = species_seq_dict(config_dict['speciesFile'])

    return True


def format_config_header(config_dict):
    """Return a # comment block summarising the run configuration."""
    skip = {'seq_species_dict'}
    lines = [
        '# MEDIVAL run configuration',
        f'# date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
    ]
    for key, val in config_dict.items():
        if key not in skip:
            lines.append(f'# {key}: {val}')
    lines.append('#')
    return '\n'.join(lines) + '\n'