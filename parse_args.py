import yaml
from dataclasses import dataclass
from pathlib import Path
import argparse

@dataclass
class Config:
    input_fasta: str #query fasta
    chunk_size: int # how we want to split up the query to run in blat
    output_dir: str # output dir name
    blat_db: Path # blat db path
    skani_db: Path # skani db path
    index: str #index that maps sequences in the database to their species and locations
    max_threads: int # number of threads the user wants to use
    kraken_db: Path # path to the kraken database
    tree: str # the folder that contains the tree and its preprocessing files
    remove: bool # is false when I am debugging and dont want files to be deleted
    species: str # user defined species of all of the sequences given in a query
    minScore: int # Is the minimum alignment score threshold, sets blat's parameter -minScore, default 30
    minIdentity: int # Is the threshold of the minimum percent identity a hit must have. Default: 95. Percent identity = ( match / Q_end - Q_start )*100
    overlap_filter: bool # Is false when you do not want to also do an overlap filter
    overlap_div_filter: bool # is false when you do not want to do an overlap and divergence filter

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
        "-i", "--index", type = str, help="(Optional) Index of sequences in your db, of their species, length, and tree leaf names. Default is the one built in the medival database. Can specify a different index here if you want to manually define each sequence's species.")
    parser.add_argument(
        "-t", "--threads", type=int, help="(Optional) Number of threads to use. Highly recommend a large number of threads (default: 1)")
    parser.add_argument(
        "-c", "--chunk", type=int, help="(Optional) Size for splitting sequences to feed to blat (default: 100000)")
    parser.add_argument(
        "-s", "--species", type = str, help = "(Optional) Species of the sequence(s) in your input fasta (replace spaces with '_'). If your file has multiple species, we recommend doing separate jobs for each species.")
    parser.add_argument(
        "-minScore", type = int, help = "(Optional) Minimum alignment score threshold for blat (default: 30)")
    parser.add_argument(
        "-minIdentity", type = int, help = "(Optional) Minimum percent identity threshold (default: 0). Percent identity = ( match / Q_end - Q_start )*100")

    # BooleanOptionalAction
    class BooleanOptionalAction(argparse.Action):
        def __init__(self, option_strings, dest, default=None, type=None, choices=None,required=False,help=None, metavar=None):
            _option_strings = []
            for option_string in option_strings:
                _option_strings.append(option_string)

                if option_string.startswith('--'):
                    option_string = '--no-' + option_string[2:]
                    _option_strings.append(option_string)
                elif option_string.startswith('-'):
                    raise ValueError('BooleanOptionalAction not supported for short options')

            super().__init__(option_strings=_option_strings,dest=dest,nargs=0,default=default,type=type,choices=choices,required=required,help=help,metavar=metavar)

        def __call__(self, parser, namespace, values, option_string=None):
            if option_string and option_string.startswith('--no-'):
                setattr(namespace, self.dest, False)
            else:
                setattr(namespace, self.dest, True)

        def format_usage(self):
            return ' | '.join(self.option_strings)

    # Boolean arguments using the custom BooleanOptionalAction
    parser.add_argument(
        "--remove",
        action=BooleanOptionalAction,
        default=None,
        help="(Optional), Clean up temporary files (default: False). Use --remove to enable, --no-remove to disable."
    )
    
    parser.add_argument(
        "--overlap-filter",
        dest="overlap_filter",
        action=BooleanOptionalAction,
        default=None,
        help="(Optional), Enable overlap filter (default: False). Use --overlap-filter to enable, --no-overlap-filter to disable."
    )
    
    parser.add_argument(
        "--overlap-div-filter", 
        dest="overlap_div_filter",
        action=BooleanOptionalAction,
        default=None,
        help="Enable overlap and divergence filter (default: False). Use --overlap-div-filter to enable, --no-overlap-div-filter to disable."
    )

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
        'threads', 'chunk', 'remove', 'species', 'minScore', 'minIdentity',
        'overlap_filter', 'overlap_div_filter'
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
    
    # Optionally, validate boolean values in config
    bool_keys = {'remove', 'overlap_filter', 'overlap_div_filter'}
    for key in bool_keys:
        if key in config_data:
            value = config_data[key]
            if not isinstance(value, (bool, int, str)):
                raise ValueError(f"Config parameter '{key}' must be a boolean, integer, or string, got {type(value).__name__}")
            if isinstance(value, str) and value.lower() not in ('true', 'false', 'yes', 'no', '1', '0', 'on', 'off'):
                raise ValueError(f"Config parameter '{key}' has invalid boolean value: '{value}'. Use true/false, yes/no, 1/0, or on/off")
    
    # Validate integer parameters
    int_keys = {'threads', 'chunk', 'minScore', 'minIdentity'}
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
        'remove': ('remove', False),
        'species': ('species', None),
        'minScore': ('minScore', 30),
        'minIdentity': ('minIdentity', 0),
        'overlap_filter': ('overlap_filter', False),
        'overlap_div_filter': ('overlap_div_filter', False)
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
    
    bool_params = ['remove', 'overlap_filter', 'overlap_div_filter']
    for param in bool_params:
        value = merged_config[param]
        if isinstance(value, str):
            # Handle string representations from config files
            merged_config[param] = value.lower() in ('true', 'yes', '1', 'on')
        elif isinstance(value, int):
            # Handle integer representations
            merged_config[param] = bool(value)

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
    return True