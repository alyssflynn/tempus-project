import argparse
from .varanno import VCFProcessor


def parse_args():
    """Parse args when provided via command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", 
        "--file", 
        required=True, 
        dest="infile",
        help="Input VCF file to process."
    )
    parser.add_argument(
        "-o", 
        "--out", 
        required=False, 
        dest="outdest",
        help="Output destination for annotation results", 
        default="varanno_output"
    )
    return parser.parse_args()


def run_annotation():
    args = parse_args()
    VCFProcessor(args.infile, args.outdest).process()
 