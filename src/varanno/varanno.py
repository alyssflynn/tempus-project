import csv
import os
import sys
import logging
from dataclasses import asdict
from .fileio import write_metadata_json, write_logs
from .record import VariantAnnotation
from .vcf import Reader, Record


__all__ = [
    "VCFProcessor",
    "Reader",
    "Record"
]

stdout_handler = logging.StreamHandler(stream=sys.stdout)
handlers = [stdout_handler]

logging.basicConfig(
    level=logging.DEBUG, 
    format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s',
    handlers=handlers
)

log = logging.getLogger(__name__)


class VCFProcessor:
    def __init__(self, infile: str, outdir: str, allow_overrides: bool = True):
        self.infile = infile
        self.outdir = outdir
        self.allow_overrides = allow_overrides
        
        # Ensure output dir exists, override if allowed
        os.makedirs(self.outdir, exist_ok=self.allow_overrides)

        # Set output file paths
        self.annotation_file = os.path.join(self.outdir, 'annotations.csv')
        self.metadata_file = os.path.join(self.outdir, 'metadata.json')
        self.error_file = os.path.join(self.outdir, 'errors.log')
        self.log_file = os.path.join(self.outdir, 'tmp.log')
        log.addHandler(logging.FileHandler(filename=self.log_file))

    def process(self):
        self.validate_input_file()

        log.info(f"Annotating VCF! {self.infile} -> {self.outdir}")
        self.reader = Reader(self.infile)

        # Generate variant annotations and write to file
        annotation_gen = self.reader.annotation_generator()
        self.write_record_annotations(annotation_gen, self.annotation_file)

        # Write metadata to file
        log.info(f"Writing metadata JSON -> {self.metadata_file}")
        write_metadata_json(self.reader.metadata, self.metadata_file)

        # Write errors to file
        if self.reader.errors:
            log.warning(f"Logging {len(self.reader.errors)} errors -> {self.error_file}")
            write_logs([err.logstr() for err in self.reader.errors], self.error_file)
        
    def validate_input_file(self):
        if not os.path.exists(self.infile):
            raise FileExistsError("Input file not found")
        
    def write_record_annotations(self, annotation_gen):
        log.info(f"Generating record annotations -> {self.annotation_file}")

        with open(self.annotation_file, "wt") as fle:
            writer = csv.DictWriter(fle, fieldnames=VariantAnnotation.__slots__)
            writer.writeheader()

            for variant in map(asdict, annotation_gen):
                writer.writerow(variant)

