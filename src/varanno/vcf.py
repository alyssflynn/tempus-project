import re
import logging
from .record import Record, annotate_batch
from .parse import VCF_META_KEYVAL, VCF_META_STRUCT


log = logging.getLogger(__name__)


class ReaderError(Exception):
    def __init__(self, error: str, text: str = None, line_no: int = None):
        self.error = error
        self.text = text
        self.line_no = line_no
        super().__init__(error, text, line_no)
        log.error(self)
    
    def logstr(self):
        return f"{self.error}: {self.text} [{self.line_no}]"
    

class Reader:
    _meta_multi = ('INFO', 'FILTER', 'FORMAT', 'ALT')
    _head_required = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}

    def __init__(self, infile: str = None):
        self.infile = infile
        self._init_meta()

    def _init_meta(self):
        self.metadata = {k: [] for k in self._meta_multi}
        self.header = None
        self.records = []
        self.errors = []

    def meta_structs(self):
        for name in self._meta_multi:
            for data in self.metadata.get(name):
                yield {"key": name, **data}

    def read(self, infile: str = None):
        self.infile = infile or self.infile
        self._init_meta()

        if not self.infile:
            raise ReaderError("Input file missing")

        log.info(f"Reading VCF file: {self.infile}")
        with open(self.infile, "r") as fle:
            for i, line in enumerate(fle):
                line_no = i + 1
                line = line.strip()

                try:
                    if line.startswith("##"):
                        self.parse_metadata(line, line_no)
                    elif line.startswith("#"):
                        self.validate_head(line, line_no)
                    elif self.header:
                        record = self.build_record(line, line_no)
                        yield record
                    else:
                        raise ReaderError("Line format invalid!", line, line_no)
                    
                except ReaderError as err:
                    log.warning(err)
                    self.errors.append(err)
                    raise err
    
    def load_records(self):
        self.records = list(self.read())
        return self.records

    def validate_head(self, line: str, line_no: int = None):
        if self.header:
            raise ReaderError("Duplicate header found", line, line_no)
        
        head = self.splitrow(line[1:])
        if not self._head_required.issubset(set(head)):
            raise ReaderError("Missing required header fields")
        
        self.header = head
        
    def parse_metadata(self, line: str, line_no: int = None):
        if (m := re.match(VCF_META_STRUCT, line)):
            data = m.groupdict() 
            key = data.pop("key")
            self.metadata[key].append(data)

        elif (m := re.match(VCF_META_KEYVAL, line)):
            key = m.group('key')
            value = m.group('value') 
            self.metadata[key] = value
        else:
            raise ReaderError("Invalid metadata!", line, line_no)

    def build_record(self, line: str, line_no: int = None):
        row = self.splitrow(line)

        if len(row) != len(self.header):
            raise ReaderError("Invalid record format!", line, line_no)
        
        return Record(*row, line_no=line_no)
        
    def annotation_generator(self, batch_size: int = 50):
        """Yield batches of records annotations from the .read() record generator."""
        log.info(f"Annotating records: Batch size {batch_size}")

        batch = []
        batch_no = 1
        for item in self.read():
            batch.append(item)
            if len(batch) == batch_size:
                yield from annotate_batch(batch)
                log.info(f"Successfully processed batch #{batch_no}")
                batch = []
                batch_no += 1
                
        # Yield any remaining items in the final batch
        if batch:
            yield from annotate_batch(batch)

    @staticmethod
    def splitrow(line: str):
        """Separates columns in a VCF file row. Lines should be tab-delimited."""
        return tuple(re.split(r"\t", line))
    