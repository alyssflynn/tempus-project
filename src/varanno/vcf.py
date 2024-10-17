import re
from dataclasses import dataclass
from .utils import VCF_META_KEYVAL, VCF_META_STRUCT


@dataclass(slots=True)
class Record:
    CHROM: str
    POS: int
    ID: str 
    REF: str
    ALT: str
    QUAL: str
    FILTER: str
    INFO: str
    FORMAT: str = None
    SAMPLE: str = None


class DuplicateHeaderError(Exception):
    pass


class ReaderError(Exception):
    def __init__(self, error: str, text: str = None, line_no: int = None):
        self.error = error
        self.text = text
        self.line_no = line_no
        super().__init__(error, text, line_no)
    

class Reader:
    _meta_multi = ('INFO', 'FILTER', 'FORMAT', 'ALT')
    _head_required = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
    # _head_allowed = {"FORMAT", "sample"} # TODO: allow multiple named samples e.g. "NA12878"

    def __init__(self, input_file_path: str = None, output_file_path: str = None):
        self.infile = input_file_path
        self.outfile = output_file_path
        self.metadata = {k: [] for k in self._meta_multi}
        self.header = None
        self.records = []

    def meta_structs(self):
        for name in self._meta_multi:
            for data in self.metadata.get(name):
                yield {"key": name, **data}

    def read(self, infile: str = None):
        self.infile = infile or self.infile

        if not self.infile:
            raise ReaderError("Input file missing")

        with open(self.infile, "r") as fle:
            for i, line in enumerate(fle):
                line_no = i + 1
                line = line.strip()

                if line.startswith("##"):
                    # TODO: yield metadata_factory(line) ?
                    self.parse_metadata(line, line_no)
                elif line.startswith("#"):
                    self.validate_head(line, line_no)
                elif self.header:
                    # TODO: yield record_factory(line) ?
                    self.build_record(line, line_no)
                else:
                    raise ReaderError("Line format invalid!", line, line_no)

    def validate_head(self, line: str, line_no: int = None):
        if self.header:
            raise ReaderError("Duplicate header found", line, line_no)
        
        head = self.splitrow(line[1:])
        if not self._head_required.issubset(set(head)):
            raise ReaderError("Missing required header fields")
        
        self.header = head
        
        # TODO: remove
        # valid = self._head_required.issubset(set(head)) \
        #     and self._head_required.union(self._head_allowed) == set(head)
        # print(head)
        # if valid:
        #     self.header = head
        # else:
        #     raise ReaderError(line, line_no)
        
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
        
        record = Record(*row)
        self.records.append(record)
        return record

    @staticmethod
    def splitrow(line: str):
        """Separates columns in a VCF file row. Lines should be tab-delimited."""
        return tuple(re.split(r"\t", line))