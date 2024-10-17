import pytest
from pathlib import Path
from varanno.vcf import Reader, ReaderError


BASE_DIR = Path(__file__).resolve().parent.parent
VCF_FILE = BASE_DIR.joinpath("data", "test_vcf_data.txt")
NON_VCF_FILE = BASE_DIR.joinpath("tests", "fixtures", "not_a_vcf_file.txt")
META_INIT = {k: [] for k in ('INFO', 'FILTER', 'FORMAT', 'ALT')}


@pytest.fixture
def reader():
    new_reader = Reader(VCF_FILE)
    new_reader.read()
    return new_reader


def test_vcf_read_succeeds():
    reader = Reader(VCF_FILE)
    reader.read()
    assert len(reader.records) == 11765
    assert reader.header == (
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
        "FILTER", "INFO", "FORMAT", "sample"
    )


def test_non_vcf_read_fails():
    with pytest.raises(ReaderError) as err:
        reader = Reader(NON_VCF_FILE)
        reader.read()
        assert str(err.value) == ""


def test_reader_parsed_metadata(reader):
    assert reader.metadata["fileformat"] == "VCFv4.0"
    assert reader.metadata["fileDate"] == "2016-06-21"
    assert reader.metadata["source"] == "Platypus_Version_0.8.1"
    assert len(reader.metadata["ALT"]) == 0
    assert len(reader.metadata["INFO"]) == 25


def test_reader_parsed_records(reader):
    assert len(reader.records) == 11765


@pytest.mark.parametrize("line, result", [
    (
        '##INFO=<ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">',
        {**META_INIT, 'INFO': [{'id': 'QD', 'number': '1', 'type': 'Float', 'description': 'Variant-quality/read-depth for this variant', 'source': None, 'version': None}]}
    ),
    (
        '##INFO=<ID=Size,Number=.,Type=Integer,Description="Size of reference call block">',
        {**META_INIT,'INFO': [{'id': 'Size', 'number': '.', 'type': 'Integer', 'description': 'Size of reference call block', 'source': None, 'version': None}]}
    ),
    (
        '##fileDate=2016-06-21',
        {**META_INIT, 'fileDate': '2016-06-21'}
    )
])
def test_reader_parse_metadata(line, result):
    reader = Reader()
    reader.parse_metadata(line)
    assert reader.metadata == result


def test_reader_validate_head():
    line = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample"
    reader = Reader()
    reader.validate_head(line) 
    assert reader.header == (
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
        "FILTER", "INFO", "FORMAT", "sample"
    )
