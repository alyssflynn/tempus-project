import pytest
from pathlib import Path
from varanno.vcf import Reader, ReaderError, Record


BASE_DIR = Path(__file__).resolve().parent.parent
VCF_FILE = BASE_DIR.joinpath("data", "test_vcf_data.txt")
NON_VCF_FILE = BASE_DIR.joinpath("tests", "fixtures", "not_a_vcf_file.txt")
META_INIT: dict[str, list] = {k: [] for k in ('INFO', 'FILTER', 'FORMAT', 'ALT')}


@pytest.fixture
def reader():
    new_reader = Reader(VCF_FILE)
    new_reader.load_records()
    return new_reader


@pytest.fixture
def header_line():
    return "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample"


@pytest.fixture
def record_line():
    return "1	1158631	.	A	G	2965	PASS	BRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621	GT:GL:GOF:GQ:NR:NV	1/1:-300.0,-43.88,0.0:3:99:160:156"


def test_vcf_read_succeeds():
    reader = Reader(VCF_FILE)
    reader.load_records()
    assert len(reader.records) == 11765
    assert reader.header == (
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
        "FILTER", "INFO", "FORMAT", "sample"
    )


def test_non_vcf_read_fails():
    with pytest.raises(ReaderError) as err:
        reader = Reader(NON_VCF_FILE)
        list(reader.read())

    assert err.value.error == "Line format invalid!"
    assert err.value.text == "Ceci n'est pas un VCF file"
    assert err.value.line_no == 1


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


def test_reader_validate_head_fails_on_duplicate():
    line = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample"
    reader = Reader()
    reader.validate_head(line) 
    with pytest.raises(ReaderError) as error:
        reader.validate_head(line) 
    assert error.value.error == "Duplicate header found"


def test_reader_build_record(header_line, record_line):
    reader = Reader()
    reader.validate_head(header_line) 
    rec = reader.build_record(record_line)
    assert isinstance(rec, Record)


def test_splitrow(header_line):
    assert Reader.splitrow(header_line) == (
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
        "FILTER", "INFO", "FORMAT", "sample"
    )
