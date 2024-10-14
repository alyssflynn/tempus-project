import pytest
from pathlib import Path
from varanno import hello, read_vcf

BASE_DIR = Path(__file__).resolve().parent.parent
VCF_DATA_PATH = BASE_DIR.joinpath("data", "test_vcf_data.txt")


@pytest.fixture
def vcf_data():
    return VCF_DATA_PATH.read_text()


def test_hello():
    assert hello() == "hello"


def test_read_vcf(vcf_data):
    assert read_vcf(VCF_DATA_PATH) == vcf_data
