import json
import pytest
from pathlib import Path
from . import FIXTURES_DIR


BASE_DIR = Path(__file__).resolve().parent.parent
VCF_DATA_PATH = BASE_DIR.joinpath("data", "test_vcf_data.txt")



    

@pytest.fixture
def vcf_text():
    return VCF_DATA_PATH.read_text()


@pytest.fixture
def vep_response():
    with open(FIXTURES_DIR.joinpath("vep_response_5.33954511.json"), "r") as fle:
        return json.load(fle)
    

def test_blah():
    pass
