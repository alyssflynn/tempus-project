import json
import pytest
from varanno import VCFProcessor
from unittest.mock import patch
from . import FIXTURES_DIR


CSV_HEAD = (
    "CHROM,POS,ID,REF,ALT,hgvs,gene_id,allele_string,variant_type,variant_effect,"
    "minor_allele_frequency,depth_of_sequence_coverage,num_reads_supporting_variant,"
    "pct_reads_supporting_variant,genotype"
)


@pytest.fixture
def tcf_path():
    return FIXTURES_DIR.joinpath("test_vcf_min.txt")


@pytest.fixture
def vep_hgvs_response():
    with open(FIXTURES_DIR.joinpath("hgvs_response.json"), "r") as fle:
        return json.load(fle)
    

def test_VCFProcessor_init(tmp_path):
    processor = VCFProcessor("in", tmp_path)
    assert processor.annotation_file == f"{tmp_path}/annotations.csv"
    assert processor.metadata_file == f"{tmp_path}/metadata.json"
    assert processor.error_file == f"{tmp_path}/errors.log"
    assert processor.log_file == f"{tmp_path}/tmp.log"


def test_VCFProcessor_validate_input_file_fails_if_not_exist(tmp_path):
    with pytest.raises(FileExistsError):
        processor = VCFProcessor("in", tmp_path)
        processor.validate_input_file()


@patch("os.path.exists")
def test_VCFPRocessor_validate_input_file_succeeds_if_exists(mock_exists, tmp_path):
    mock_exists.return_value = True
    processor = VCFProcessor("in", tmp_path)
    processor.validate_input_file()


@patch("varanno.record.batch_vep_hgvs")
def test_VCFPRocessor_process(mock_batch_vep_hgvs, tcf_path, tmp_path, vep_hgvs_response):
    mock_batch_vep_hgvs.return_value = vep_hgvs_response

    processor = VCFProcessor(tcf_path, tmp_path)
    processor.process()

    csvout = tmp_path.joinpath("annotations.csv")
    assert csvout.read_text().splitlines()[0] == CSV_HEAD
