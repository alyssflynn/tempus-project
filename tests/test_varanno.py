import json
import pytest
from pathlib import Path
from varanno import AnnotatedRecord
from varanno.vcf import Record
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


@pytest.fixture
def record() -> Record:
    return Record(
        CHROM="5",
        POS="33954511",
        ID=".",
        REF="T",
        ALT="C",
        QUAL="2965",
        FILTER="PASS",
        INFO="BRF=0.18;FR=1.0000;HP=1;HapScore=1;MGOF=1;MMLQ=35;MQ=59.86;NF=69;NR=36;PP=2965;QD=20;"
             "SC=ACAGGAAGGCTGTCCATCCAA;SbPval=0.55;Source=Platypus;TC=106;TCF=69;TCR=37;TR=105;"
             "WE=33954519;WS=33954501",
        FORMAT="GT:GL:GOF:GQ:NR:NV",
        SAMPLE="1/1:-300.0,-29.78,0.0:1:99:106:105"
    )

@pytest.fixture
def annotated_record(record, vep_response) -> AnnotatedRecord:
    annot = AnnotatedRecord(record)
    annot._vep_api_response = vep_response
    return annot


def test_record_annotations(annotated_record):
    assert annotated_record.annotations == {
        "gene_id": "ENSG00000164175",
        "allele_string": "T/C",
        "variant_type": "SNV_SUB",
        "variant_effect": "splice_region_variant",
        "depth_of_sequence_coverage": 106.0,
        "minor_allele_frequency": 0.4181,
        "n_reads_supporting_variant": 105.0,
        "pct_reads_supporting_variant": 99.0566,
    }
