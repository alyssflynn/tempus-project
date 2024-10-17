import json
import pytest
from varanno.record import Record, VariantAnnotation, annotate_batch, annotation_factory
from unittest.mock import patch
from . import FIXTURES_DIR


@pytest.fixture
def vep_hgvs_response():
    with open(FIXTURES_DIR.joinpath("hgvs_response.json"), "r") as fle:
        return json.load(fle)
    

@pytest.fixture
def vep_hgvs_response_data(vep_hgvs_response):
    return vep_hgvs_response[0]


@pytest.fixture
def vep_intron_variant():
    with open(FIXTURES_DIR.joinpath("vep_hgvs_response.json"), "r") as fle:
        return json.load(fle)
    

@pytest.fixture
def vep_intergenic_variant():
    with open(FIXTURES_DIR.joinpath("vep_hgvs_response_rec1.json"), "r") as fle:
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
def variant_annotation():
    return VariantAnnotation("ENSG00000078808","A/G","SNV_SUB","synonymous_variant",0.5,160.0,156.0,97.5,"homozygous_alt")


def test_annotation_factory(record, vep_hgvs_response_data):
    result = annotation_factory(record, vep_hgvs_response_data)
    assert result == VariantAnnotation(
        gene_id='ENSG00000135744',
        allele_string='T/C',
        variant_type='SNV_SUB',
        variant_effect='missense_variant',
        minor_allele_frequency=0.7051,
        depth_of_sequence_coverage=106.0,
        num_reads_supporting_variant=105.0,
        pct_reads_supporting_variant=99.0566,
        genotype='homozygous_alt',
    )


# @patch("varanno.record.")
@patch("varanno.vep.batch_vep_hgvs", vep_hgvs_response)
def test_annotate_batch(record):
    records = [record, record]
    results = list(annotate_batch(records))
    assert results == [
        VariantAnnotation(
            gene_id='ENSG00000164175',
            allele_string='T/C',
            variant_type='SNV_SUB',
            variant_effect='missense_variant',
            minor_allele_frequency=0.8331,
            depth_of_sequence_coverage=106.0,
            num_reads_supporting_variant=105.0,
            pct_reads_supporting_variant=99.0566,
            genotype='homozygous_alt',
        ),
        VariantAnnotation(
            gene_id='ENSG00000164175',
            allele_string='T/C',
            variant_type='SNV_SUB',
            variant_effect='missense_variant',
            minor_allele_frequency=0.8331,
            depth_of_sequence_coverage=106.0,
            num_reads_supporting_variant=105.0,
            pct_reads_supporting_variant=99.0566,
            genotype='homozygous_alt',
        )
    ]
