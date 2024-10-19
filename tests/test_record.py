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
def complex_record() -> Record:
    return Record(
        CHROM="1",
        POS="91859795",
        ID=".",
        REF="TATGTGA",
        ALT="CATGTGA,CATGTGG",
        QUAL="2962",
        FILTER="PASS",
        INFO="BRF=0.23;FR=0.5000,0.5000;HP=3;HapScore=2;MGOF=1;MMLQ=37;MQ=59.76;NF=115,51;"
             "NR=93,44;PP=2962,2892;QD=20;SC=TTGCCAGCAATATGTGATAAG;SbPval=0.54;Source=Platypus;"
             "TC=209;TCF=115;TCR=94;TR=208,95;WE=91859809;WS=91859785",
        FORMAT="GT:GL:GOF:GQ:NR:NV",
        SAMPLE="1/2:-1,-1,-1:1:99:209,209:208,95"
    )


@pytest.fixture
def variant_annotation():
    return VariantAnnotation("ENSG00000078808","A/G","SNV_SUB","synonymous_variant",0.5,160.0,156.0,97.5,"homozygous_alt")


def test_record_hgvs(record):
    assert record.hgvs == "5:g.33954511T>C"


def test_complex_record_hgvs(complex_record):
    assert complex_record.hgvs == "1:g.91859795_91859802TATGTGAdelinsCATGTGA,CATGTGG"


def test_annotation_factory(record, vep_hgvs_response_data):
    result = annotation_factory(record, vep_hgvs_response_data)
    assert result == VariantAnnotation(
        CHROM="5",
        POS="33954511",
        ID=".",
        REF="T",
        ALT="C",
        hgvs="5:g.33954511T>C",
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


def test_annotation_factory_complex(complex_record):
    result = annotation_factory(complex_record, {"error": "no data"})
    assert result == VariantAnnotation(
        CHROM=complex_record.CHROM,
        POS=complex_record.POS,
        ID=complex_record.ID,
        REF=complex_record.REF,
        ALT=complex_record.ALT,
        hgvs="1:g.91859795_91859802TATGTGAdelinsCATGTGA,CATGTGG",
        gene_id=None,
        allele_string=None,
        variant_type=None,
        variant_effect=None,
        minor_allele_frequency=None,
        depth_of_sequence_coverage=209.0,
        num_reads_supporting_variant="208,95",
        pct_reads_supporting_variant=None,
        genotype="unknown",
    )


@patch("varanno.vep.batch_vep_hgvs", vep_hgvs_response)
def test_annotate_batch(record):
    records = [record, record]
    results = list(annotate_batch(records))
    assert results == [
        VariantAnnotation(
            CHROM="5",
            POS="33954511",
            ID=".",
            REF="T",
            ALT="C",
            hgvs="5:g.33954511T>C",
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
            CHROM="5",
            POS="33954511",
            ID=".",
            REF="T",
            ALT="C",
            hgvs="5:g.33954511T>C",
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
