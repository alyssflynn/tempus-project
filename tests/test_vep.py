import json
import pytest
from varanno.vep import filter_vep_result, hgvs_string
from . import FIXTURES_DIR


@pytest.fixture
def vep_intron_variant():
    with open(FIXTURES_DIR.joinpath("vep_hgvs_response.json"), "r") as fle:
        return json.load(fle)
    

@pytest.fixture
def vep_intergenic_variant():
    with open(FIXTURES_DIR.joinpath("vep_hgvs_response_rec1.json"), "r") as fle:
        return json.load(fle)


@pytest.mark.parametrize("args,result", [
    (("1", "1158631", "A", "G"), "1:g.1158631A>G"),
    (("X", "25047897", "G", "GATACA"), "X:g.25047897G>GATACA"),
])
def test_hgvs_string_succeeds_for_valid_input(args, result):
    assert hgvs_string(*args) == result


def test_hgvs_string_fails_for_invalid_input():
    with pytest.raises(RuntimeError):
        hgvs_string("1", None, "A", "G")



def test_filter_vep_result_intron_variant(vep_intron_variant):
    assert filter_vep_result(vep_intron_variant, "C") == {
        "gene_id": "ENSG00000240498",
        "allele_string": "G/C",
        "variant_type": "SNV_SUB",
        "variant_effect": "intron_variant",
        'minor_allele_frequency': 0.4181,
    }


def test_filter_vep_result_intergenic_variant(vep_intergenic_variant):
    assert filter_vep_result(vep_intergenic_variant, "G") == {
        "gene_id": None,
        "allele_string": "A/G",
        "variant_type": "SNV_SUB",
        "variant_effect": "intergenic_variant",
        "minor_allele_frequency": None
    }
