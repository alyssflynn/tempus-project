import json
import pytest
from varanno.vep import (
    hgvs_string, first_element, find_vep_gene_id, 
    find_vep_allele_string, find_vep_variant_effect, find_vep_maf,
    realign_hgvs_inputs_outputs
)
from . import FIXTURES_DIR


@pytest.fixture
def vep_hgvs_response():
    with open(FIXTURES_DIR.joinpath("hgvs_response.json"), "r") as fle:
        return json.load(fle)


@pytest.fixture
def hgvs_response_missing_result_output():
    with open(FIXTURES_DIR.joinpath("hgvs_response_missing_result.json"), "r") as fle:
        return json.load(fle)


@pytest.fixture
def vep_hgvs_response_data(vep_hgvs_response):
    return vep_hgvs_response[0]


@pytest.mark.parametrize("args,result", [
    (("1", "1158631", "A", "G"), "1:g.1158631A>G"),
    (("X", "25047897", "G", "GATACA"), "X:g.25047897G>GATACA"),
])
def test_hgvs_string_succeeds_for_valid_input(args, result):
    assert hgvs_string(*args) == result


def test_hgvs_string_fails_for_invalid_input():
    with pytest.raises(RuntimeError):
        hgvs_string("1", None, "A", "G")


@pytest.mark.parametrize("input,result", [
    ([{"a": 1, "b": 2}], {"a": 1, "b": 2}),
    ({"a": 1, "b": 2}, {"a": 1, "b": 2}),
])
def test_first_element(input, result):
    assert first_element(input) == result


def test_find_vep_gene_id(vep_hgvs_response_data):
    assert find_vep_gene_id(vep_hgvs_response_data) == "ENSG00000135744"


def test_find_vep_allele_string(vep_hgvs_response_data):
    assert find_vep_allele_string(vep_hgvs_response_data) == "T/C"


def test_find_vep_variant_effect(vep_hgvs_response_data):
    assert find_vep_variant_effect(vep_hgvs_response_data) == "missense_variant"


def test_find_vep_maf(vep_hgvs_response_data):
    assert find_vep_maf(vep_hgvs_response_data, "C") == 0.7051
    assert find_vep_maf(vep_hgvs_response_data, "A") is None


def test_realign_hgvs_inputs_outputs(hgvs_response_missing_result_output):
    hgvs_notations = [
        "1:g.1246004A>G", 
        "1:g.91859795_91859802TATGTGAdelinsCATGTGA,CATGTGG",
        "1:g.1647983_1647991TGGCTTACdelinsAGGCTTAT"
    ]
    result = list(realign_hgvs_inputs_outputs(hgvs_response_missing_result_output, hgvs_notations))
    assert result[0]["input"] == "1:g.1246004A>G"
    assert result[1]["error"] == "No data returned"
    assert result[2]["input"] == "1:g.1647983_1647991TGGCTTACdelinsAGGCTTAT"
