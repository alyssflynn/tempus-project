import pytest
from varanno.allele import variant_type, parse_genotype


@pytest.mark.parametrize("allelestr,result", [
    ("A/G", "SNV_SUB"),
    ("A/-", "SNV_DEL"),
    ("A/.", "SNV_DEL"),
    ("./A", "SNV_INS"),
    ("AC/GT", "MNV_SUB"),
    ("AC/..", "MNV_DEL"),
    ("../CG", "MNV_INS"),
    ("A/A", "CNV_DUP"),
    ("A/A/A", "CNV_MUL"),
    ("A/AG", "COMPLEX"),
])
def test_variant_type(allelestr, result):
    assert variant_type(allelestr) == result


@pytest.mark.parametrize("gtstr, result", [
    ("0/0", "homozygous_ref"),
    ("0/1", "heterozygous"),
    ("1/1", "homozygous_alt")
])
def test_parse_genotype(gtstr, result):
    assert parse_genotype(gtstr) == result
