import pytest
from varanno.allele import variant_type


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
