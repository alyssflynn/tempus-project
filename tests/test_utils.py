import pytest 
from varanno.utils import cast_float


@pytest.mark.parametrize("input, result", [
    ("4", 4.0),
    ("A", "A"),
    (["A"], ["A"])
])
def test_cast_float(input, result):
    assert cast_float(input) == result
