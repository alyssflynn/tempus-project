import pytest
from varanno import VCFProcessor
from . import FIXTURES_DIR



def test_VCFProcessor_init():
    processor = VCFProcessor("in", "out")
    assert processor.annotation_file == "out/annotations.csv"
    assert processor.metadata_file == "out/metadata.json"
    assert processor.error_file == "out/errors.log"
    assert processor.log_file == "out/tmp.log"
    # TODO add fileIO tests
