__all__ = [
    "hello",
    "read_vcf"
]

def hello():
    return "hello"

def read_vcf(filepath: str):
    with open(filepath, "r") as fle:
        return fle.read()

class VCFFileReader:
    pass
