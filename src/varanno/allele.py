import re


ALLELE_STR = r"^([ATCG.-]+)(?:/([ATCG.+]+)){1,4}$"


VAR_TYPES = (
    ("SNV_SUB", r"^([ATCG])/(?!\1)([ATCG])$"),
    ("SNV_DEL", r"^([ATCG])/([.-])$"),
    ("SNV_INS", r"^([.-])/([ATCG])$"),
    ("MNV_SUB", r"^([ATCG]{2,})/(?!\1)([ATCG]{2,})$"),
    ("MNV_DEL", r"^([ATCG]+)/([.-]+)$"),
    ("MNV_INS", r"^([.-]+)/([ATCG]+)$"),
    ("CNV_DUP", r"^([ATCG])/(\1)$"),
    ("CNV_MUL", r"^([ATCG])(?:/(\1))+$"),
)


class InvalidAlleleStringError(Exception):
    pass


def variant_type(allele_string: str):
    """Determines the variant type from an allele string. Raises an error if unable.

    Possible outputs:
    - Single Nucleotide Variants (SNVs)
        - SNV_SUB: SNV substitution, e.g. "A/G"
        - SNV_DEL: SNV deletion, e.g. "A/-" or "A/."
        - SNV_INS: SNV insertion, e.g. "-/A" or "./A"
    - Multi-nucleotide Variants (MNVs)
        - MNV_SUB: MNV substitution, e.g. "AC/GT"
        - MNV_DEL: MNV deletion, e.g. "AC/.." (TODO: check validity)
        - MNV_INS: MNV insertion, e.g. "../AG" (TODO: check validity)
    - Copy Number Variants (CNVs)
        - CNV_DUP: duplication, e.g. "A/A"
        - CNV_MUL: multiplication, e.g. "A/A/A" or more
    - Complex Variants:
        - COMPLEX: more complex changes, e.g. "A/CTG"
    """
    for vartype, regex in VAR_TYPES:
        if re.match(regex, allele_string):
            return vartype

    if re.match(ALLELE_STR, allele_string):
        return "COMPLEX"

    raise InvalidAlleleStringError(f"Invalid allele string: {allele_string}")


def parse_genotype(gtstr: str):
    """Evaluates the genotype of the sample.

    - 0/0 : the sample is homozygous reference
    - 0/1 : the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
    - 1/1 : the sample is homozygous alternate

    Source: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format
    """
    return {
        "0/0": "homozygous_ref",
        "0/1": "heterozygous",
        "1/1": "homozygous_alt",
    }.get(gtstr)
