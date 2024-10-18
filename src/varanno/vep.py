import requests
import pydash

HGVSString = str

# VEP_API_BASE_URL_LATEST = "https://rest.ensembl.org"
VEP_API_BASE_URL_GRCh37 = "https://grch37.rest.ensembl.org"


def hgvs_string(chromosome: str, pos: str | int, ref: str, alt: str) -> HGVSString:
    """Builds variant string in HGVS notation."""
    if not all((chromosome, pos, ref, alt)):
        raise RuntimeError
    return f"{chromosome}:g.{pos}{ref}>{alt}"


def first_element(data) -> dict:
    """Returns the first element in a list, otherwise returns the unchanged object."""
    if type(data) is list and data:
        return data[0]
    return data


def vep_api_hgvs_get(hgvs_string: HGVSString, species: str = "human") -> dict:
    """Fetch variant consequences from VEP API based on a HGVS notation.

    Endpoint: `GET vep/:species/hgvs/:hgvs_notation`

    Example: [/vep/human/hgvs/9:g.22125504G>C?](https://rest.ensembl.org/vep/human/hgvs/1:g.1158631A%3EG?content-type=application/json)

    ([VEP API docs](https://rest.ensembl.org/documentation/info/vep_hgvs_get))
    """
    url = f"{VEP_API_BASE_URL_GRCh37}/vep/{species}/hgvs/{hgvs_string}?"

    res = requests.get(url, headers={"Content-Type": "application/json"})
    if not res.ok:
        res.raise_for_status()

    return first_element(res.json())


def batch_vep_hgvs(
    hgvs_strings: list[HGVSString], species: str = "human"
) -> list[dict]:
    """Fetch variant consequences for multiple HGVS notations.

    Endpoint: `POST vep/:species/hgvs`
    ([VEP API docs](https://grch37.rest.ensembl.org/documentation/info/vep_hgvs_post))
    """
    url = f"{VEP_API_BASE_URL_GRCh37}/vep/{species}/hgvs"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = {"hgvs_notations": hgvs_strings}

    res = requests.post(url, headers=headers, json=data)
    if not res.ok:
        res.raise_for_status()

    return res.json()


def find_vep_gene_id(data: dict) -> str | None:
    """
    Finds the first instance of "gene_id" in the VEP API response.
    Gene info is located within the "transcript_consequences" array.
    """
    for tc in data.get("transcript_consequences", []):
        if gene_id := tc.get("gene_id"):
            return gene_id
    return None


def find_vep_allele_string(data: dict) -> str | None:
    """Retrieves "allele_string" from the VEP API data."""
    return data.get("allele_string")


def find_vep_variant_effect(data: dict) -> str | None:
    """Retrieves "most_severe_consequence" from the VEP API data."""
    return data.get("most_severe_consequence")


def find_vep_maf(data: dict, alt: str) -> float | None:
    """
    Finds the minor allele frequency within the VEP API response.
    Located at "colocated_variants[0].frequences.{alt}.af".
    """
    maf = None
    for covar in data.get("colocated_variants", []):
        if maf := pydash.get(covar, f"frequencies.{alt}.af"):
            return maf

    return maf
