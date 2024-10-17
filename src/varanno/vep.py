import requests
from .allele import variant_type

HGVSString = str

VEP_API_BASE_URL = "https://rest.ensembl.org"


def hgvs_string(chromosome: str, pos: str | int, ref: str, alt: str) -> HGVSString:
    """Builds variant string in HGVS notation."""
    if not all((chromosome, pos, ref, alt)):
        raise RuntimeError
    return f"{chromosome}:g.{pos}{ref}>{alt}"


def first_element(data) -> dict:
    """Returns the first element in a list, otherwise returns the unchanged object."""
    if type(data) is list:
        return data[0]
    return data


def vep_api_hgvs_get(hgvs_string: HGVSString) -> list[dict]:
    """Fetch variant consequences from VEP API based on a HGVS notation.

    Endpoint: `GET vep/:species/hgvs/:hgvs_notation`

    Example: [/vep/human/hgvs/9:g.22125504G>C?](https://rest.ensembl.org/vep/human/hgvs/1:g.1158631A%3EG?content-type=application/json)

    ([VEP API docs](https://rest.ensembl.org/documentation/info/vep_hgvs_get))
    """
    url = f"{VEP_API_BASE_URL}/vep/human/hgvs/{hgvs_string}?"
    
    res = requests.get(url, headers={"Content-Type": "application/json"})
    if not res.ok:
        res.raise_for_status()  
    
    return res.json()


def find_vep_gene_id(data: dict) -> str | None:
    """
    Finds the first instance of "gene_id" in the VEP API response.
    Gene info is located within the "transcript_consequences" array. 
    """
    for tc in data.get("transcript_consequences", []):
        if (gene_id := tc.get("gene_id")):
            return gene_id
    

def find_vep_maf(data: dict, alt: str) -> float:
    """
    Finds the minor allele frequency within the VEP API response. 
    Located "colocated_variants[0].frequences.{alt}.af".
    """
    maf = None
    if (colocated_variant := first_element(data.get("colocated_variants", []))):
        if (frequencies := colocated_variant.get("frequencies", {})):
            maf = frequencies.get(alt, {}).get("af") 

    return maf


def filter_vep_result(api_output: dict, alt: str):
    """
    Using the VEP hgvs API output, gets the gene of the variant, type of variation (substitution,
    insertion, CNV, etc.) and their effect (missense, silent, intergenic, etc.).

    API docs: https://rest.ensembl.org/documentation/info/vep_region_get
    """
    data = first_element(api_output)
            
    # get gene of variant
    gene_id = find_vep_gene_id(data)

    # type of variation
    allele_string = data.get("allele_string")
    vtype = variant_type(allele_string) 

    # Variant effect
    effect = data.get("most_severe_consequence")

    # Minor allele frequency
    maf = find_vep_maf(data, alt)

    return {
        "gene_id": gene_id,
        "allele_string": allele_string,
        "variant_type": vtype,
        "variant_effect": effect,
        "minor_allele_frequency": maf
    }

