from dataclasses import dataclass
from .allele import variant_type, parse_genotype
from .utils import cast_float, parse_record_info, parse_format_sample
from .vep import (
    batch_vep_hgvs, hgvs_string, vep_api_hgvs_get,
    find_vep_gene_id, find_vep_maf, find_vep_allele_string, 
    find_vep_variant_effect
)


@dataclass(slots=True)
class Record:
    CHROM: str
    POS: int
    ID: str 
    REF: str
    ALT: str
    QUAL: str
    FILTER: str
    INFO: str
    FORMAT: str = None
    SAMPLE: str = None

    line_no: int = None
    hgvs: str = None

    def __post_init__(self):
        self.hgvs = hgvs_string(self.CHROM, self.POS, self.REF, self.ALT)
        

def pct_reads_supporting_variant(n_reads_supporting_variant: float, total_coverage: float):
    """
    Calculates percentage of reads supporting the variant versus those supporting reference reads.
    """
    try:
        var_reads = float(n_reads_supporting_variant)
        tot_reads = float(total_coverage)
        return round((var_reads / tot_reads) * 100, 4)
    
    except (ValueError, ZeroDivisionError):
        return None
    

@dataclass(slots=True)
class VariantAnnotation:
    CHROM: str
    POS: int
    ID: str 
    REF: str
    ALT: str
    gene_id: str
    allele_string: str
    variant_type: str
    variant_effect: str
    minor_allele_frequency: float
    depth_of_sequence_coverage: float
    num_reads_supporting_variant: float
    pct_reads_supporting_variant: float
    genotype: str
    
    
def annotation_factory(record: Record, vep_data: dict = None):
    """Generate annotations for a given variant record. 
    """
    info = parse_record_info(record.INFO)
    sample = parse_format_sample(record.FORMAT, record.SAMPLE)
    gt = parse_genotype(sample.get("GT"))

    # Depth of sequence coverage at the site of variation.
    total_coverage = cast_float(info.get("TC"))

    # Number of reads supporting the variant.
    num_var_reads = cast_float(sample.get("NV"))

    # Percentage of reads supporting the variant versus those supporting reference reads.
    pct_supporting = pct_reads_supporting_variant(num_var_reads, total_coverage)

    # VEP data
    if vep_data is None:
        vep_data = vep_api_hgvs_get(record.hgvs)

    # get gene of variant
    gene_id = find_vep_gene_id(vep_data)

    # type of variation
    allele_string = find_vep_allele_string(vep_data)
    vtype = variant_type(allele_string) 

    # Variant effect
    effect = find_vep_variant_effect(vep_data)

    # Minor allele frequency
    maf = find_vep_maf(vep_data, record.ALT)

    return VariantAnnotation(
        CHROM=record.CHROM,
        POS=record.POS,
        ID=record.ID,
        REF=record.REF,
        ALT=record.ALT,
        gene_id=gene_id,
        allele_string=allele_string,
        variant_type=vtype,
        variant_effect=effect,
        minor_allele_frequency=maf,
        depth_of_sequence_coverage=total_coverage,
        num_reads_supporting_variant=num_var_reads,
        pct_reads_supporting_variant=pct_supporting,
        genotype=gt
    )


def annotate_batch(records: list[Record]):
    hgvs_strings = [rec.hgvs for rec in records]
    for record, result in zip(records, batch_vep_hgvs(hgvs_strings)):
        yield annotation_factory(record, result)
