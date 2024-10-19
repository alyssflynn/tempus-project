import logging
from dataclasses import dataclass, field
from typing import Any
from .allele import variant_type, parse_genotype
from .utils import cast_float
from .parse import parse_record_info, parse_format_sample
from .vep import (
    HGVSString,
    batch_vep_hgvs,
    hgvs_string,
    vep_api_hgvs_get,
    find_vep_gene_id,
    find_vep_maf,
    find_vep_allele_string,
    find_vep_variant_effect,
    realign_hgvs_inputs_outputs,
)


log = logging.getLogger(__name__)


@dataclass(slots=True)
class Record:
    CHROM: str
    POS: str | int
    ID: str
    REF: str
    ALT: str
    QUAL: str
    FILTER: str
    INFO: str
    FORMAT: str | None = None
    SAMPLE: str | None = None

    line_no: int | None = None
    hgvs: HGVSString = field(init=False)

    def __post_init__(self):
        self.hgvs = hgvs_string(self.CHROM, self.POS, self.REF, self.ALT)


def pct_reads_supporting_variant(
    n_reads_supporting_variant: float | Any, total_coverage: float | Any
):
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
    POS: str | int
    ID: str
    REF: str
    ALT: str
    hgvs: str
    gene_id: str | None
    allele_string: str | None
    variant_type: str | None
    variant_effect: str | None
    minor_allele_frequency: float | None
    depth_of_sequence_coverage: float | None
    num_reads_supporting_variant: float | None
    pct_reads_supporting_variant: float | None
    genotype: str | None


def annotation_factory(record: Record, vep_data: dict | None = None):
    """Generate annotations for a given variant record."""
    # log.info(f"Annotating record: {record}")

    info = parse_record_info(record.INFO)
    sample = (
        parse_format_sample(record.FORMAT, record.SAMPLE)
        if record.FORMAT and record.SAMPLE
        else {}
    )
    rec_gt = sample.get("GT")
    genotype = parse_genotype(rec_gt) if rec_gt else None
    if genotype is None:
        log.warning(
            f"Failed to parse genotype for record (pos: {record.POS}, GT:{rec_gt}), setting to 'unknown'"
        )
        genotype = "unknown"

    # Depth of sequence coverage at the site of variation.
    rec_tc = info.get("TC")
    total_coverage = cast_float(rec_tc) if rec_tc else None

    # Number of reads supporting the variant.
    rec_nv = sample.get("NV")
    num_var_reads = cast_float(rec_nv) if rec_nv else None

    # Percentage of reads supporting the variant versus those supporting reference reads.
    pct_supporting = pct_reads_supporting_variant(num_var_reads, total_coverage)

    # VEP data
    if vep_data is None:
        vep_data = vep_api_hgvs_get(record.hgvs)

    # get gene of variant
    gene_id = find_vep_gene_id(vep_data)

    # type of variation
    allele_string = find_vep_allele_string(vep_data)
    vtype = variant_type(allele_string) if allele_string else None

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
        hgvs=record.hgvs,
        gene_id=gene_id,
        allele_string=allele_string,
        variant_type=vtype,
        variant_effect=effect,
        minor_allele_frequency=maf,
        depth_of_sequence_coverage=total_coverage,
        num_reads_supporting_variant=num_var_reads,
        pct_reads_supporting_variant=pct_supporting,
        genotype=genotype,
    )


def annotate_batch(records: list[Record]):
    """Generates annotations for a batch of VCF records.

    Realigns outputs to appropriate inputs if the VEP API response
    contains fewer items than requested (it appears to filter out
    queries it can't process).
    """
    hgvs_strings = [rec.hgvs for rec in records]
    hgvs_results = batch_vep_hgvs(hgvs_strings)

    if len(hgvs_results) != len(hgvs_strings):
        log.warning(
            "VEP API result length does not match input length!",
            extra={
                "num_hgvs_notations": len(hgvs_strings),
                "num_results": len(hgvs_results),
            },
        )
        hgvs_results = list(realign_hgvs_inputs_outputs(hgvs_results, hgvs_strings))

    for record, result in zip(records, hgvs_results):
        yield annotation_factory(record, result)
