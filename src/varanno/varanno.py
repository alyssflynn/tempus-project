from functools import cached_property
from .vcf import Reader, Record
from .utils import cast_float, parse_record_info, parse_format_sample, parse_genotype
from .vep import hgvs_string, filter_vep_result, vep_api_hgvs_get


__all__ = [
    "AnnotatedRecord"
]


class AnnotatedRecord:
    """
    Generates the following pieces of information for each variant:

    1. Depth of sequence coverage at the site of variation.
    2. Number of reads supporting the variant.
    3. Percentage of reads supporting the variant versus those supporting reference reads.
    4. Using the VEP hgvs API, get the gene of the variant, type of variation (substitution,
       insertion, CNV, etc.) and their effect (missense, silent, intergenic, etc.). The API
       documentation is available here: https://rest.ensembl.org/#VEP
    5. The minor allele frequency of the variant if available.
    6. Any additional annotations that you feel might be relevant.
    """
    def __init__(self, record: Record):
        self.record = record
        self.info = parse_record_info(self.record.INFO)
        self.sample = parse_format_sample(self.record.FORMAT, self.record.SAMPLE)
        self.genotype = parse_genotype(self.info.get("GT"))
        self.hgvs = hgvs_string(self.record.CHROM, self.record.POS, self.record.REF, self.record.ALT)
    
    @cached_property
    def _vep_api_response(self):
        return vep_api_hgvs_get(self.hgvs)
    
    @cached_property
    def vep_data(self):
        return filter_vep_result(self._vep_api_response, self.record.ALT)  

    @cached_property
    def total_coverage(self):
        """
        Depth of sequence coverage at the site of variation.

        Relevant metadata:
        ##INFO=<ID=TC..."Total coverage at this locus">
        ##INFO=<ID=TR..."Total number of reads containing this variant">
        ##FORMAT=<ID=NR..."Number of reads covering variant location in this sample">

        TODO: verify this is correct
        """
        return cast_float(self.info.get("TC"))

    @cached_property
    def n_reads_supporting_variant(self):
        """
        Number of reads supporting the variant.

        Relevant metadata:
        `##FORMAT=<ID=NV..."Number of reads containing variant in this sample">`
        `##INFO=<ID=TR..."Total number of reads containing this variant">`

        TODO: verify this is correct
        """
        return cast_float(self.sample.get("NV"))

    @cached_property
    def pct_reads_supporting_variant(self):
        """
        Percentage of reads supporting the variant versus those supporting reference reads.
        """
        n_var_reads = self.n_reads_supporting_variant
        total_reads_at_position = self.total_coverage

        try:
            n_supporting = float(n_var_reads)
            tot_reads = float(total_reads_at_position)
            return round((n_supporting / tot_reads) * 100, 4)
        
        except (ValueError, ZeroDivisionError):
            return None

    @cached_property
    def annotations(self):
        return { 
            **self.vep_data,
            "depth_of_sequence_coverage": self.total_coverage,
            "n_reads_supporting_variant": self.n_reads_supporting_variant,
            "pct_reads_supporting_variant": self.pct_reads_supporting_variant,
        }
