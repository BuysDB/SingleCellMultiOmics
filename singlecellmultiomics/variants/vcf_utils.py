import pysam
from multiprocessing import Pool
from collections import defaultdict,OrderedDict
from itertools import product
from singlecellmultiomics.utils.sequtils import reverse_complement

def conversion_dict():
    conversions_single_nuc = ("CA", "CG", "CT", "TA", "TC", "TG")
    pattern_counts = OrderedDict()
    for ref, to in conversions_single_nuc:
        for context in product('ACGT',repeat=2 ):
            pattern_counts[(f'{context[0]}{ref}{context[1]}', to)] = 0
    return pattern_counts


def vcf_to_position_set(path: str):
    """
    Create a set of (contig, position) tuples from a vcf file

    Args:
        path(str): path to vcf file

    Returns:
        variant_locations(set): set of (contig (str), pos (int) ) tuples

    """
    with pysam.VariantFile(path) as v:
        contigs = v.header.contigs
    vs = set()
    with Pool() as workers:
        for r in workers.imap_unordered(_extract_all_variant_locations, (
                (contig, path)
                for contig in contigs)):
            vs.update(r)
    return vs


def _extract_all_variant_locations(args):
    """
    Extract all variant locations from a vcf file.
    Wrapped by vcf_to_position_set

    Args:
        args: (contig, variants_path) tuple

    Returns:
        locations (set)
    """
    contig, variants_path = args
    locations = set()
    with pysam.VariantFile(variants_path) as vcf:
        for record in vcf.fetch(contig):
            locations.add((record.chrom, record.pos))

    return locations


def _vcf_to_variant_contexts_detect(args):
    contig, ref_path, detected_variants_path, blacklist, whitelist = args
    pattern_obs = defaultdict(conversion_dict)  # Cell -> patterncounts

    with pysam.FastaFile(ref_path) as reference, pysam.VariantFile(detected_variants_path) as detected_vcf:
        for record in detected_vcf.fetch(contig):
            if len(record.ref) != 1:
                continue

            if blacklist is not None and (record.chrom, record.pos) in blacklist:
                continue
            if whitelist is not None and (record.chrom, record.pos) not in whitelist:
                continue

            origin_context = reference.fetch(record.contig, record.pos - 2, record.pos + 1).upper()
            for sample in record.samples:
                for allele in record.samples[sample].alleles:
                    if allele is None:
                        continue
                    if len(allele) != 1:
                        continue
                    if allele == record.ref:
                        continue

                    if not (record.ref + allele in set( ("CA", "CG", "CT", "TA", "TC", "TG"))):
                        context = reverse_complement(origin_context)
                        allele = reverse_complement(allele)
                    else:
                        context = origin_context
                    # print(allele,record.ref, context)
                    pattern_obs[sample][context, allele] += 1
    return pattern_obs


def vcf_to_variant_contexts(vcf_to_extract_contexts: str, reference_path: str, blacklist: set = None, whitelist: set = None):
    pattern_obs = defaultdict(conversion_dict)  # Cell -> { ('ACA','T') : obs, ... }

    with Pool() as workers, pysam.VariantFile(vcf_to_extract_contexts) as v:
        contigs = v.header.contigs
        for r in workers.imap_unordered(_vcf_to_variant_contexts_detect, (
                (contig, reference_path, vcf_to_extract_contexts, blacklist, whitelist) for contig in contigs)):
            for sample, sample_counts in r.items():
                for (context, allele), obs in sample_counts.items():
                    pattern_obs[sample][(context, allele)] += obs
    return pattern_obs

