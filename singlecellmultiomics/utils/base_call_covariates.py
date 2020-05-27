from multiprocessing import Pool

from pysam import AlignmentFile, FastaFile

from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from singlecellmultiomics.utils.sequtils import reverse_complement, get_context
from pysamiterators import CachedFasta
from array import array


def get_covariate_key(read, qpos, refpos, reference, refbase, cycle_bin_size=3, k_rad=1):

    qual = read.query_qualities[qpos]
    qbase = read.query_sequence[qpos]

    if qbase is 'N':
        return None

    context = get_context(read.reference_name, refpos, reference, qbase, k_rad)

    if 'N' in context:
        return None

    if read.is_reverse:
        cycle = len(read.query_sequence) - qpos  # -1
        context = reverse_complement(context)
    else:
        cycle = qpos + 1

    return qual, read.is_read2, int(round(cycle / cycle_bin_size)) * cycle_bin_size, context


def extract_covariates(bam_path: str,
                       reference_path: str,
                       contig: str,
                       start: int,
                       end: int,
                       start_fetch: int,
                       end_fetch: int,
                       filter_kwargs: dict,
                       covariate_kwargs: dict):
    """
    Count mismatches and matches for similar base-calls

    Returns:
        match_mismatch(dict) : dictionary ( covariate_key: [mismatches, matches], .. )
    """
    # known is a set() containing locations of known variation (snps)
    # @todo: extend to indels
    global known  # <- Locations, set of (contig, position) tuples to ignore

    joined = dict()

    # Filters which select which reads are used to estimate covariates:
    min_mapping_quality = filter_kwargs.get('min_mapping_quality', 0)
    deduplicate = filter_kwargs.get('deduplicate', False)
    filter_qcfailed = filter_kwargs.get('filter_qcfailed', False)

    with AlignmentFile(bam_path) as alignments,  FastaFile(reference_path) as fa:
        reference = CachedFasta(fa)  # @todo: prefetch selected region
        for read in alignments.fetch(contig, start_fetch, end_fetch):
            if (deduplicate and read.is_duplicate) or \
                    (read.is_qcfail and filter_qcfailed) or \
                    (read.mapping_quality < min_mapping_quality):
                continue

            for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):

                if refpos > end or refpos < start:  # Prevent the same location to be counted multiple times
                    continue

                refbase = refbase.upper()
                if refbase is 'N' or (read.reference_name, refpos) in known:
                    continue

                key = get_covariate_key(read, qpos, refpos, reference, refbase, **covariate_kwargs)
                if key is None:
                    continue

                matched = (refbase == read.query_sequence[qpos])
                try:
                    joined[key][matched] += 1
                except KeyError:
                    if matched:
                        joined[key] = array('l', [0, 1])
                    else:
                        joined[key] = array('l', [1, 0])
    return joined


def extract_covariates_wrapper(kwargs):
    extract_covariates(**kwargs)


def extract_covariates_from_bam(bam_path, reference_path, known_variants, n_processes=None, bin_size=10_000_000,
            min_mapping_quality = 40,
            deduplicate = True,
            filter_qcfailed = True  ):

    known = known_variants

    joined = dict()

    job_generation_args = {
        'contig_length_resource': bam_path,
        'bin_size': bin_size,
        'fragment_size': 0}

    filter_kwargs = {
        'min_mapping_quality': 40,
        'deduplicate': True,
        'filter_qcfailed': True
    }

    covariate_kwargs = {
        'cycle_bin_size': 3,
        'k_rad' : 1
    }

    jobs_total = sum(1 for _ in (blacklisted_binning_contigs(**job_generation_args)))

    with AlignmentFile(bam_path) as alignments, Pool(n_processes) as workers:

        for i, r in enumerate(
                        workers.imap_unordered(extract_covariates, (
                                   {
                                         'bam_path': bam_path,
                                         'reference_path': reference_path,
                                         'contig': contig,
                                         'start': start,
                                         'end': end,
                                         'start_fetch': start_fetch,
                                         'end_fetch': end_fetch,
                                         'filter_kwargs': filter_kwargs,
                                         'covariate_kwargs': covariate_kwargs
                                    }
                                   for contig, start, end, start_fetch, end_fetch in
                                             blacklisted_binning_contigs(**job_generation_args)))):
            print(round(100 * (i / jobs_total), 1), end='\r')

            for key, tf in r.items():
                try:
                    joined[key][0] += tf[0]
                    joined[key][1] += tf[1]
                except KeyError:
                    joined[key] = array('l', [0, 0])
                    joined[key][0] += tf[0]
                    joined[key][1] += tf[1]
    return joined
