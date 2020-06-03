from multiprocessing import Pool
from pysam import AlignmentFile, FastaFile
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from singlecellmultiomics.utils.sequtils import reverse_complement, get_context
from pysamiterators import CachedFasta
from array import array
from uuid import uuid4
from singlecellmultiomics.bamProcessing import merge_bams, get_contigs_with_reads


def get_covariate_key(read, qpos, refpos, reference, refbase, cycle_bin_size=3, k_rad=1):

    qual = read.query_qualities[qpos]
    qbase = read.query_sequence[qpos]

    if qbase == 'N':
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
                if refbase == 'N' or (read.reference_name, refpos) in known:
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
    return extract_covariates(**kwargs)


def extract_covariates_from_bam(bam_path, reference_path, known_variants, n_processes=None, bin_size=10_000_000,
            min_mapping_quality = 40,
            deduplicate = True,
            filter_qcfailed = True  ):

    global known
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

    with Pool(n_processes) as workers:

        for i, r in enumerate(
                        workers.imap_unordered(extract_covariates_wrapper, (
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


def recalibrate_base_calls(read, reference, joined_prob, covariate_kwargs):
    # @todo: make copy to save phred scores of soft-clipped bases
    # This array will contain all recalibrated phred scores:
    new_qualities = array('B', [0] * len(read.query_qualities))

    # Iterate all aligned pairs and replace phred score:

    for qpos, refpos in read.get_aligned_pairs(matches_only=True, with_seq=False):

        key = get_covariate_key(read, qpos, refpos, reference, None, **covariate_kwargs)
        try:
            phred = joined_prob[key]
        except KeyError:
            phred = joined_prob[None]
        new_qualities[qpos] = phred

    read.query_qualities = new_qualities


def _recalibrate_reads(bam_path, reference_path, contig, start, end, covariate_kwargs, **kwargs):
    # Recalibrate the reads in bam_path

    global joined_prob  # Global to share over multiprocessing
    # joined_prob contains  P(error| d), where d is a descriptor generated by get_covariate_key

    o_path = f'out_{uuid4()}.bam'

    # Open source bam file:
    with AlignmentFile(bam_path) as alignments, FastaFile(reference_path) as fa:
        # @todo: extract only selected region from fasta file:
        reference = CachedFasta(fa)
        # Open target bam file:
        with AlignmentFile(o_path, header=alignments.header, mode='wb') as out:
            # Iterate all reads in the source bam file:
            for read in alignments.fetch(contig, start, end):
                recalibrate_base_calls(read, reference, joined_prob, covariate_kwargs)
                out.write(read)

    pysam.index(o_path)
    return o_path


def __recalibrate_reads(kwargs):
    return _recalibrate_reads(**kwargs)


def recalibrate_reads(bam_path, target_bam_path, reference_path, n_processes, covariates,  covariate_kwargs, intermediate_bam_size=20_000_000):
    job_generation_args = {
        'contig_length_resource': bam_path,
        'bin_size': intermediate_bam_size,
        'fragment_size': 0
    }
    global joined_prob
    joined_prob = covariates


    with Pool(n_processes) as workers:
        intermediate_bams = list( workers.imap_unordered(__recalibrate_reads, (
                        {
                            'bam_path': bam_path,
                            'reference_path': reference_path,
                            'contig': contig,
                            'start': None,
                            'end': None,
                           # 'start_fetch': start_fetch,
                            #'end_fetch': end_fetch,
                            'covariate_kwargs': covariate_kwargs
                        }
                        for contig in list(get_contigs_with_reads(bam_path)))))
                        #for contig, start, end, start_fetch, end_fetch in
                        #blacklisted_binning_contigs(**job_generation_args))))

        merge_bams(intermediate_bams, target_bam_path)
