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
