import math
from pysam import FastaFile, AlignmentFile
from singlecellmultiomics.utils.prefetch import Prefetcher
from collections import Counter
import numpy as np
from pysamiterators import CachedFasta
from array import array
import os


class Reference(Prefetcher):
    """ This is a picklable wrapper to pass reference handles """

    def __init__(self):
        self.args = locals().copy()
        del self.args['self']

    def instance(self, arg_update):
        if 'self' in self.args:
            del self.args['self']
        clone = Reference(**self.args)
        return clone

        # Todo: exit statements
    def prefetch(self, contig, start, end):
        return FastaFile(**self.args)


def invert_strand_f(s):
    if s=='+':
        return '-'
    elif s=='-':
        return '+'
    return '.'

def get_contig_lengths_from_resource(resource ) -> dict:
    """
    Extract contig lengts from the supplied resouce (Fasta file or Bam/Cram/Sam )
    Returns:
        lengths(dict)
    """
    if type(resource) is AlignmentFile:
        return dict(zip(resource.references, resource.lengths))

    elif type(resource) is str:
        est_type = get_file_type(resource)

        if est_type in (AlignmentFile,FastaFile):
            with est_type(resource) as f:
                lens = dict(zip(f.references, f.lengths))

            return lens


    raise NotImplementedError('Unable to extract contig lengths from this resource')


def get_file_type(s: str):
    """Guess the file type of the input string, returns None when the file type can not be determined"""
    if s.endswith('.bam') or s.endswith('.cram') or s.endswith('.sam'):
        return AlignmentFile
    if s.endswith('.fa') or s.endswith('.fasta') or s.endswith('.fa.gz') or s.endswith('.fasta.gz'):
        return FastaFile

    return None

def create_fasta_dict_file(refpath: str, skip_if_exists=True):
    """Create index dict file for the reference fasta at refpath

    Args:
        refpath : path to fasta file

        skip_if_exists : do not generate the index if it exists

    Returns:
        dpath (str) : path to the dict index file
    """
    dpath = refpath.replace('.fa','').replace('.fasta','')+'.dict'
    if os.path.exists(dpath):
        return dpath

    with FastaFile(refpath) as reference, open(dpath,'w') as o:
        for ref, l in zip(reference.references, reference.lengths ):
            o.write(f'{ref}\t{l}\n')
    return dpath


def get_chromosome_number(chrom: str) -> int:
    """
    Get chromosome number (index) of the supplied chromosome:
        '1' -> 1, chr1 -> 1, returns -1 when not available, chrM -> -1
    """
    try:
        return int(chrom.replace('chr',''))
    except Exception as e:
        return -1

def is_autosome(chrom: str) -> bool:
    """ Returns True when the chromsome is an autosomal chromsome,
    not an alternative allele, mitochrondrial or sex chromosome

    Args:
        chrom(str) : chromosome name

    Returns:
        is_main(bool) : True when the chromsome is an autosome

    """
    return is_main_chromosome(chrom) and get_chromosome_number(chrom)!=-1



def is_main_chromosome(chrom: str) -> bool:
    """ Returns True when the chromsome is a main chromsome,
    not an alternative or other

    Args:
        chrom(str) : chromosome name

    Returns:
        is_main(bool) : True when the chromsome is a main chromsome

    """
    if chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith('JH') or chrom.startswith('GL') or chrom.startswith(
            'KI') or chrom.startswith('chrUn') or chrom.endswith('_random') or 'ERCC' in chrom or chrom.endswith('_alt') or "HLA-" in chrom or chrom.startswith('Un_') or 'decoy' in chrom:
        return False
    return True

def get_contig_list_from_fasta(fasta_path: str, with_length: bool=False) -> list:
    """Obtain list of contigs froma  fasta file,
        all alternative contigs are pooled into the string MISC_ALT_CONTIGS_SCMO

    Args:
        fasta_path (str or pysam.FastaFile) : Path or handle to fasta file

        with_length(bool): return list of lengths

    Returns:
        contig_list (list ) : List of contigs + ['MISC_ALT_CONTIGS_SCMO'] if any alt contig is present in the fasta file
        """

    contig_list = []
    has_alt = False
    if with_length:
        lens = []

    if type(fasta_path) is str:
        fa = FastaFile(fasta_path)
    elif type(fasta_path) is FastaFile:
        fa = fasta_path
    else:
        raise TypeError('Supply pysam.FastaFile or str')

    for reference, length in zip(fa.references, fa.lengths):
        if is_main_chromosome(reference):
            contig_list.append(reference)
            if with_length:
                lens.append(length)
        else:
            has_alt = True

    # Close handle if we just opened one
    if type(fasta_path) is str:
        fa.close()

    if has_alt:
        contig_list.append('MISC_ALT_CONTIGS_SCMO')
        if with_length:
            lens.append(None)

    if with_length:
        return contig_list, lens

    return contig_list

def phred_to_prob(phred):
    """Convert a phred score (ASCII) or integer to a numeric probability
    Args:
        phred (str/int) : score to convert
    returns:
        probability(float)
    """

    try:
        if isinstance(phred, int):
            return math.pow(10, -(phred) / 10)
        return math.pow(10, -(ord(phred) - 33) / 10)
    except ValueError:
        return 1


def hamming_distance(a, b):
    return sum((i != j and i != 'N' and j != 'N' for i, j in zip(a, b)))


complement_translate = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')


def reverse_complement(seq):
    """Obtain reverse complement of seq

    returns:
        reverse complement (str)
    """
    return seq.translate(complement_translate)[::-1]


def complement(seq):
    """Obtain complement of seq

    returns:
        complement (str)
    """
    return seq.translate(complement_translate)


def split_nth(seq, separator, n):
    """
    Split sequence at the n-th occurence of separator

    Args:
        seq(str) : sequence to split
        separator(str): separator to split on
        n(int) : split at the n-th occurence
    """
    pos = 0
    for i in range(n):
        pos = seq.index(separator, pos + 1)

    return seq[:pos], seq[pos + 1:]


def create_MD_tag(reference_seq, query_seq):
    """Create MD tag
    Args:
        reference_seq (str) : reference sequence of alignment
        query_seq (str) : query bases of alignment
    Returns:
        md_tag(str) : md description of the alignment
    """
    no_change = 0
    md = []
    for ref_base, query_base in zip(reference_seq.upper(), query_seq):
        if ref_base.upper() == query_base:
            no_change += 1
        else:
            if no_change > 0:
                md.append(str(no_change))
            md.append(ref_base)
            no_change = 0
    if no_change > 0:
        md.append(str(no_change))
    return ''.join(md)


def prob_to_phred(prob: float):
    """
    Convert probability of base call being correct into phred score
    Values are clipped to stay within 0 to 60 phred range

    Args:
        prob  (float): probability of base call being correct

    Returns:
        phred_score (byte)
    """
    return np.rint(-10 * np.log10(np.clip(1-prob, 1-0.999999, 0.999999))).astype('B')


def get_context(contig: str, position: int, reference: FastaFile, ibase: str = None, k_rad: int = 1):
    """
    Args:
        contig: contig of the location to extract context
        position: zero based position
        reference: pysam.FastaFile handle or similar object which supports .fetch()
        ibase: single base to inject into the middle of the context
        k_rad: radius to extract

    Returns:
        context(str) : extracted context with length k_rad*2 + 1

    """
    if ibase is not None:
        ctx = reference.fetch(contig, position-k_rad, position+k_rad+1).upper()
        return ctx[:k_rad]+ibase+ctx[1+k_rad:]
    else:
        return reference.fetch(contig, position-k_rad, position+k_rad+1).upper()

def base_probabilities_to_likelihood(probs: dict):
    probs['N'] = [1-p  for base, ps in probs.items() for p in ps if base != 'N' ]
    return {base:np.product(v)/np.power(0.25, len(v)-1) for base,v in probs.items() }

def likelihood_to_prob(likelihoods):
    total_likelihood = sum(likelihoods.values())
    return {key: value / total_likelihood
    for key, value in likelihoods.items()}


def phredscores_to_base_call(probs: dict):
    """
    Perform base calling on a observation dictionary.
    Returns N when there are multiple options with the same likelihood

    Args:
        probs: dictionary with confidence scores probs = {
            'A':[0.95,0.99,0.9],
            'T':[0.1],
        }

    Returns:
        base(str) : Called base
        phred(float) : probability of the call to be correct
    """
    # Add N:
    likelihood_per_base = base_probabilities_to_likelihood(probs)
    total_likelihood = sum(likelihood_per_base.values())
    base_probs = Counter({base:p/total_likelihood for base, p in likelihood_per_base.items() }).most_common()

    # We cannot make a base call when there are no observations or when the most likely bases have the same prob
    if len(base_probs) == 0 or (len(base_probs) >= 2 and base_probs[0][1] == base_probs[1][1]):
        return 'N', 0

    return (base_probs[0][0],  base_probs[0][1])


def pick_best_base_call( *calls ) -> tuple:
    """ Pick the best base-call from a list of base calls

    Example:
        >>> pick_best_base_call( ('A',32), ('C',22) ) )
        ('A', 32)

        >>> pick_best_base_call( ('A',32), ('C',32) ) )
        None

    Args:
        calls (generator) : generator/list containing tuples

    Returns:
        tuple (best_base, best_q) or ('N',0) when there is a tie
    """
    # (q_base, quality, ...)
    best_base, best_q = None, -1
    tie = False

    for call in calls:
        if call is None:
            continue
        if call[1]>best_q:
            best_base= call[0]
            best_q=call[1]
            tie=False
        elif call[1]==best_q and call[0]!=best_base:
            tie=True

    if tie or best_base is None:
        return ('N',0)

    return best_base, best_q


def read_to_consensus_dict(read, start: int = None, end: int = None, only_include_refbase: str = None, skip_first_n_cycles:int = None, skip_last_n_cycles: int = None, min_phred_score: int  = None):
    """
    Obtain consensus calls for read, between start and end
    """

    if read is None:
        return dict()

    return { (read.reference_name, refpos):
                    (read.query_sequence[qpos],
                     read.query_qualities[qpos],
                     refbase
                    )

         for qpos, refpos, refbase in read.get_aligned_pairs(
                                             matches_only=True,
                                             with_seq=True)

         if (start is None or refpos>=start) and \
            (end is None or refpos<=end) and \
            (min_phred_score is None or read.query_qualities[qpos]>=min_phred_score) and \
            (skip_last_n_cycles is None or ( read.is_reverse and qpos>skip_last_n_cycles) or (not read.is_reverse and qpos<read.infer_query_length()-skip_last_n_cycles)) and \
            (skip_first_n_cycles is None or ( not read.is_reverse and qpos>skip_first_n_cycles) or ( read.is_reverse and qpos<read.infer_query_length()-skip_first_n_cycles)) and \
            (only_include_refbase is None or refbase.upper()==only_include_refbase)
           }


def get_consensus_dictionaries(R1, R2, only_include_refbase=None, dove_safe=False, min_phred_score=None, skip_first_n_cycles_R1=None, skip_last_n_cycles_R1=None,skip_first_n_cycles_R2=None, skip_last_n_cycles_R2=None, dove_R2_distance=0, dove_R1_distance=0  ):

    assert (R1 is None or R1.is_read1) and (R2 is None or R2.is_read2)

    if dove_safe:
        if R1 is None or R2 is None:
            raise ValueError(
                'Its not possible to determine a safe region when the alignment of R1 or R2 is not specified')


        if R1.is_reverse and not R2.is_reverse:
            start, end = R2.reference_start + dove_R2_distance, R1.reference_end - dove_R1_distance -1
        elif not R1.is_reverse and R2.is_reverse:
            start, end = R1.reference_start + dove_R1_distance, R2.reference_end - dove_R2_distance -1
        else:
            raise ValueError('This method only works for inwards facing reads')
    else:
        start, end = None, None

    return read_to_consensus_dict(R1, start, end, only_include_refbase=only_include_refbase, skip_last_n_cycles=skip_last_n_cycles_R1, skip_first_n_cycles=skip_first_n_cycles_R1,min_phred_score=min_phred_score), \
           read_to_consensus_dict(R2, start, end, only_include_refbase=only_include_refbase, skip_last_n_cycles=skip_last_n_cycles_R2, skip_first_n_cycles=skip_last_n_cycles_R2, min_phred_score=min_phred_score)
