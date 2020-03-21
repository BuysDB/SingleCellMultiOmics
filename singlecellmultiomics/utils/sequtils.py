import math
from pysam import FastaFile


def is_main_chromosome(chrom):
    """ Returns True when the chromsome is a main chromsome,
    not an alternative or other

    Args:
        chrom(str) : chromosome name

    Returns:
        is_main(bool) : True when the chromsome is a main chromsome

    """
    if chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith('JH') or chrom.startswith('GL') or chrom.startswith(
            'KI') or chrom.startswith('chrUn') or chrom.endswith('_random') or 'ERCC' in chrom or chrom.endswith('_alt') or "HLA-" in chrom:
        return False
    return True

def get_contig_list_from_fasta(fasta_path, with_length=False):
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
    for ref_base, query_base in zip(reference_seq, query_seq):
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
