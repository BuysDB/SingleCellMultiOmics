import math

def phred_to_prob(phred):
    """Convert a phred score (ASCII) to a numeric probability

    returns:
        probability(float)
    """

    try:
        return math.pow(10,-(ord(phred)-33)/10 )
    except ValueError:
        return 1

def hamming_distance(a,b):
    return sum((i!=j and i!='N' and j!='N' for i,j in zip(a,b)))

complement_translate = str.maketrans('ATCGNatcgn','TAGCNtagcn')
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
