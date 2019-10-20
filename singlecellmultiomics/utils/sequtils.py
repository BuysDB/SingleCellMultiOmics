import math

def phred_to_prob(phred):
    """Convert a phred score (ASCII) or integer to a numeric probability
    Args:
        phred (str/int) : score to convert
    returns:
        probability(float)
    """

    try:
        if type(phred)==int:
            return math.pow(10,-(phred)/10 )
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

def split_nth(seq,separator,n):
    """
    Split sequence at the n-th occurence of separator

    Args:
        seq(str) : sequence to split
        separator(str): separator to split on
        n(int) : split at the n-th occurence
    """
    pos = 0
    for i in range(n):
        pos = seq.index(separator,pos+1)

    return seq[:pos],seq[pos+1:]
