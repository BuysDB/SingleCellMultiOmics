import math

def phred_to_prob(phred):
    try:
        return math.pow(10,-(ord(phred)-33)/10 )
    except ValueError:
        return 1

def hamming_distance(a,b):
    return sum((i!=j and i!='N' and j!='N' for i,j in zip(a,b)))
