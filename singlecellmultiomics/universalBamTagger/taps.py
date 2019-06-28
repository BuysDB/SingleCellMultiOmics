import collections
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
import pysamiterators.iterators

class TAPSFlagger(DigestFlagger ):

    def __init__(self, reference, **kwargs):
        DigestFlagger.__init__(self, **kwargs )
        self.overlap_tag = 'XM'
        self.reference = reference

    def digest(self, reads):

        fragment_contig, fragment_start, fragment_end = pysamiterators.iterators.getListSpanningCoordinates([R1,R2])
        fragment_size = fragment_end-fragment_start

        for qpos, rpos, ref_base in read.get_aligned_pairs(with_seq=True):
            pass

    """
    z unmethylated C in CpG context (CG)
    Z methylated C in CpG context (CG)
    x unmethylated C in CHG context ( C[ACT]G )
    X methylated C in CHG context   ( C[ACT]G )
    h unmethylated C in CHH context ( C[ACT][ACT] )
    H methylated C in CHH context ( C[ACT][ACT] )
    """
