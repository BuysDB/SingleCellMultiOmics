import collections
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
import pysamiterators.iterators
import itertools
complement = str.maketrans('ATGC', 'TACG')

class TAPSFlagger(DigestFlagger ):

    def __init__(self, reference, **kwargs):
        DigestFlagger.__init__(self, **kwargs )
        self.overlap_tag = 'XM'
        self.reference = reference


        """
        z unmethylated C in CpG context (CG)
        Z methylated C in CpG context (CG)
        x unmethylated C in CHG context ( C[ACT]G )
        X methylated C in CHG context   ( C[ACT]G )
        h unmethylated C in CHH context ( C[ACT][ACT] )
        H methylated C in CHH context ( C[ACT][ACT] )
        """
        self.context_mapping={}
        self.context_mapping[False] = {
            'CGA':'z',
            'CGT':'z',
            'CGC':'z',
            'CGG':'z',
            # CHG:
            'CAG':'x',
            'CCG':'x',
            'CTG':'x',
        }

        self.context_mapping[True] = { context:letter.upper() for context,letter in self.context_mapping[False].items() }

        #CHH:
        # h unmethylated C in CHH context ( C[ACT][ACT] )

        for x in list( itertools.product('ACT',repeat=2) ):
            self.context_mapping[True][ ''.join( ['C']+list(x)) ] =  'H'
            self.context_mapping[False][ ''.join( ['C']+list(x)) ] =  'h'



    def digest(self, reads):

        #fragment_contig, fragment_start, fragment_end = pysamiterators.iterators.getListSpanningCoordinates([r for r in reads if r is not None])
        #if fragment_contig is None or fragment_end is None or fragment_start is None:
        #    return

        #fragment_size = fragment_end-fragment_start

        for read in reads:
            if read is None:
                continue
            if read.is_unmapped:
                continue
            methylationBlocks = []
            for qpos, rpos, ref_base in read.get_aligned_pairs(with_seq=True):

                methylationStateString = '.'
                if qpos is None:
                    continue
                if ref_base is None:
                    continue
                if rpos is None:
                    continue

                ref_base = ref_base.upper()
                qbase = read.seq[qpos]
                strand = read.has_tag('RS')

                methylated = False
                if ref_base=='C' and strand=='+':
                    context = self.reference.fetch(read.reference_name, rpos, rpos+3)
                    if qbase=='T':
                        methylated=True
                    methylationStateString = self.context_mapping[methylated].get(context)

                elif ref_base=='G' and strand=='-':
                    context = self.reference.fetch(read.reference_name, rpos-3, rpos).translate(complement)[::-1]
                    if qbase=='A':
                        methylated=True
                    methylationStateString = self.context_mapping[methylated].get(context)

                methylationBlocks.append(methylationStateString)
            read.set_tag('XM',''.join(methylationBlocks))
