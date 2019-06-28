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

        print(self.context_mapping)

    def digest(self, reads):

        #fragment_contig, fragment_start, fragment_end = pysamiterators.iterators.getListSpanningCoordinates([r for r in reads if r is not None])
        #if fragment_contig is None or fragment_end is None or fragment_start is None:
        #    return

        #fragment_size = fragment_end-fragment_start

        modified_contexts  = []
        unmodified_contexts  = []

        modified_quad_contexts  = []
        unmodified_quad_contexts  = []


        fragment_methylated_count = 0
        for read in reads:

            if read is None:
                continue
            if read.is_unmapped:
                continue
            if not read.has_tag('RS'):
                continue
            strand = read.get_tag('RS')
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
                qbase = read.seq[qpos].upper()

                methylated = False
                if ref_base=='C' and strand=='+':
                    context = self.reference.fetch(read.reference_name, rpos, rpos+3).upper()
                    quad_context = self.reference.fetch(read.reference_name, rpos-1, rpos+3).upper()
                    #print(ref_base, qbase, context)
                    if qbase=='T':
                        methylated=True
                    methylationStateString = self.context_mapping[methylated].get(context,'uU'[methylated])
                    if methylated:
                        modified_contexts.append(context)
                        modified_quad_contexts.append(quad_context)
                    else:
                        unmodified_contexts.append(context)
                        unmodified_quad_contexts.append(quad_context)

                elif ref_base=='G' and strand=='-':
                    origin = self.reference.fetch(read.reference_name, rpos-2, rpos+1).upper()
                    context = origin.translate(complement)[::-1]
                    quad_context= origin.translate(self.reference.fetch(read.reference_name, rpos-2, rpos+2))[::-1].upper()
                    if qbase=='A':
                        methylated=True
                    methylationStateString = self.context_mapping[methylated].get(context,'uU'[methylated])
                    if methylated:
                        modified_contexts.append(context)
                        modified_quad_contexts.append(quad_context)
                    else:
                        unmodified_contexts.append(context)
                        unmodified_quad_contexts.append(quad_context)

                if methylated:
                    fragment_methylated_count+=1
                methylationBlocks.append(methylationStateString)


            if len(methylationBlocks)>0:
                read.set_tag('XM',''.join(methylationBlocks))
        for read in reads:
            if read is None:
                continue
            if len(modified_contexts):
                read.set_tag('CM',','.join(list(sorted(set(modified_contexts)))))
                read.set_tag('QM',','.join(list(sorted(set(modified_quad_contexts)))))
            if len(unmodified_contexts):
                read.set_tag('CU',','.join(list(sorted(set(unmodified_contexts)))))
                read.set_tag('QU',','.join(list(sorted(set(unmodified_quad_contexts)))))
            read.set_tag('MC',fragment_methylated_count)
