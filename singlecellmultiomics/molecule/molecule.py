from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
from singlecellmultiomics.fragment import Fragment
import collections
import singlecellmultiomics.tagtools as tagtools
import itertools
import numpy as np

def molecule_to_random_primer_dict(molecule, primer_length=6, primer_read=2): #1: read1 2: read2
    rp = collections.defaultdict(list)
    for fragment in molecule:
        if fragment[primer_read-1] is not None:
            hstart, hseq = tagtools.getRandomPrimerHash(fragment[primer_read-1], onStart=True, primerLength=6)
            rp[hstart, hseq].append(fragment)
    return rp

class Molecule():
    def __init__(self, fragments=None, cache_size=10_000):
        self.fragments  = []
        self.spanStart = None
        self.spanEnd = None
        self.chromosome = None
        self.cache_size = cache_size

        if fragments is not None:
            if type(fragments) is list:
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)


    def __len__(self):
        return len(self.fragments)

    def __repr__(self):
        frag_repr = '\n\t'.join([str(fragment) for fragment in self.fragments])
        return f"""Molecule
        with {len(self.fragments)} assinged fragments
        """ + frag_repr

    def get_umi(self):
        umi_abundance = collections.Counter()
        for fragment in self.fragments:
            umi_abundance[fragment.get_umi()]+=1
        return umi_abundance.most_common(1)[0][0]

    def get_sample(self):
        for fragment in self.fragments:
            return fragment.get_sample()

    def is_multimapped(self):
        for fragment in self.fragments:
            if not fragment.is_multimapped:
                return False
        return True

    def _add_fragment(self, fragment):
        self.fragments.append(fragment)
        # Update span:
        add_span = fragment.get_span()
        self.spanStart = add_span[1] if self.spanStart is None else min(add_span[1], self.spanStart)
        self.spanEnd = add_span[2] if self.spanEnd is None else max(add_span[2], self.spanEnd)
        self.chromosome = add_span[0]

    def add_fragment(self, fragment):
        if len(self.fragments)==0:
            self._add_fragment(fragment)
            return True

        for f in self.fragments:
            if f == fragment:
                # it matches this molecule:
                self._add_fragment(fragment)
                return True
        return False

    def can_be_yielded(self, chromsome, position):
        if chromsome!=self.chromosome:
            return False
        return position < (self.spanStart-self.cache_size*0.5) or position > (self.spanEnd+self.cache_size*0.5)

    def get_mean_rt_fragment_size(self):
        rt_reactions = molecule_to_random_primer_dict(self)
        amount_of_rt_reactions = len(rt_reactions)

        #this obtains the maximum fragment size:
        frag_chrom, frag_start, frag_end = pysamiterators.iterators.getListSpanningCoordinates([v for v in itertools.chain.from_iterable(self) if v is not None])

        #Obtain the fragment sizes of all RT reactions:
        rt_sizes = []
        for (rt_end,hexamer), fragment in rt_reactions.items():
            rt_chrom, rt_start, rt_end = pysamiterators.iterators.getListSpanningCoordinates(itertools.chain.from_iterable(fragment))
            rt_sizes.append([rt_end-rt_start])
        return np.mean(rt_sizes)


    def __getitem__(self, key):
        return self.fragments[key]

    def check_variants(self, variants, exclude_other_calls=True ): #when enabled other calls (non ref non alt will be set None)
        variant_dict = {}
        for variant in variants.fetch( self.chromosome, self.spanStart, self.spanEnd ):
            variant_dict[ (variant.chrom, variant.pos-1)] = (variant.ref, variant.alts)


        variant_calls = collections.defaultdict( collections.Counter )
        for fragment in self:

            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                allow_unsafe=(R1 is None))
            if start is None or end is None:
                continue

            for read in fragment:
                if read is None:
                    continue

                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(read):

                    if query_pos is None or ref_pos is None or ref_pos<start or ref_pos>end:
                        continue
                    query_base = read.seq[query_pos]

                    k =  (read.reference_name, ref_pos)
                    if k in variant_dict:
                        call = None
                        ref, alts = variant_dict[k]
                        if query_base == ref:
                            call = ('ref',query_base)
                        elif query_base in alts:
                            call = ('alt',query_base)

                        if not exclude_other_calls or call is not None:
                            variant_calls[k][call]+=1

        return variant_calls

    def iter_reads(self):
        for fragment in self.fragments:
            for read in fragment:
                if read is not None:
                    yield read

    def __iter__(self):
        for fragment in self.fragments:
            yield fragment

    """Get the total amount of methylated bases
    Parameters
    -------
    context : int
        3 or 4 base context

    Returns
    -------
    collections.Counter
    sum of methylated bases in contexts

    """
    def get_methylated_count(self, context=3):
        r = collections.Counter()


    def get_methylation_dict(self):
        methylated_positions =  collections.Counter () #chrom-pos->count
        methylated_state = dict()#chrom-pos->1, 0, -1
        for fragment in self:
            for read in fragment:
                if read is None or not read.has_tag('XM'):
                    continue
                methylation_status_string = read.get_tag('XM')
                i = 0
                for qpos, rpos, ref_base in read.get_aligned_pairs(with_seq=True):
                    if qpos is None:
                        continue
                    if ref_base is None:
                        continue
                    if rpos is None:
                        continue
                    methylation_status = methylation_status_string[i]
                    if methylation_status.isupper():
                        methylated_positions[(read.reference_name, rpos)]+=1
                        if methylated_state.get( (read.reference_name, rpos),1)==1 :
                            methylated_state[(read.reference_name, rpos)] = 1
                        else:
                            methylated_state[(read.reference_name, rpos)] = -1
                    else:
                        if methylation_status=='.':
                            pass
                        else:
                            if methylated_state.get( (read.reference_name, rpos),0)==0:
                                methylated_state[(read.reference_name, rpos)] = 0
                            else:
                                methylated_state[(read.reference_name, rpos)] = -1
                    i+=1
        return methylated_positions,methylated_state

def MoleculeIterator( alignments, moleculeClass=Molecule, fragmentClass=Fragment, check_eject_every=1000, **pysamArgs):
    molecules = []

    added_fragments = 0
    for R1,R2 in pysamiterators.iterators.MatePairIterator(alignments,performProperPairCheck=False,**pysamArgs):

        fragment = fragmentClass([R1,R2])
        if not fragment.is_valid() :
            continue

        added = False
        for molecule in molecules:
            if molecule.add_fragment(fragment):
                added = True
                break
        if not added:
            molecules.append(moleculeClass(fragment ))

        added_fragments+=1

        if added_fragments>check_eject_every:
            current_chrom, _, current_position = fragment.get_span()
            to_pop = []
            for i,m in enumerate(molecules):
                if m.can_be_yielded(current_chrom,current_position):
                    to_pop.append(i)

            for i,j in enumerate(to_pop):
                yield molecules.pop(i-j)

    # Yield remains
    return iter(molecules)
