from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
from singlecellmultiomics.fragment import Fragment

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


    def check_variants(self, variants):
        variant_dict = {}
        for variant in variants.fetch( self.chromosome, self.spanStart, self.spanEnd ):
            variant_dict[ (variant.chrom, variant.pos)] = (variant.ref, variant.alts)


    def __iter__(self):
        for fragment in self.fragments:
            yield fragment

def MoleculeIterator( alignments, moleculeClass=Molecule, fragmentClass=Fragment, check_eject_every=1000):
    molecules = []

    added_fragments = 0
    for R1,R2 in pysamiterators.iterators.MatePairIterator(alignments,performProperPairCheck=False):

        fragment = fragmentClass([R1,R2])
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
