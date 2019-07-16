from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
from singlecellmultiomics.fragment import Fragment
import collections
import itertools
import numpy as np

def molecule_to_random_primer_dict(molecule, primer_length=6, primer_read=2, max_N_distance=0): #1: read1 2: read2
    rp = collections.defaultdict(list)

    # First add all reactions without a N in the sequence:
    for fragment in molecule:
        if fragment[primer_read-1] is not None:
            hstart, hseq = fragment.get_random_primer_hash()
            if not 'N' in hseq:
                rp[hstart, hseq].append(fragment)

    # Try to match reactions with N with reactions without a N
    for fragment in molecule:
        if fragment[primer_read-1] is not None:
            hstart, hseq = fragment.get_random_primer_hash()
            if 'N' in hseq:
                # find nearest
                for other_start, other_seq in rp:
                    if other_start!=hstart:
                        continue

                    if 'N' in other_seq:
                        continue

                    if hamming_distance(hseq,other_seq)<=max_N_distance:
                        rp[other_start, other_seq].append(fragment)
    return rp

class Molecule():
    """Molecule class, contains one or more associated fragments

    Attributes:
        fragments (list): associated fragments

        spanStart (int): starting coordinate of molecule, None if not available

        spanEnd (int): ending coordinate of molecule, None if not available

        chromosome (str): mapping chromosome of molecule, None if not available

        cache_size (int): radius of molecule assignment cache

        strand (bool): mapping strand.
            True when strand is REVERSE
            False when strand is FORWARD
            None when strand is not determined
    """

    def __init__(self, fragments=None, cache_size=10_000):
        """Initialise Molecule

        Parameters
        ----------
        fragments :  list(singlecellmultiomics.fragment.Fragment)
            Fragment to assign to Molecule. More fragments can be added later

        cache_size (int): radius of molecule assignment cache
        """

        self.fragments  = []
        self.spanStart = None
        self.spanEnd = None
        self.chromosome = None
        self.cache_size = cache_size
        self.strand = None

        if fragments is not None:
            if type(fragments) is list:
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)


    def __len__(self):
        return len(self.fragments)



    def get_strand(self):
        """Obtain mapping strand of molecule
        Returns:
            strand : True,False,None
                True when strand is REVERSE
                False when strand is FORWARD
                None when strand is not determined
        """
        return self.strand

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
        """Obtain sample

        Returns:
            sample (str):
                Sample associated with the molecule. Usually extracted from SM tag.
                Calls fragment.get_sample() to obtain the sample
        """
        for fragment in self.fragments:
            return fragment.get_sample()


    def get_cut_site(self):
        """For restriction based protocol data, obtain genomic location of cut site

        Returns:
            None if site is not available

            chromosome (str)
            position (int)
            strand (bool)
        """

        for fragment in self.fragments:
            site = fragment.get_site_location()
            if site is not None:
                return tuple( (*site, fragment.get_strand()))
        return None

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
        if fragment.strand is not None:
            self.strand = fragment.strand

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

    def can_be_yielded(self, chromosome, position):

        """Check if the molecule is far enough away from the supplied location to be ejected from a buffer.

        Args:
            chromosome (str) : chromosome / contig of location to test
            position (int) : genomic location of location to test

        Returns:
            can_be_yielded (bool)
        """

        if chromosome is None:
            return False
        if chromosome!=self.chromosome:
            return True
        return position < (self.spanStart-self.cache_size*0.5) or position > (self.spanEnd+self.cache_size*0.5)

    def get_rt_reaction_fragment_sizes(self):

        """Obtain all RT reaction fragment sizes
        Returns:
            rt_sizes (list of ints)
        """

        rt_reactions = molecule_to_random_primer_dict(self)
        amount_of_rt_reactions = len(rt_reactions)

        #this obtains the maximum fragment size:
        frag_chrom, frag_start, frag_end = pysamiterators.iterators.getListSpanningCoordinates([v for v in itertools.chain.from_iterable(self) if v is not None])

        #Obtain the fragment sizes of all RT reactions:
        rt_sizes = []
        for (rt_end,hexamer), fragment in rt_reactions.items():
            rt_chrom, rt_start, rt_end = pysamiterators.iterators.getListSpanningCoordinates(itertools.chain.from_iterable(fragment))
            rt_sizes.append([rt_end-rt_start])
        return rt_sizes

    def get_mean_rt_fragment_size(self):
        """Obtain the mean RT reaction fragment size

        Returns:
            mean_rt_size (float)
        """
        return np.nanmean(
            self.get_rt_reaction_fragment_sizes()
        )



    def __getitem__(self, index):
        """Obtain a fragment belonging to this molecule.

        Args:
            index (int):
                index of the fragment [0 ,1 , 2 ..]

        Returns:
            fragment (singlecellmultiomics.fragment.Fragment)
        """
        return self.fragments[key]



    def get_base_observation_dict(self, return_refbases=False):
        '''
        Obtain observed bases at reference aligned locations

        Args:
            return_refbases ( bool ):
                return both observed bases and reference bases

        Returns:
            { genome_location (tuple) : base (string) : obs (int) }
            and
            { genome_location (tuple) : base (string) if return_refbases is True }
        '''
        base_obs = collections.defaultdict(collections.Counter)
        if return_refbases:
            ref_bases = {}
        for fragment in self:
            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                allow_unsafe=(R1 is None))

            for read in [R1,R2]:
                if read is None:
                    continue
                for cycle, query_pos, ref_pos, ref_base in pysamiterators.iterators.ReadCycleIterator(read,with_seq=True):
                    if query_pos is None or ref_pos is None or ref_pos<start or ref_pos>end:
                        continue
                    query_base = read.seq[query_pos]
                    if query_base=='N':
                        continue
                    base_obs[(read.reference_name,ref_pos)][query_base]+=1
                    if return_refbases:
                        ref_bases[(read.reference_name,ref_pos)]=ref_base.upper()
        if return_refbases:
            return base_obs, ref_bases
        return base_obs



    def get_consensus(self, base_obs=None):
        """Get dictionary containing consensus calls in respect to reference

        Args:
            base_obs (collections.defaultdict(collections.Counter)) :
                { genome_location (tuple) : base (string) : obs (int) }

        Returns:
            consensus (dict)  :  {location : base}
        """
        consensus = {} # postion -> base , key is not set when not decided
        if base_obs is None:
            base_obs = self.get_base_observation_dict()

        for location, obs in base_obs.items():
            votes = obs.most_common()
            if len(votes)==1 or votes[1][1]<votes[0][1]:
                consensus[location] = votes[0][0]
        return consensus

    def check_variants(self, variants, exclude_other_calls=True ): #when enabled other calls (non ref non alt will be set None)
        """Verify variants in molecule

        Args:
            variants (pysam.VariantFile) : Variant file handle to extract variants from

        Returns:
            dict (collections.defaultdict( collections.Counter )) : { (chrom,pos) : ( call (str) ): observations  (int) }
        """
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
        """Iterate over all associated reads
        Returns:
            generator (pysam.AlignedSegment)
        """

        for fragment in self.fragments:
            for read in fragment:
                if read is not None:
                    yield read

    def __iter__(self):
        """Iterate over all associated fragments

        Yields:
            singlecellmultiomics.fragment.Fragment
        """
        for fragment in self.fragments:
            yield fragment


    def get_methylated_count(self, context=3):
        """Get the total amount of methylated bases
        Args:
            context (int) : 3 or 4 base context

        Returns:
            r (collections.Counter) : sum of methylated bases in contexts
        """

        r = collections.Counter()


    def get_TAPS_methylation_calls(self, capture_context=None, reference=None):
        """Get dictionary of TAPS methylated bases
        Args:
            capture_context (bool)

            reference (pysam.FastaFile)

        Returns:
            methylation_status (dict): { tuple(chrom,pos):  converted (bool) }

        """

        if self.is_multimapped():
            return None

        strand =  self.get_strand() #molecule.fragments[0][0].get_tag('RS')=='-'
        if strand is None:
            return None

        try:
            base_obs, ref_bases = self.get_base_observation_dict(return_refbases=True)
        except ValueError:
            # We cannot determine a reliable consensus sequence
            return None

        consensus = self.get_consensus(base_obs=base_obs)
        methylation_status = {} # location-> call
        for location, ref_base in ref_bases.items():
            if not location in consensus:
                continue
            query_base = consensus[location]
            k = (*location, strand)
            if ref_base=='C' and strand: # fwd
                if query_base=='T':
                    methylation_status[k] = True #  modified
                elif query_base=='C':
                    methylation_status[k] = False # not modified

            elif ref_base=='G' and not strand: # rev

                if query_base=='A':
                    methylation_status[k] = True #  modified
                elif query_base=='G':
                    methylation_status[k] = False # not modified

        return methylation_status


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


def MoleculeIterator( alignments, moleculeClass=Molecule, fragmentClass=Fragment, check_eject_every=1000, molecule_class_args={}, fragment_class_args={}, **pysamArgs):
    """Iterate over molecules in pysam.AlignmentFile

    Args:
        alignments (pysam.AlignmentFile): Alignments to extract molecules form

        moleculeClass (pysam.FastaFile): Class to use for molecules.

        fragmentClass (pysam.FastaFile): Class to use for fragments.

        check_eject_every (int): Check for yielding every N reads.

        molecule_class_args (dict): arguments to pass to moleculeClass.

        fragment_class_args (dict): arguments to pass to fragmentClass.

        **kwargs: arguments to pass to the pysam.AlignmentFile.fetch function

    Yields:
        methylation_status (dict): { tuple(chrom,pos):  converted (bool) }
    """

    molecules = []

    added_fragments = 0
    for R1,R2 in pysamiterators.iterators.MatePairIterator(alignments,performProperPairCheck=False,**pysamArgs):

        fragment = fragmentClass([R1,R2], **fragment_class_args)
        if not fragment.is_valid() :
            continue

        added = False
        for molecule in molecules:
            if molecule.add_fragment(fragment):
                added = True
                break
        if not added:
            molecules.append(moleculeClass(fragment, **molecule_class_args ))

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
    yield from iter(molecules)
