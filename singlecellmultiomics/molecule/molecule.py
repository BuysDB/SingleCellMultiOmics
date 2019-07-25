from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
import singlecellmultiomics.universalBamTagger
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

        reference (pysam.FastaFile) : reference file, used to obtain base contexts and correct aligned_pairs iteration when the MD tag is not correct

        strand (bool): mapping strand.
            True when strand is REVERSE
            False when strand is FORWARD
            None when strand is not determined
    """

    def __init__(self, fragments=None, cache_size=10_000,reference=None):
        """Initialise Molecule

        Parameters
        ----------
        fragments :  list(singlecellmultiomics.fragment.Fragment)
            Fragment to assign to Molecule. More fragments can be added later

        cache_size (int): radius of molecule assignment cache
        """

        self.reference = reference
        self.fragments  = []
        self.spanStart = None
        self.spanEnd = None
        self.chromosome = None
        self.cache_size = cache_size
        self.strand = None
        self.umi = None
        self.umi_hamming_distance = None
        self.fragment_match = None # when set, when comparing to a fragment the fragment to be added has to match this hash
        self.umi_counter = collections.Counter() # Observations of umis
        if fragments is not None:
            if type(fragments) is list:
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)


    def __len__(self):
        return len(self.fragments)

    def get_consensus_base_frequencies(self):
        """Obtain the frequency of bases in the molecule consensus sequence

        Returns:
            base_frequencies (collections.Counter) : Counter containing base frequecies, for example: { 'A':10,'T':3, C:4 }
        """
        return collections.Counter( self.get_consensus().values() )

    def get_consensus_gc_ratio(self):
        """Obtain the GC ratio of the molecule consensus sequence

        Returns:
            gc_ratio(float) : GC ratio
        """
        bf = self.get_consensus_base_frequencies()
        return (bf['G']+bf['C'])/sum(bf.values())

    def get_umi_error_rate(self):
        """Obtain fraction of fragments that are associated
        to the molecule with a exact matching UMI vs total amount of associated fragments
        Returns:
            exact_matching_fraction (float)
        """
        mc = 0
        other = 0
        for i,(umi,obs) in enumerate( m.umi_counter.most_common() ):
            if i==0:
                mc = obs
            else:
                other += obs

        return mc/(other+mc)

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

    def update_umi(self):
        """Set UMI
        Sets:
            self.umi (str):
                Returns the most common umi associated to the molecule
        """
        self.umi = self.umi_counter.most_common(1)[0][0]

    def get_umi(self):
        return self.umi

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

    def get_mean_mapping_qual(self):
        """Get mean mapping quality of the molecule

        Returns:
            mean_mapping_qual (float)
        """
        return np.mean( [fragment.mapping_quality for fragment in self] )

    def is_multimapped(self):
        """Check if the molecule is multimapping

        Returns:
            is_multimapped (bool) : True when multimapping
        """
        for fragment in self.fragments:
            if not fragment.is_multimapped:
                return False
        return True

    def _add_fragment(self, fragment):
        self.match_hash = fragment.match_hash
        self.fragments.append(fragment)
        # Update span:
        add_span = fragment.get_span()
        self.spanStart = add_span[1] if self.spanStart is None else min(add_span[1], self.spanStart)
        self.spanEnd = add_span[2] if self.spanEnd is None else max(add_span[2], self.spanEnd)
        self.chromosome = add_span[0]
        self.span = (self.chromosome, self.spanStart, self.spanEnd)
        if fragment.strand is not None:
            self.strand = fragment.strand
        self.umi_counter[fragment.umi]+=1
        self.umi_hamming_distance = fragment.umi_hamming_distance
        self.saved_base_obs = None
        self.update_umi()

    def add_fragment(self, fragment):
        """Associate a fragment with this Molecule

        Args:
            fragment (singlecellmultiomics.fragment.Fragment) : Fragment to associate
        Returns:
            has_been_added (bool) : Returns False when the fragments which have already been associated to the molecule refuse the fragment
        """
        if len(self.fragments)==0:
            self._add_fragment(fragment)
            self.sample=fragment.sample
            return True

        #for f in self.fragments:
        #    if f == fragment:
        #        # it matches this molecule:
        #        self._add_fragment(fragment)
        #        return True
        if self==fragment:
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

    def get_rt_reactions(self):
        return molecule_to_random_primer_dict(self)

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


    def set_methylation_call_tags(self,
                              call_dict, bismark_call_tag='XM',
                              total_methylated_tag='MC',
                              total_unmethylated_tag='uC',
                             ):
        """Set methylation call tags given a methylation dictionary

        This method sets multiple tags in every read associated to the molecule.
        The tags being set are the bismark_call_tag, every aligned base is annotated
        with a zZxXhH or ".", and a tag for both the total methylated C's and unmethylated C's

        Args:
            call_dict (dict) : Dictionary containing bismark calls (chrom,pos) : letter

            bismark_call_tag (str) : tag to write bismark call string

            total_methylated_tag (str) : tag to write total methylated bases

            total_unmethylated_tag (str) : tag to write total unmethylated bases

        Returns:
            can_be_yielded (bool)
        """

        # Contruct XM strings
        for fragment in self:
            for read in fragment:
                if read is None:
                    continue
                read.set_tag(
                    # Write the methylation tag to the read
                    bismark_call_tag,
                    ''.join([
                        call_dict.get(
                            (read.reference_name,rpos),'.' )  # Obtain all aligned positions from the call dict
                        for qpos, rpos in read.get_aligned_pairs(matches_only=True) # iterate all positions in the alignment
                        if qpos is not None and rpos is not None]) # make sure to ignore non matching positions ? is this neccesary?
                )

                # Set total methylated bases
                read.set_tag(
                    total_methylated_tag,
                    sum( x.isupper() for x in read.get_tag('XM') )
                )

                # Set total unmethylated bases
                read.set_tag(
                    total_unmethylated_tag,
                    sum( x.islower() for x in read.get_tag('XM') )
                )


    def __getitem__(self, index):
        """Obtain a fragment belonging to this molecule.

        Args:
            index (int):
                index of the fragment [0 ,1 , 2 ..]

        Returns:
            fragment (singlecellmultiomics.fragment.Fragment)
        """
        return self.fragments[index]


    def get_mean_base_quality(self, chromosome, position, base):
        """Get the mean phred score at the supplied coordinate and base-call

        Args:
            chromosome (str)
            position (int)
            base (str)

        Returns:
            mean_phred_score (float)
        """
        for fragment in self:
            if fragment.span[0]!=chromosome:
                continue
            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            try:
                start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                allow_unsafe=(R1 is None))
            except ValueError as e:
                continue

            qualities = []
            for read in (R1,R2):
                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                    read,with_seq=False):
                    if query_pos is None or ref_pos != position or read.seq[query_pos]!=base:
                        continue

                    qualities.append( ord( read.qual[query_pos] ) )
        if len(qualities)==0:
            raise IndexError("There are no observations if the supplied base/location combination")
        return np.mean(qualities)

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

        # Check if cached is available
        if self.saved_base_obs is not None:
            if not return_refbases:
                return self.saved_base_obs[0]
            else:
                if self.saved_base_obs[1] is not None:
                    return self.saved_base_obs

        base_obs = collections.defaultdict(collections.Counter)

        if return_refbases:
            ref_bases = {}
        used = 0 #  some alignments yielded valid calls
        ignored = 0
        for fragment in self:
            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            try:
                start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                allow_unsafe=(R1 is None))
            except ValueError as e:
                ignored+=1
                continue
            used+=1
            for read in [R1,R2]:
                if read is None:
                    continue
                for cycle, query_pos, ref_pos, ref_base in pysamiterators.iterators.ReadCycleIterator(
                    read,with_seq=True,reference=self.reference):

                    if query_pos is None or ref_pos is None or ref_pos<start or ref_pos>end:
                        continue
                    query_base = read.seq[query_pos]
                    if query_base=='N':
                        continue
                    base_obs[(read.reference_name,ref_pos)][query_base]+=1

                    if return_refbases:
                        ref_bases[(read.reference_name,ref_pos)]=ref_base.upper()

        if used==0 and ignored>0:
            raise ValueError('Could not extract any safe data from molecule')

        self.saved_base_obs = (base_obs, ref_bases)

        if return_refbases:
            return base_obs, ref_bases

        return base_obs

    def get_match_mismatch_frequency(self, ignore_locations=None):
        """Get amount of base-calls matching and mismatching the reference sequence,
           mismatches in every read are counted

        Args:
            ignore_locations (iterable(tuple([chrom(str),pos(int)])) ) :
                Locations not to take into account for the match and mismatch frequency

        Returns:
            matches(int), mismatches(int)
        """
        matches = 0
        mismatches = 0

        base_obs, ref_bases = self.get_base_observation_dict(return_refbases=True)
        for location, obs in base_obs.items():
            if ignore_locations is not None and location in ignore_locations:
                continue

            if location in ref_bases:
                ref = ref_bases[location]
                if ref not in 'ACTG': # don't count weird bases in the reference @warn
                    continue
                matches += obs[ref]
                mismatches += sum((base_obs for base,base_obs in obs.most_common() if base != ref))


        return matches, mismatches

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
            try:
                base_obs, ref_bases = self.get_base_observation_dict(return_refbases=True)
            except ValueError as e:
                # We cannot determine safe regions
                raise



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


    def get_html(self,reference=None):
        """Get html representation of the molecule

        Returns:
            html_rep(str) : Html representation of the molecule
        """

        span_len = self.spanEnd-self.spanStart
        if span_len > 1000:
            raise ValueError('The molecule is too long to display')

        visualized = ['.']  * span_len
        for location,base in self.get_consensus().items():
            visualized[location-self.spanStart] = base
        return visualized


    def get_methylation_dict(self):
        """Obtain methylation dictionary

        Returns:
            methylated_positions (collections.Counter):
                (read.reference_name, rpos) : times seen methylated

            methylated_state (dict):
                {(read.reference_name, rpos) : 1/0/-1 }
                1 for methylated
                0 for unmethylated
                -1 for unknown

        """
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




class MoleculeIterator():

    def __init__(self, alignments, moleculeClass=Molecule,
        fragmentClass=Fragment,
        check_eject_every=10_000, #bigger sizes are very speed benificial
        # because the relative amount of molecules which can be ejected will be much higher
        molecule_class_args={},
        fragment_class_args={},
        perform_qflag=True,
        pooling_method=1,
        **pysamArgs):
        """Iterate over molecules in pysam.AlignmentFile

        Args:
            alignments (pysam.AlignmentFile): Alignments to extract molecules form

            moleculeClass (pysam.FastaFile): Class to use for molecules.

            fragmentClass (pysam.FastaFile): Class to use for fragments.

            check_eject_every (int): Check for yielding every N reads.

            molecule_class_args (dict): arguments to pass to moleculeClass.

            fragment_class_args (dict): arguments to pass to fragmentClass.

            perform_qflag (bool):  Make sure the sample/umi etc tags are copied
                from the read name into bam tags

            pooling_method(int) : 0: no  pooling, 1: only compare molecules with the same sample id.
            **kwargs: arguments to pass to the pysam.AlignmentFile.fetch function

        Yields:
            molecule (Molecule): Molecule
        """
        self.alignments = alignments
        self.moleculeClass = moleculeClass
        self.fragmentClass = fragmentClass
        self.check_eject_every = check_eject_every
        self.molecule_class_args = molecule_class_args
        self.fragment_class_args = fragment_class_args
        self.perform_qflag = perform_qflag
        self.pysamArgs = pysamArgs
        self.matePairIterator=None
        self.pooling_method=pooling_method
        self._clear_cache()

    def _clear_cache(self):
        """Clear cache containing non yielded molecules"""
        self.waiting_fragments = 0
        self.yielded_fragments = 0
        self.yielded_molecules = 0
        self.check_ejection_iter = 0
        if self.pooling_method==0:
            self.molecules=[]
        elif self.pooling_method==1:
            self.molecules_per_cell =collections.defaultdict(list)
        else:
            raise NotImplementedError()

    def __repr__(self):
        return f"""Molecule Iterator, generates fragments from {self.fragmentClass} into molecules based on {self.moleculeClass}.
        Yielded {self.yielded_fragments} fragments, {self.waiting_fragments} fragments are waiting to be ejected.
        {self.get_molecule_cache_size()} molecules cached.
        Mate pair iterator: {str(self.matePairIterator)}"""


    def get_molecule_cache_size(self):
        if self.pooling_method==0:
            return len(self.molecules)
        elif self.pooling_method==1:
            return sum( len(cell_molecules) for cell, cell_molecules in self.molecules_per_cell.items() )

        else:
            raise NotImplementedError()

    def __iter__(self):
        if self.perform_qflag:
            qf = singlecellmultiomics.universalBamTagger.QueryNameFlagger()

        self._clear_cache()

        self.waiting_fragments = 0
        self.matePairIterator = pysamiterators.iterators.MatePairIterator(self.alignments,performProperPairCheck=False,**self.pysamArgs)
        for R1,R2 in self.matePairIterator:
            # Make sure the sample/umi etc tags are placed:
            if self.perform_qflag:
                qf.digest([R1,R2])

            fragment = self.fragmentClass([R1,R2], **self.fragment_class_args)
            if not fragment.is_valid() :
                continue

            added = False
            if self.pooling_method==0:
                for molecule in molecules:
                    if molecule.add_fragment(fragment):
                        added = True
                        break
            elif self.pooling_method==1:
                for molecule in self.molecules_per_cell[fragment.sample]:
                    if molecule.add_fragment(fragment):
                        added = True
                        break

            if not added:
                if self.pooling_method==0:
                    molecules.append(self.moleculeClass(fragment, **self.molecule_class_args ))
                else:
                    self.molecules_per_cell[fragment.sample].append(
                        self.moleculeClass(fragment, **self.molecule_class_args )
                    )

            self.waiting_fragments+=1
            self.check_ejection_iter += 1
            if self.check_ejection_iter>self.check_eject_every:
                current_chrom, _, current_position = fragment.get_span()
                to_pop = []
                self.check_ejection_iter=0
                if self.pooling_method==0:
                    for i,m in enumerate(molecules):
                        if m.can_be_yielded(current_chrom,current_position):
                            to_pop.append(i)
                            self.waiting_fragments-=len(m)
                            self.yielded_fragments+=len(m)

                    for i,j in enumerate(to_pop):
                        yield molecules.pop(i-j)
                else:
                    for cell, cell_molecules in self.molecules_per_cell.items():
                        for i,m in enumerate(cell_molecules):
                            if m.can_be_yielded(current_chrom,current_position):
                                to_pop.append((cell,i))
                                self.waiting_fragments-=len(m)
                                self.yielded_fragments+=len(m)
                    for i,(cell, index) in enumerate(to_pop):
                        yield molecules[cell][index]
                        del molecules[cell][index]

        # Yield remains
        if self.pooling_method==0:
            yield from iter(molecules)
        else:
            for cell, cell_molecules in self.molecules_per_cell.items():
                for i,m in enumerate(cell_molecules):
                    yield m
        self._clear_cache()
