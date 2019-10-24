from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
import singlecellmultiomics.universalBamTagger
from singlecellmultiomics.fragment import Fragment
import collections
import itertools
import numpy as np
from singlecellmultiomics.utils import style_str
import textwrap
import  singlecellmultiomics.alleleTools
import functools
import typing
import pysam
import pysamiterators
from singlecellmultiomics.utils import find_ranges, create_MD_tag


@functools.lru_cache(maxsize=1000)
def might_be_variant(chrom,pos, variants, ref_base=None):
    """Returns True if a variant exists at the given coordinate"""
    if ref_base=='N':
        return False
    for record in  variants.fetch(chrom,pos,pos+1):
        return True
    return False

def molecule_to_random_primer_dict(molecule, primer_length=6, primer_read=2, max_N_distance=0): #1: read1 2: read2
    rp = collections.defaultdict(list)

    # First add all reactions without a N in the sequence:
    for fragment in molecule:

        hstart, hseq = fragment.get_random_primer_hash()
        if hseq == None:
            # This should really not happen with freshly demultiplexed data, it means we cannot extract the random primer sequence
            # which should be present as a tag (rP) in the record
            rp[None, None].append(fragment)
        elif not 'N' in hseq:
            rp[hstart, hseq].append(fragment)

    # Try to match reactions with N with reactions without a N
    for fragment in molecule:

        hstart, hseq = fragment.get_random_primer_hash()
        if hseq is not None and 'N' in hseq:
            # find nearest
            for other_start, other_seq in rp:
                if other_start!=hstart:
                    continue

                if hseq.count('N')>max_N_distance:
                    continue

                if 'N' in other_seq:
                    continue

                if hamming_distance(hseq,other_seq)==0:
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

    def __init__(self,
        fragments : typing.Optional[typing.Iterable] = None,
        cache_size : int = 10_000,
        reference : typing.Union[pysam.FastaFile, pysamiterators.CachedFasta] = None,
        min_max_mapping_quality :  typing.Optional[int] = None,# When all fragments have a mappin quality below this value the is_valid method will return False
        allele_resolver : typing.Optional[singlecellmultiomics.alleleTools.AlleleResolver] = None ,
        **kwargs):
        """Initialise Molecule

        Parameters
        ----------
        fragments :  list(singlecellmultiomics.fragment.Fragment)
            Fragment to assign to Molecule. More fragments can be added later

        min_max_mapping_quality :  When all fragments have a mappin quality below this value the is_valid method will return False

        allele_resolver :  alleleTools.AlleleResolver or None. Supply an allele resolver in order to assign an allele to the molecule

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
        self.min_max_mapping_quality = min_max_mapping_quality
        self.umi_counter = collections.Counter() # Observations of umis
        if fragments is not None:
            if type(fragments) is list:
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)

        # Obtain allele if available
        self.allele_resolver = allele_resolver
        self.allele = None
        if self.allele_resolver is not None:
            try:
                hits = self.get_allele(allele_resolver)
                # Only store when we have a unique single hit:
                if len(hits)==1:
                    self.allele = list(hits)[0]
            except ValueError as e:
                # This happens when a consensus can not be obtained
                pass

    def get_a_reference_id(self):
        for read in self.iter_reads():
            if not read.is_unmapped:
                return read.reference_id
        return -1

    def get_consensus_read(self, target_file,
            read_name,consensus=None,
            phred_scores=None,
            cigarstring=None,
            mdstring=None,
            start=None,
            supplementary=False

            ):
        """get pysam.AlignedSegment containing aggregated molecule information

        Args:
            target_file(pysam.AlignmentFile) : File to create the read for

            read_name(str) : name of the read to write
        Returns:
            read(pysam.AlignedSegment)
        """
        if start is None:
            start = self.spanStart
        if consensus is None:
            try: # Obtain the consensus sequence
                consensus = self.get_consensus()
            except Exception as e:
                raise
        if type(consensus)==str:
            sequence = consensus
        else:
            sequence = ''.join(( consensus.get((self.chromosome, ref_pos),'N')
             for ref_pos in range(self.spanStart, self.spanEnd+1)
            ))

        if type(phred_scores)==dict:
            phred_score_array = list( phred_scores.get((self.chromosome, ref_pos),0)
             for ref_pos in range(self.spanStart, self.spanEnd+1)
            )
        else:
            phred_score_array = phred_scores

        # Construct consensus - read
        cread = pysam.AlignedSegment(header=target_file.header)
        cread.reference_name = self.chromosome
        cread.reference_start = start
        cread.query_name = read_name
        cread.query_sequence = sequence
        cread.query_qualities = phred_score_array
        cread.is_supplementary = supplementary
        if cigarstring is not None:
            cread.cigarstring  = cigarstring
        else:
            cread.cigarstring  = f'{len(sequence)}M'
        cread.mapping_quality = self.get_max_mapping_qual()
        cread.is_reverse = self.strand
        if mdstring is not None:
            cread.set_tag('MD',mdstring)

        self.write_tags_to_psuedoreads( (cread,) )

        return cread

    """ method = 1
        sequence = []
        cigar = []
        if method==0:
            prev_end = None
            for block_start,block_end in molecule.get_aligned_blocks():
                if molecule.strand:
                    print(block_end>block_start,block_start, block_end)
                if prev_end is not None:
                    cigar.append(f'{block_start - prev_end}D')

                block_len = block_end-block_start+1
                cigar.append(f'{block_len}M')
                for ref_pos in range(block_start,block_end+1):
                    call = consensus.get((molecule.chromosome, ref_pos),'N')
                    sequence.append(call)
                prev_end = block_end+1

            cigarstring = ''.join(cigar)
        """

    def get_feature_vector(self, window_size=90):
        """ Obtain a feature vector representation of the molecule

        Returns:
            feature_vector(np.array)
        """

        return np.array([
            self.get_strand(),
            self.has_valid_span(),
            self.get_umi_error_rate(),
            self.get_consensus_gc_ratio(),
            len(get_raw_barcode_sequences),
            self.get_safely_aligned_length(),
            self.get_max_mapping_qual(),
            (self.alleles is None),
            self.contains_valid_fragment(),
            self.is_multimapped(),
            self.get_feature_window(window_size=window_size)
        ])

    def get_tag_counter(self):
        """
        Obtain a dictionary with tag -> value -> frequency

        Returns:
            tag_obs (collections.defaultdict(collections.Counter)):
                { tag(str) : { value(int/str): frequency:(int) }


        """
        tags_obs = collections.defaultdict(collections.Counter)
        for tag, value in itertools.chain(*[r.tags for r in self.iter_reads()]) :
            tags_obs[tag][value] += 1
        return tag_obs

    def write_tags_to_psuedoreads(self,reads):
        """
        Write molecule information to the supplied reads as BAM tags
        """
        # write methylation tags to new reads if applicable:
        if self.methylation_call_dict is not None:
            self.set_methylation_call_tags(self.methylation_call_dict, reads=reads)

        for read in reads:
            read.set_tag('SM', self.sample)
            if hasattr(self, 'get_cut_site'):
                read.set_tag('DS', self.get_cut_site()[1])

            if self.umi is not None:
                read.set_tag('RX',self.umi)
                bc = list(self.get_barcode_sequences())[0]
                read.set_tag('BC',bc)
                read.set_tag('MI',bc+self.umi)

            # Store total amount of RT reactions:
            read.set_tag('TR',len(self.get_rt_reactions()))

            if self.allele is not None:
                read.set_tag('DA', self.allele)

        if self.allele_resolver is not None:
            self.write_allele_phasing_information_tag(self.allele_resolver,reads=reads )



    def deduplicate_to_single(self, target_bam, read_name, classifier,reference=None):
        """
        Deduplicate all reads associated to this molecule to a single pseudoread

        Args:
            target_bam (pysam.AlignmentFile) : file to associate the read with
            read_name (str) : name of the pseudoread
            classifier (sklearn classifier) : classifier for consensus prediction

        Returns:
            read (pysam.AlignedSegment) : Pseudo-read containing aggregated information
        """
        # Set all associated reads to duplicate
        for read in self.iter_reads():
            read.is_duplicate = True

        features = self.get_base_calling_feature_matrix(reference=reference)
        predicted_sequence = classifier.predict(features)
        predicted_sequence[ features[:, [ x*8 for x in range(4) ] ].sum(1)==0 ] ='N'
        phred_scores = np.rint(
                -10*np.log10( np.clip(1-classifier.predict_proba(features).max(1), 0.000000001, 0.999999999 )
            )).astype('B')

        read = self.get_consensus_read(
                    read_name=read_name,
                    target_file = target_bam,
                    consensus=''.join(predicted_sequence),
                    phred_scores=phred_scores)
        read.is_read1 = True
        return read

    def deduplicate_to_single_CIGAR_spaced(self, target_bam, read_name, classifier, max_N_span = 300,reference=None ):
        """
        Deduplicate all associated reads to a single pseudoread, when the span is larger than max_N_span
        the read is split up in multi-segments. Uncovered locations are spaced using N's in the CIGAR.

        Args:
            target_bam (pysam.AlignmentFile) : file to associate the read with
            read_name (str) : name of the pseudoread
            classifier (sklearn classifier) : classifier for consensus prediction
        Returns:
            reads( list [ pysam.AlignedSegment ] )
        """
        # Set all associated reads to duplicate
        for read in self.iter_reads():
            read.is_duplicate = True

        features, reference_bases, CIGAR, alignment_start, alignment_end = self.get_base_calling_feature_matrix_spaced(
                True,
                reference=reference)

        predicted_sequence =  classifier.predict(features)
        reference_sequence = ''.join([base for chrom, pos, base  in reference_bases])
        predicted_sequence[ features[:, [ x*8 for x in range(4) ] ].sum(1)==0 ] ='N'
        predicted_sequence = ''.join( predicted_sequence )

        phred_scores = np.rint(
                -10*np.log10( np.clip(1-classifier.predict_proba(features).max(1),
                                      0.000000001,
                                      0.999999999 )
            )).astype('B')

        reads = []

        query_index_start = 0
        query_index_end = 0
        reference_position = alignment_start # pointer to current position
        reference_start = alignment_start # pointer to alignment start of current read
        supplementary = False
        partial_CIGAR = []
        partial_MD = []

        for operation, amount in CIGAR:
            if operation=='M': # Consume query and reference
                query_index_end+=amount
                reference_position+=amount
                partial_CIGAR.append(f'{amount}{operation}')

            if operation=='N':
                # Consume reference:
                reference_position+=amount
                if amount>max_N_span: # Split up in supplementary alignment
                    # Eject previous
                    #reference_seq =

                    consensus_read = self.get_consensus_read(
                                read_name=read_name,
                                target_file = target_bam,
                                consensus=predicted_sequence[query_index_start:query_index_end],
                                phred_scores=phred_scores[query_index_start:query_index_end],
                                cigarstring=''.join(partial_CIGAR),
                                mdstring = create_MD_tag(
                                    reference_sequence[query_index_start:query_index_end],
                                    predicted_sequence[query_index_start:query_index_end]
                                ),
                                start = reference_start,
                                supplementary=supplementary
                    )
                    reads.append( consensus_read )
                    if not supplementary:
                        consensus_read.is_read1 = True

                    supplementary= True
                    # Start new:
                    query_index_start = query_index_end
                    reference_start = reference_position
                    partial_CIGAR = []
                else:
                    partial_CIGAR.append(f'{amount}{operation}')


        reads.append( self.get_consensus_read(
                    read_name=read_name,
                    target_file = target_bam,
                    consensus=''.join(predicted_sequence[query_index_start:query_index_end]),
                    phred_scores=phred_scores[query_index_start:query_index_end],
                    cigarstring=''.join(partial_CIGAR),
                    mdstring = create_MD_tag(
                                        reference_sequence[query_index_start:query_index_end],
                                        predicted_sequence[query_index_start:query_index_end]

                                ),
                    start = reference_start,
                    supplementary=supplementary
        ))

        # Write last index tag to last read ..
        if supplementary:
            reads[-1].is_read2 = True

        # Write NH tag (the amount of records with the same query read):
        for read in reads:
            read.set_tag('NH', len(reads))

        return reads

    def get_base_calling_feature_matrix(self, return_ref_info=False, start=None, end=None, reference=None, NUC_RADIUS = 1, USE_RT=True, select_read_groups=None):
        """
        Obtain feature matrix for base calling

        Args:
            return_ref_info (bool) : return both X and array with feature information
            start (int) : start of range, genomic position
            end (int) : end of range (inclusive), genomic position
            reference(pysam.FastaFile) : reference to fetch reference bases from, if not supplied the MD tag is used
            NUC_RADIUS(int) : generate kmer features target nucleotide
            USE_RT(bool) : use RT reaction features
            select_read_groups(set) : only use reads from these read groups to generate features
        """
        if start is None:
            start = self.spanStart
        if end is None:
            end = self.spanEnd

        with np.errstate(divide='ignore', invalid='ignore'):
            BASE_COUNT = 5
            RT_INDEX = 7 if USE_RT else None
            STRAND_INDEX = 0
            PHRED_INDEX = 1
            RC_INDEX = 2
            ALIGNED_WP_INDEX = 3
            CYCLE_INDEX = 4
            MQ_INDEX = 5
            FS_INDEX = 6

            COLUMN_OFFSET = 0
            features_per_block = 8 - (not USE_RT)

            origin_start = start
            origin_end = end

            end += NUC_RADIUS
            start -= NUC_RADIUS

            features = np.zeros( (end - start + 1, (features_per_block*BASE_COUNT) + COLUMN_OFFSET ) )

            if return_ref_info:
                ref_bases = {}

            for rt_id,fragments in self.get_rt_reactions().items():
                # we need to keep track what positions where covered by this RT reaction
                RT_reaction_coverage = set() # (pos, base_call)
                for fragment in fragments:
                    for read in fragment:
                        if select_read_groups is not None:
                            if not read.has_tag('RG'):
                                raise ValueError("Not all reads in the BAM file have a read group defined.")
                            if not read.get_tag('RG') in select_read_groups:
                                continue
                        # Skip reads outside range
                        if read is None or read.reference_start > (end+1) or read.reference_end < start:
                            continue
                        for cycle, q_pos, ref_pos, ref_base in  pysamiterators.ReadCycleIterator(read, matches_only=True,with_seq=True, reference=reference):

                            row_index = ref_pos-start
                            if row_index<0 or row_index>=features.shape[0]:
                                continue

                            query_base = read.seq[q_pos]
                            # Base index block:
                            block_index = 'ACGTN'.index(query_base)

                            # Update rt_reactions
                            if USE_RT:
                                if not (ref_pos,query_base) in RT_reaction_coverage:
                                    features[row_index][RT_INDEX + COLUMN_OFFSET +features_per_block*block_index] += 1
                                RT_reaction_coverage.add( (ref_pos,query_base) )

                            # Update total phred score
                            features[row_index][PHRED_INDEX + COLUMN_OFFSET +features_per_block*block_index] += read.query_qualities[q_pos]

                            # Update total reads

                            features[row_index][RC_INDEX + COLUMN_OFFSET +features_per_block*block_index] += 1


                            # Update primer mp
                            if fragment.safe_span:
                                features[row_index][ALIGNED_WP_INDEX + COLUMN_OFFSET +features_per_block*block_index] +=  (ref_pos<fragment.span[1] or ref_pos>fragment.span[2] )
                            else:
                                features[row_index][ALIGNED_WP_INDEX + COLUMN_OFFSET +features_per_block*block_index] +=  (ref_pos<fragment.span[1] or ref_pos>fragment.span[2] )
                                #fragment_sizes[key].append( abs( fragment.span[1] - fragment.span[2] ) )

                            # Update fragment sizes:
                            features[row_index][FS_INDEX + COLUMN_OFFSET +features_per_block*block_index] += abs( fragment.span[1] - fragment.span[2] )

                            # Update cycle
                            features[row_index][CYCLE_INDEX + COLUMN_OFFSET +features_per_block*block_index] += cycle

                            # Update MQ:
                            features[row_index][MQ_INDEX + COLUMN_OFFSET +features_per_block*block_index] += read.mapping_quality

                            # update strand:
                            features[row_index][STRAND_INDEX + COLUMN_OFFSET +features_per_block*block_index] += read.is_reverse

                            if return_ref_info:
                                row_index_in_output = ref_pos-origin_start
                                if row_index_in_output<0 or row_index_in_output>=origin_end-origin_start+1:
                                    continue

                                ref_bases[ref_pos] = ref_base.upper()

            # Normalize all and return

            for block_index in range(BASE_COUNT): #ACGTN
                for index in (PHRED_INDEX, ALIGNED_WP_INDEX, CYCLE_INDEX, MQ_INDEX, FS_INDEX, STRAND_INDEX  ):
                    features[:,index + COLUMN_OFFSET +features_per_block*block_index] /=  features[:,RC_INDEX + COLUMN_OFFSET +features_per_block*block_index]
            #np.nan_to_num( features, nan=-1, copy=False )
            features[np.isnan(features)] = -1

            if NUC_RADIUS>0:
                # duplicate columns in shifted manner
                x = features
                features = np.zeros( (x.shape[0]-NUC_RADIUS*2, x.shape[1]*(1+NUC_RADIUS*2)) )
                for offset in range(0,NUC_RADIUS*2+1):
                    slice_start = offset
                    slice_end = -(NUC_RADIUS*2)  + offset
                    if slice_end == 0:
                        features[:,features_per_block*BASE_COUNT*offset:features_per_block*BASE_COUNT*(offset+1)] = x[slice_start:,:]
                    else:
                        features[:,features_per_block*BASE_COUNT*offset:features_per_block*BASE_COUNT*(offset+1)] = x[slice_start:slice_end,:]

            if return_ref_info:
                ref_info = [
                    (self.chromosome, ref_pos, ref_bases.get(ref_pos,'N'))
                    for ref_pos in range(origin_start, origin_end+1)]
                return  features, ref_info
            return features

    @functools.lru_cache(maxsize=4)
    def get_base_calling_feature_matrix_spaced(self,return_ref_info=False, reference=None, **feature_matrix_args):
        """
        Obtain a base-calling feature matrix for all reference aligned bases.

        Returns:
            X : feature matrix
            y : reference bases
            CIGAR : alignment of feature matrix to reference tuples (operation, count)
            reference(pysam.FastaFile) : reference to fetch reference bases from, if not supplied the MD tag is used
        """

        X = None
        if return_ref_info:
            y = []
        CIGAR = []
        prev_end = None
        alignment_start = None
        alignment_end = None
        for start,end in self.get_aligned_blocks():
            if return_ref_info:
                x,y_ = self.get_base_calling_feature_matrix(
                        return_ref_info=return_ref_info, start=start, end=end,
                        reference=reference, **feature_matrix_args
                        )
                y+=y_
            else:
                x = self.get_base_calling_feature_matrix(
                        return_ref_info=return_ref_info, start=start, end=end,reference=reference, **feature_matrix_args
                        )
            if X is None:
                X = x
            else:
                X = np.append(X,x,axis=0)

            if prev_end is not None:
                CIGAR.append( ('N', start-prev_end-1) )
            CIGAR.append( ( 'M', (end-start+1) ) )
            prev_end = end

            if alignment_start is None:
                alignment_start = start
                alignment_end = end
            else:
                alignment_start=min(alignment_start,start)
                alignment_end=max(alignment_end,end)

        if return_ref_info:
            return X,y,CIGAR,alignment_start, alignment_end
        else:
            return X,CIGAR,alignment_start, alignment_end

    def get_base_calling_training_data(self,mask_variants=None,might_be_variant_function=None,reference=None, **feature_matrix_args):
        if mask_variants is not None and  might_be_variant_function is None:
            might_be_variant_function = might_be_variant

        features, feature_info, _CIGAR, _alignment_start, _alignment_end  = self.get_base_calling_feature_matrix_spaced(True,reference=reference, **feature_matrix_args)
        # check which bases should not be used
        use_indices = [
            mask_variants is None or
            not might_be_variant_function(chrom,pos, mask_variants, base)
            for chrom, pos, base in feature_info ]

        X_molecule = features[use_indices]
        y_molecule = [
            base for use,(chrom, pos, base) in
            zip(use_indices,feature_info) if use
            ]
        return X_molecule, y_molecule

    def has_valid_span(self):
        """Check if the span of the molecule is determined

        Returns:
            has_valid_span (bool)
        """
        if self.spanStart is not None and self.spanEnd is not None:
            return True
        return False

    def get_strand_repr(self, unknown='?'):
        """Get string representation of mapping strand
        Args:
            unknown (str) :  set what character/string to return
                             when the strand is not available

        Returns:
            strand_repr (str) : + forward, - reverse, ? unknown
        """
        s = self.get_strand()
        if s is None:
            return unknown
        if s:
            return '-'
        else:
            return '+'

    def write_tags(self):
        """ Write BAM tags to all reads associated to this molecule

        This function sets the following tags:
            - mI : most common umi
            - DA : allele
            - af : amount of associated fragments
            - rt : rt_reaction_index
            - rd : rt_duplicate_index
            - ap : phasing information (if allele_resolver is set)
        """
        self.is_valid(set_rejection_reasons=True)
        if  self.umi is not None:
            self.set_meta('mI', self.umi)
        if self.allele is not None:
            self.set_meta('DA', str(self.allele))

        # associatedFragmentCount :
        self.set_meta('af', len(self))

        # Write RT reaction tags (rt: rt reaction index, rd rt duplicate index)
        for rt_reaction_index,(_,frags) in enumerate(self.get_rt_reactions().items()):
            for rt_duplicate_index,frag in enumerate(frags):
                frag.set_meta('rt', rt_reaction_index)
                frag.set_meta('rd', rt_duplicate_index)

        if self.allele_resolver is not None:
            self.write_allele_phasing_information_tag()

    def set_rejection_reason(self,reason):
        """ Add rejection reason to all fragments associated to this molecule

        Args:
            reason (str) : rejection reason to set
        """
        for fragment in self:
            fragment.set_rejection_reason(reason)

    def is_valid(self, set_rejection_reasons=False):

        if self.is_multimapped():
            if set_rejection_reasons:
                self.set_rejection_reason('multimapping')
            return False

        if self.min_max_mapping_quality is not None and \
            self.get_max_mapping_qual()<self.min_max_mapping_quality:
            if set_rejection_reasons:
                self.set_rejection_reason('MQ')
            return False

        if not self.contains_valid_fragment():
            if set_rejection_reasons:
                self.set_rejection_reason('invalid_fragments')
            return False

        return True


    def get_aligned_blocks(self):
        """ get all consecutive blocks of aligned reference positions

        Returns:
            sorted list of aligned blocks (list) : [ (start, end), (start, end) ]
        """
        return find_ranges(
            sorted(list(set(
                (ref_pos
                for read in self.iter_reads()
                for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False) ))))
            )

    def __len__(self):
        return len(self.fragments)

    def get_consensus_base_frequencies(self):
        """Obtain the frequency of bases in the molecule consensus sequence

        Returns:
            base_frequencies (collections.Counter) : Counter containing base frequecies, for example: { 'A':10,'T':3, C:4 }
        """
        return collections.Counter( self.get_consensus().values() )

    def get_feature_vector(self):
        """ Obtain a feature vector representation of the molecule

        Returns:
            feature_vector(np.array)
        """

        return np.array([
            len(self),
            self.get_strand(),
            self.has_valid_span(),
            self.get_umi_error_rate(),
            self.get_consensus_gc_ratio(),
            len(self.get_raw_barcode_sequences()),
            self.get_safely_aligned_length(),
            self.get_max_mapping_qual(),
            (self.allele is None),
            self.contains_valid_fragment(),
            self.is_multimapped(),
            self.get_undigested_site_count(),
            self.is_valid()
        ])

    def get_alignment_tensor(self,
            max_reads,
            window_radius=20,
            centroid=None,
            mask_centroid=False,
            refence_backed=False,
            skip_missing_reads=False
            ):
        """ Obtain a tensor representation of the molecule alignment around the given centroid
        Args:
            max_reads (int) : maximum amount of reads returned in the tensor, this will be the amount of rows/4 of the returned feature matrix

            window_radius (int) : radius of bp around centroid

            centroid(int) : center of extracted window, when not specified the cut location of the molecule is used

            mask_centroid(bool) : when True, mask reference base at centroid with N

            refence_backed(bool) : when True the molecules reference is used to emit reference bases instead of the MD tag
        Returns:
            tensor_repr(np.array) : (4*window_radius*2*max_reads) dimensional feature matrix
        """
        reference = None
        if refence_backed:
            reference  = self.reference
            if self.reference is None:
                raise ValueError("refence_backed set to True, but the molecule has no reference assigned. Assing one using pysam.FastaFile()")

        height = max_reads
        chromosome = self.chromosome
        if centroid is None:
            _,centroid,strand = self.get_cut_site()
        span_start = centroid-window_radius
        span_end = centroid+window_radius
        span_len = abs(span_start-span_end)
        base_content_table = np.zeros( (height,span_len))
        base_mismatches_table= np.zeros( (height,span_len))
        base_indel_table =np.zeros( (height,span_len))
        base_qual_table =np.zeros( (height,span_len))
        base_clip_table =np.zeros( (height,span_len))
        pointer = 0

        mask = None
        if mask_centroid:
            mask = set((chromosome,centroid))

        for _,frags in self.get_rt_reactions().items() :
            for frag in frags:
                pointer = frag.write_tensor(chromosome, span_start, span_end,
                                            index_start=pointer,
                                           base_content_table=base_content_table,
                                            base_mismatches_table=base_mismatches_table,
                                            base_indel_table=base_indel_table,
                                            base_qual_table=base_qual_table,
                                            base_clip_table=base_clip_table,
                                            height=height,
                                            mask_reference_bases=mask,
                                            reference= reference,
                                            skip_missing_reads=skip_missing_reads
                                           )
        x = np.vstack(
            [
                base_content_table,
                base_mismatches_table,
                base_indel_table,
                base_qual_table,
                base_clip_table
            ])

        return x

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
        for i,(umi,obs) in enumerate( self.umi_counter.most_common() ):
            if i==0:
                mc = obs
            else:
                other += obs

        return mc/(other+mc)


    def get_barcode_sequences(self):
        """Obtain (Cell) barcode sequences associated to molecule

        Returns:
            barcode sequences (set) : barcode sequence(s)
        """
        return set(read.get_tag('BC') for read in self.iter_reads())

    def get_raw_barcode_sequences(self):
        """Obtain (Cell) barcode sequences associated to molecule, not hamming corrected

        Returns:
            barcode sequences (set) : barcode sequence(s)
        """
        return set(read.get_tag('bc') for read in self.iter_reads())


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

        max_show = 6 # maximum amount of fragments to show
        frag_repr = '\n\t'.join([
            textwrap.indent(str(fragment),' '*4) for fragment in self.fragments[:max_show]]
        )

        return f"""{self.__class__.__name__}
        with {len(self.fragments)} assinged fragments
        { "Allele :" +  (self.allele if self.allele is not None else "No allele assigned")}
        """ + frag_repr + ('' if len(self.fragments)<max_show else f'... {len(self.fragments)-max_show} fragments not shown')


    def update_umi(self):
        """Set UMI
        Sets:
            self.umi (str):
                Returns the most common umi associated to the molecule
        """
        self.umi = self.umi_counter.most_common(1)[0][0]

    def get_umi(self):
        """Obtain umi of molecule

        Returns:
            umi (str):
                return main umi associated with this molecule
        """

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

    def get_max_mapping_qual(self):
        """Get max mapping quality of the molecule
        Returns:
            max_mapping_qual (float)
        """
        return max( [fragment.mapping_quality for fragment in self] )


    def contains_valid_fragment(self):
        """Check if an associated fragment exists which returns True for is_valid()

        Returns:
            contains_valid_fragment (bool) : True when any associated fragment is_valid()
        """
        return any(
            (hasattr(fragment,'is_valid') and fragment.is_valid()
            for fragment in self.fragments ) )


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

        # if we already had a fragment, this fragment is a duplicate:
        if len(self.fragments)>1:
            fragment.set_duplicate(True)
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

    def get_safely_aligned_length(self):
        """Get the amount of safely aligned bases (excludes primers)
        Returns:
            aligned_bases (int) : Amount of safely aligned bases
             or None when this cannot be determined because both mates are not mapped
        """
        if self.spanStart is None or self.spanEnd is None:
            return None
        return abs( self.spanEnd - self.spanStart )


    def add_fragment(self, fragment, use_hash=True):
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

        if use_hash:
            if self==fragment:
                self._add_fragment(fragment)
                return True
        else:
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
            can_be_yielded (bool) : True when the molecule is far enough away from the supplied location to be ejected from a buffer.
        """

        if chromosome is None:
            return False
        if chromosome!=self.chromosome:
            return True
        return position < (self.spanStart-self.cache_size*0.5) or position > (self.spanEnd+self.cache_size*0.5)

    def get_rt_reactions(self):
        """Obtain RT reaction dictionary

        returns:
            rt_dict (dict):  {(primer,pos) : [fragment, fragment..] }
        """
        return molecule_to_random_primer_dict(self)

    def get_rt_reaction_fragment_sizes(self, max_N_distance=1):

        """Obtain all RT reaction fragment sizes
        Returns:
            rt_sizes (list of ints)
        """

        rt_reactions = molecule_to_random_primer_dict(self,max_N_distance=max_N_distance)
        amount_of_rt_reactions = len(rt_reactions)

        #this obtains the maximum fragment size:
        frag_chrom, frag_start, frag_end = pysamiterators.iterators.getListSpanningCoordinates([v for v in itertools.chain.from_iterable(self) if v is not None])

        #Obtain the fragment sizes of all RT reactions:
        rt_sizes = []
        for (rt_end,hexamer), fragments in rt_reactions.items():

            if rt_end is None:
                continue

            rt_chrom, rt_start, rt_end = pysamiterators.iterators.getListSpanningCoordinates(
                 itertools.chain.from_iterable([fragment for fragment in fragments if fragment is not None and fragment.get_random_primer_hash()[0] is not None] )
                )

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

    def write_pysam(self, target_file):
        """Write all associated reads to the target file

        Args:
            target_file (pysam.AlignmentFile) : Target file
        """
        for fragment in self:
            fragment.write_pysam(target_file)

    def set_methylation_call_tags(self,
                              call_dict, bismark_call_tag='XM',
                              total_methylated_tag='MC',
                              total_unmethylated_tag='uC',
                              total_methylated_CPG_tag='sZ',
                              total_unmethylated_CPG_tag='sz',
                              total_methylated_CHH_tag='sH',
                              total_unmethylated_CHH_tag='sh',
                              total_methylated_CHG_tag='sX',
                              total_unmethylated_CHG_tag='sx',
                              reads = None
                             ):
        """Set methylation call tags given a methylation dictionary

        This method sets multiple tags in every read associated to the molecule.
        The tags being set are the bismark_call_tag, every aligned base is annotated
        with a zZxXhH or ".", and a tag for both the total methylated C's and unmethylated C's

        Args:
            call_dict (dict) : Dictionary containing bismark calls (chrom,pos) :
                        {'context':letter,'reference_base': letter   , 'consensus': letter }

            bismark_call_tag (str) : tag to write bismark call string

            total_methylated_tag (str) : tag to write total methylated bases

            total_unmethylated_tag (str) : tag to write total unmethylated bases

            reads (iterable) : reads to write the tags to, when not supplied, the tags are written to all associated reads
        Returns:
            can_be_yielded (bool)
        """
        self.methylation_call_dict = call_dict

        # molecule_XM dictionary containing count of contexts
        molecule_XM = collections.Counter(
            list(
                d.get('context','.') for d in self.methylation_call_dict.values() )
                )
        # Contruct XM strings
        if reads is None:
            reads = self.iter_reads()
        for read in reads:
            read.set_tag(
                # Write the methylation tag to the read
                bismark_call_tag,
                ''.join([
                    call_dict.get(
                        (read.reference_name,rpos),{} ).get('context','.')  # Obtain all aligned positions from the call dict
                    for qpos, rpos in read.get_aligned_pairs(matches_only=True) # iterate all positions in the alignment
                    if qpos is not None and rpos is not None]) # make sure to ignore non matching positions ? is this neccesary?
            )



            # Set total methylated bases
            read.set_tag(
                total_methylated_tag,
                molecule_XM['Z']+molecule_XM['X']+molecule_XM['H'] )

            # Set total unmethylated bases
            read.set_tag(
                total_unmethylated_tag,
                molecule_XM['z']+molecule_XM['x']+molecule_XM['h'] )

            # Set total CPG methylated and unmethylated:
            read.set_tag(
                total_methylated_CPG_tag,
                molecule_XM['Z'])

            read.set_tag(
                total_unmethylated_CPG_tag,
                molecule_XM['z'])

            # Set total CHG methylated and unmethylated:
            read.set_tag(
                total_methylated_CHG_tag,
                molecule_XM['X'])

            read.set_tag(
                total_unmethylated_CHG_tag,
                molecule_XM['x'])

            # Set total CHH methylated and unmethylated:
            read.set_tag(
                total_methylated_CHH_tag,
                molecule_XM['H'])

            read.set_tag(
                total_unmethylated_CHH_tag,
                molecule_XM['h'])

    def set_meta(self,tag,value):
        """Set meta information to all fragments

        Args:
            tag (str):
                2 letter tag
            value: value to set

        """
        for f in self:
            f.set_meta(tag,value)

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
        qualities = []
        for fragment in self:
            if fragment.span[0]!=chromosome:
                continue
            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            try:
                start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                R1PrimerLength=fragment.R1_primer_length,
                R2PrimerLength=fragment.R2_primer_length,
                allow_unsafe=(R1 is None))
            except ValueError as e:
                continue

            for read in (R1,R2):
                if read is None:
                    continue
                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                    read,with_seq=False):
                    if query_pos is None or ref_pos != position or read.seq[query_pos]!=base:
                        continue

                    qualities.append( ord( read.qual[query_pos] ) )
        if len(qualities)==0:
            raise IndexError("There are no observations if the supplied base/location combination")
        return np.mean(qualities)

    def get_allele(self, allele_resolver=None, return_allele_informative_base_dict=False):
        """Obtain the allele(s) this molecule maps to

        Args:
            allele_resolver(singlecellmultiomics.alleleTools.AlleleResolver)  : resolver used
            return_allele_informative_base_dict(bool) : return dictionary containing the bases used for allele determination
            defaultdict(list,
            {'allele1': [('chr18', 410937, 'T'),
              ('chr18', 410943, 'G'),
              ('chr18', 410996, 'G'),
              ('chr18', 411068, 'A')]})

        Returns:
            alleles(set( str )) : Set of strings containing associated alleles
        """

        if allele_resolver is None:
            if self.allele_resolver is not None:
                allele_resolver = self.allele_resolver
            else:
                raise ValueError("Supply allele resolver or set it to molecule.allele_resolver")

        alleles = set()
        if return_allele_informative_base_dict:
            aibd = collections.defaultdict(list)
        try:
            for (chrom,pos),base in self.get_consensus(base_obs = self.get_base_observation_dict_NOREF()).items():
                c = allele_resolver.getAllelesAt(chrom,pos,base)
                if c is not None and len(c)==1:
                    alleles.update(c)
                    if return_allele_informative_base_dict:
                        aibd[list(c)[0]].append((chrom,pos,base))

        except Exception as e:
            if return_allele_informative_base_dict:
                return dict()
            else:
                return {}

        if return_allele_informative_base_dict:
            return aibd
        return alleles

    def write_allele_phasing_information_tag(self,allele_resolver=None,tag='ap', reads=None):
        """
        Write allele phasing information to ap tag

        For every associated read a tag wil be written containing:
        chromosome,postion,base,allele_name|chromosome,postion,base,allele_name|...
        for all variants found by the AlleleResolver
        """
        if reads is None:
            reads = self.iter_reads()

        haplotype = self.get_allele(
                return_allele_informative_base_dict=True,
                allele_resolver=allele_resolver)

        phased_locations = [
                    (allele,chromosome, position, base)
                    for allele, bps in haplotype.items()
                    for chromosome, position, base in bps   ]

        phase_str = '|'.join( [f'{chromosome},{position},{base},{allele}' for allele,chromosome, position, base in phased_locations] )

        if len(phase_str)>0:
            for read in reads:
                read.set_tag(tag,phase_str)


    def get_base_observation_dict_NOREF(self):
        '''
        identical to get_base_observation_dict but does not obtain reference bases,
        has to be used when no MD tag is present
        Args:
            return_refbases ( bool ):
                return both observed bases and reference bases

        Returns:
            { genome_location (tuple) : base (string) : obs (int) }
            and
            { genome_location (tuple) : base (string) if return_refbases is True }
        '''

        base_obs = collections.defaultdict(collections.Counter)


        used = 0 #  some alignments yielded valid calls
        ignored = 0
        for fragment in self:
            R1 = fragment.get_R1()
            R2 = fragment.get_R2()
            try:
                start, end = pysamiterators.iterators.getPairGenomicLocations(
                R1=R1,
                R2=R2,
                allow_unsafe=(R1 is None or fragment.unsafe_trimmed))
            except ValueError as e:
                ignored+=1
                continue
            used+=1
            for read in [R1,R2]:
                if read is None:
                    continue
                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                    read,with_seq=False):

                    if query_pos is None or ref_pos is None or ref_pos<start or ref_pos>end:
                        continue
                    query_base = read.seq[query_pos]
                    if query_base=='N':
                        continue
                    base_obs[(read.reference_name,ref_pos)][query_base]+=1

        if used==0 and ignored>0:
            raise ValueError('Could not extract any safe data from molecule')

        return base_obs

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
                R1PrimerLength=fragment.R1_primer_length,
                R2PrimerLength=fragment.R2_primer_length,
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
                    #query_qual = read.qual[query_pos]
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

    def get_consensus(self, base_obs=None, classifier=None, store_consensus=True, reuse_cached_consensus=True):
        """Get dictionary containing consensus calls in respect to reference.
        By default mayority voting is used to determine the consensus base. If a classifier is supplied the classifier is used to determine the consensus base.

        Args:
            base_obs (collections.defaultdict(collections.Counter)) :
                { genome_location (tuple) : base (string) : obs (int) }

            classifier : fitted classifier to use for consensus calling. When no classifier is provided the consensus is determined by majority voting
            store_consensus (bool) : Store the generated consensus for re-use

        Returns:
            consensus (dict)  :  {location : base}
        """
        consensus = {} # postion -> base , key is not set when not decided

        if classifier is not None:

            if reuse_cached_consensus and hasattr( self, 'classifier_consensus' ) and self.classifier_consensus is not None:
                return self.classifier_consensus

            features,reference_bases,CIGAR,alignment_start, alignment_end = self.get_base_calling_feature_matrix_spaced(True)

            predicted_sequence =  classifier.predict(features)
            reference_sequence = ''.join([base for chrom, pos, base  in reference_bases])
            predicted_sequence[ features[:, [ x*8 for x in range(4) ] ].sum(1)==0 ] ='N'

            phred_scores = np.rint(
                    -10*np.log10( np.clip(1-classifier.predict_proba(features).max(1),
                                          0.000000001,
                                          0.999999999 )
                )).astype('B')


            consensus = { (chrom,pos):consensus_base for (chrom, pos, ref_base),consensus_base  in zip( reference_bases,predicted_sequence)  }

            if store_consensus:
                self.classifier_consensus = consensus
                self.classifier_phred_scores = phred_scores
            return consensus


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

        if store_consensus:
            self.majority_consensus = consensus

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
                R1PrimerLength=fragment.R1_primer_length,
                R2PrimerLength=fragment.R2_primer_length,

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

    def get_aligned_reference_bases_dict(self):
        """Get dictionary containing all reference bases to which this molecule aligns
        Returns:
            aligned_reference_positions (dict) :  { (chrom,pos) : 'A', (chrom,pos):'T', .. }
        """
        aligned_reference_positions={}
        for read in self.iter_reads():
            for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True, matches_only=True):
                aligned_reference_positions[(read.reference_name,ref_pos)] = ref_base.upper()
        return aligned_reference_positions


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


    def get_html(self,reference=None,consensus=None, show_reference_sequence=True, show_consensus_sequence=True, reference_bases=None):
        """Get html representation of the molecule
        Returns:
            html_rep(str) : Html representation of the molecule
        """

        html = f"""<h3>{self.chromosome}:{self.spanStart}-{self.spanEnd}
            sample:{self.get_sample()}  {'valid molecule' if self[0].is_valid() else 'Non valid molecule'}</h3>
            <h5>UMI:{self.get_umi()} Mapping qual:{round(self.get_mean_mapping_qual(),1)} Cut loc: {"%s:%s" % self[0].get_site_location()} </h5>
            <div style="white-space:nowrap; font-family:monospace; color:#888">"""
        #undigested:{self.get_undigested_site_count()}
        consensus = self.get_consensus()

        # Obtain reference bases dictionary:
        if reference_bases is None:
            if reference is None:
                reference_bases= self.get_aligned_reference_bases_dict()

            else:
                # obtain reference_bases from reference file
                raise NotImplementedError()

        for fragment in itertools.chain(*self.get_rt_reactions().values()):
            html+= f'<h5>{fragment.get_R1().query_name}</h5>'
            for readid,read in [
                    (1,fragment.get_R1()),
                    (2, fragment.get_R2())]: # go over R1 and R2:
                # This is just the sequence:
                if read is None:
                    continue
                html+= fragment.get_html(
                    self.chromosome,
                    self.spanStart,
                    self.spanEnd,
                    show_read1=(readid==1),
                    show_read2=(readid==2)
                    ) + '<br />'

        # Obtain reference sequence and consensus sequence:
        if consensus is None:
            consensus = self.get_consensus()

        span_len = self.spanEnd - self.spanStart
        visualized = ['.']  * span_len
        reference_vis = ['?']  * span_len
        for location,query_base in consensus.items():
            try:
                if reference_bases is None or reference_bases.get(location,'?')==query_base:
                    visualized[location[1]-self.spanStart] = query_base
                    if reference_bases is not None:
                        reference_vis[location[1]-self.spanStart] = query_base # or reference_bases.get(location,'?')
                else:
                    visualized[location[1]-self.spanStart] = style_str(query_base,color='red',weight=800)
                    if reference_bases is not None:
                        reference_vis[location[1]-self.spanStart] = style_str(reference_bases.get(location,'?'),color='black',weight=800)
            except IndexError as e:
                pass # Tried to visualize a base outside view

        if show_consensus_sequence:
            html+=''.join(visualized) + '<br />'

        if show_reference_sequence:
            html+=''.join(reference_vis) + '<br />'


        html+="</div>"
        return html


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
