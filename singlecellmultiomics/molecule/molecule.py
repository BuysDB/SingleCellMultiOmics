from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
import singlecellmultiomics.bamProcessing
from singlecellmultiomics.fragment import Fragment
from array import array
import itertools
import numpy as np
from singlecellmultiomics.utils import style_str, prob_to_phred, phredscores_to_base_call, base_probabilities_to_likelihood, likelihood_to_prob
import textwrap
import singlecellmultiomics.alleleTools
import functools
import typing
import pysam
import pysamiterators
from singlecellmultiomics.utils import find_ranges, create_MD_tag
import pandas as pd
from uuid import uuid4
from cached_property import cached_property
from collections import Counter, defaultdict


###############

# Variant validation function
def detect_alleles(molecules,
                   contig,
                   position,
                   min_cell_obs=3,
                   base_confidence_threshold=None,
                   classifier=None):  # [ alleles]
    """
    Detect the alleles (variable bases) present at the selected location

    Args:

        molecules : generator to extract molecules from

        variant_location(tuple) : (contig, position) zero based location of the location to test

        min_cell_obs (int) : minimum amount of cells containing the allele to be emitted

        confidence_threshold(float) : minimum confidence of concensus base-call to be taken into account

        classifier (obj) : classifier used for consensus call, when no classifier is supplied a mayority vote is used

    """
    observed_alleles = defaultdict(set)  # cell -> { base_call , .. }
    for molecule in molecules:
        base_call = molecule.get_consensus_base(contig, position, classifier=classifier)

        # confidence = molecule.get_mean_base_quality(*variant_location, base_call)
        if base_call is not None:
            observed_alleles[base_call].add(molecule.sample)

    return [allele for allele, cells in observed_alleles.items() if len(cells) >= min_cell_obs]


def get_variant_phase(molecules, contig, position, variant_base, allele_resolver,
                      phasing_ratio_threshold=None):  # (location,base) -> [( location, base, idenfifier)]
    alleles = [variant_base]
    phases = defaultdict(Counter)  # Allele_id -> variant->obs
    for molecule in molecules:
        # allele_obs = molecule.get_allele(return_allele_informative_base_dict=True,allele_resolver=allele_resolver)
        allele = list(molecule.get_allele(allele_resolver))
        if allele is None or len(allele) > 1 or len(allele) == 0:
            continue
        allele = allele[0]

        base = molecule.get_consensus_base(contig, position)
        if base in alleles:
            phases[base][allele] += 1
        else:
            pass

    if len(phases[variant_base]) == 0:
        raise ValueError("Phasing not established, no gSNVS available")

    phased_allele_id = phases[variant_base].most_common(1)[0][0]

    # Check if the phasing noise is below the threshold:
    if phasing_ratio_threshold is not None:
        correct = phases[variant_base].most_common(1)[0][1]
        total = sum(phases[variant_base].values())
        phasing_ratio = correct / total
        if correct / total < phasing_ratio_threshold:
            raise ValueError(f'Phasing ratio not met. ({phasing_ratio}) < {phasing_ratio_threshold}')
    # Check if the other allele i

    return phased_allele_id


###############

@functools.lru_cache(maxsize=1000)
def might_be_variant(chrom, pos, variants, ref_base=None):
    """Returns True if a variant exists at the given coordinate"""
    if ref_base == 'N':
        return False
    try:
        for record in variants.fetch(chrom, pos, pos + 1):
            return True
    except ValueError as e:
        return False # Happens when the contig does not exists, return False
    return False

def consensii_default_vector():
    """a numpy vector with 5 elements initialsed as zeros"""
    return np.zeros(5)

# 1: read1 2: read2
def molecule_to_random_primer_dict(
        molecule,
        primer_length=6,
        primer_read=2,
        max_N_distance=0):
    rp = defaultdict(list)

    # First add all reactions without a N in the sequence:
    for fragment in molecule:

        h_contig, hstart, hseq = fragment.get_random_primer_hash()
        if hseq is None:
            # This should really not happen with freshly demultiplexed data, it means we cannot extract the random primer sequence
            # which should be present as a tag (rS) in the record
            rp[None, None, None].append(fragment)
        elif 'N' not in hseq:
            rp[h_contig, hstart, hseq].append(fragment)

    # Try to match reactions with N with reactions without a N
    for fragment in molecule:

        h_contig, hstart, hseq = fragment.get_random_primer_hash()
        if hseq is not None and 'N' in hseq:
            # find nearest
            for other_contig, other_start, other_seq in rp:
                if other_contig != h_contig or other_start != hstart:
                    continue

                if hseq.count('N') > max_N_distance:
                    continue

                if 'N' in other_seq:
                    continue

                if hamming_distance(hseq, other_seq) == 0:
                    rp[other_contig, other_start, other_seq].append(fragment)
    return rp


class Molecule():
    """Molecule class, contains one or more associated fragments

    Attributes:
        fragments (list): associated fragments

        spanStart (int): starting coordinate of molecule, None if not available

        spanEnd (int): ending coordinate of molecule, None if not available

        chromosome (str): mapping chromosome of molecule, None if not available

        cache_size (int): radius of molecule assignment cache

        reference (pysam.FastaFile) : reference file, used to obtain base contexts
            and correct aligned_pairs iteration when the MD tag is not correct

        strand (bool): mapping strand.
            True when strand is REVERSE
            False when strand is FORWARD
            None when strand is not determined
    """
    def get_empty_clone(self, fragments=None):
        return type(self)(fragments,
                          cache_size = self.cache_size,
                          reference=self.reference,
                          min_max_mapping_quality=self.min_max_mapping_quality,
                          allele_assingment_method=self.allele_assingment_method,
                          allele_resolver=self.allele_resolver,
                          mapability_reader=self.mapability_reader,
                          max_associated_fragments=self.max_associated_fragments,
                          **self.kwargs)

    def __init__(self,
                 fragments: typing.Optional[typing.Iterable] = None,
                 cache_size: int = 10_000,
                 reference: typing.Union[pysam.FastaFile, pysamiterators.CachedFasta] = None,
                 # When all fragments have a mapping quality below this value
                 # the is_valid method will return False,
                 min_max_mapping_quality: typing.Optional[int] = None,
                 mapability_reader: typing.Optional[singlecellmultiomics.bamProcessing.MapabilityReader] = None,
                 allele_resolver: typing.Optional[singlecellmultiomics.alleleTools.AlleleResolver] = None,
                 max_associated_fragments=None,
                 allele_assingment_method=1, # 0: all variants from the same allele, 1: likelihood
                 **kwargs

                 ):
        """Initialise Molecule

        Parameters
        ----------
        fragments :  list(singlecellmultiomics.fragment.Fragment)
            Fragment to assign to Molecule. More fragments can be added later

        min_max_mapping_quality :  When all fragments have a mapping quality below this value the is_valid method will return False

        allele_resolver :  alleleTools.AlleleResolver or None. Supply an allele resolver in order to assign an allele to the molecule

        mapability_reader : singlecellmultiomics.bamProcessing.MapabilityReader, supply a mapability_reader to set mapping_quality of 0 to molecules mapping to locations which are not mapping uniquely during in-silico library generation.

        cache_size (int): radius of molecule assignment cache

        max_associated_fragments(int) : Maximum amount of fragments associated to molecule. If more fragments are added using add_fragment() they are not added anymore to the molecule

        """
        self.kwargs = kwargs
        self.reference = reference
        self.fragments = []
        self.spanStart = None
        self.spanEnd = None
        self.chromosome = None
        self.cache_size = cache_size
        self.strand = None
        self.umi = None
        self.overflow_fragments = 0
        self.umi_hamming_distance = None
        # when set, when comparing to a fragment the fragment to be added has
        # to match this hash
        self.fragment_match = None
        self.min_max_mapping_quality = min_max_mapping_quality
        self.umi_counter = Counter()  # Observations of umis
        self.max_associated_fragments = max_associated_fragments
        if fragments is not None:
            if isinstance(fragments, list):
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)

        self.allele_resolver = allele_resolver
        self.mapability_reader = mapability_reader
        self.allele_assingment_method = allele_assingment_method
        self.methylation_call_dict = None
        self.finalised = False
        self.obtained_allele_likelihoods = None

    @cached_property
    def can_be_split_into_allele_molecules(self):
        l = self.allele_likelihoods
        if l is None or len(l)<=1:
            return False
        return True

    def split_into_allele_molecules(self):
        """
        Split this molecule into multiple molecules, associated to multiple alleles
        Returns:
            list_of_molecules: list
        """
        # Perform allele based clustering
        allele_clustered_frags = {}
        for fragment in self:
            n = self.get_empty_clone(fragment)
            if n.allele not in allele_clustered_frags:
                allele_clustered_frags[n.allele] = []
            allele_clustered_frags[n.allele].append(n)

        allele_clustered = {}
        for allele, assigned_molecules in allele_clustered_frags.items():
            for i, m in enumerate(assigned_molecules):
                if i == 0:
                    allele_clustered[allele] = m
                else:
                    allele_clustered[allele].add_molecule(m)

        if len(allele_clustered)>1:
            for m in allele_clustered.values():
                m.set_meta('cr', 'SplitUponAlleleClustering')

        return list(allele_clustered.values())

    @cached_property
    def allele(self):
        if self.allele_resolver is None:
            return None
        if self.allele_assingment_method == 0:
            # Obtain allele if available
            if self.allele_resolver is not None:
                try:
                    hits = self.get_allele(self.allele_resolver)
                    # Only store when we have a unique single hit:
                    if len(hits) == 1:
                        self.allele = list(hits)[0]
                except ValueError as e:
                    # This happens when a consensus can not be obtained
                    pass
        elif self.allele_assingment_method == 1:
            al = Counter(self.allele_likelihoods)
            if al is None or len(al)<1:
                return None
            return al.most_common(1)[0][0]

        raise NotImplementedError(f'allele_assingment_method {self.allele_assingment_method} is not defined')

    def __finalise__(self):
        """This function is called when all associated fragments have been gathered"""

        # Perfom allele assignment based on likelihood:
        # this is now only generated upon demand, see .allele method

        if self.mapability_reader is not None:
            self.update_mapability()

        self.finalised = True

    def update_mapability(self, set_mq_zero=False):
        """ Update mapability of this molecule.
        mapping qualities are set to 0 if the mapability_reader returns False
        for site_is_mapable

        The mapability_reader can be set when initiating the molecule, or added later.

        Args:
            set_mq_zero(bool) : set mapping quality of associated reads to 0 when the
            mappability reader returns a bad verdict

        Tip:
            Use `createMapabilityIndex.py` to create an index to feed to the mapability_reader
        """

        mapable = None
        try:
            mapable = self.mapability_reader.site_is_mapable(
                *self.get_cut_site())
        except TypeError:
            pass
        except Exception as e:
            raise

        if mapable is False:
            self.set_meta('mp', 'bad')
            if set_mq_zero:
                for read in self.iter_reads():
                    read.mapping_quality = 0
        elif mapable is True:
            self.set_meta('mp', 'unique')
        else:
            self.set_meta('mp', 'unknown')

    def calculate_consensus(self, consensus_model, molecular_identifier, out, **model_kwargs):
        """
        Create consensus read for molecule

        Args:

            consensus_model

            molecular_identifier (str) : identier for this molecule, will be suffixed to the reference_id

            out(pysam.AlingmentFile) : target bam file

            **model_kwargs : arguments passed to the consensus model

        """
        try:
            consensus_reads = self.deduplicate_to_single_CIGAR_spaced(
                out,
                f'c_{self.get_a_reference_id()}_{molecular_identifier}',
                consensus_model,
                NUC_RADIUS=model_kwargs['consensus_k_rad']
            )
            for consensus_read in consensus_reads:
                consensus_read.set_tag('RG', self[0].get_read_group())
                consensus_read.set_tag('mi', molecular_identifier)
                out.write(consensus_read)

        except Exception as e:

            self.set_rejection_reason('CONSENSUS_FAILED', set_qcfail=True)
            self.write_pysam(out)

    def get_a_reference_id(self):
        """
        Obtain a reference id for a random associated mapped read
        """
        for read in self.iter_reads():
            if not read.is_unmapped:
                return read.reference_id
        return -1

    def get_consensus_read(self, target_file,
                           read_name, consensus=None,
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
            try:  # Obtain the consensus sequence
                consensus = self.get_consensus()
            except Exception as e:
                raise
        if isinstance(consensus, str):
            sequence = consensus
        else:
            sequence = ''.join(
                (consensus.get(
                    (self.chromosome, ref_pos), 'N') for ref_pos in range(
                    self.spanStart, self.spanEnd + 1)))

        if isinstance(phred_scores, dict):
            phred_score_array = list(
                phred_scores.get(
                    (self.chromosome, ref_pos), 0) for ref_pos in range(
                    self.spanStart, self.spanEnd + 1))
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
            cread.cigarstring = cigarstring
        else:
            cread.cigarstring = f'{len(sequence)}M'
        cread.mapping_quality = self.get_max_mapping_qual()

        cread.is_reverse = self.strand
        if mdstring is not None:
            cread.set_tag('MD', mdstring)

        self.write_tags_to_psuedoreads((cread,))

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
            len(self.get_raw_barcode_sequences()),
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
            tag_obs (defaultdict(Counter)):
                { tag(str) : { value(int/str): frequency:(int) }


        """
        tags_obs = defaultdict(Counter)
        for tag, value in itertools.chain(
                *[r.tags for r in self.iter_reads()]):
            try:
                tags_obs[tag][value] += 1
            except TypeError:
                # Dont count arrays for example
                pass
        return tags_obs

    def write_tags(self):
        """ Write BAM tags to all reads associated to this molecule

        This function sets the following tags:
            - mI : most common umi
            - DA : allele
            - af : amount of associated fragments
            - rt : rt_reaction_index
            - rd : rt_duplicate_index
            - TR : Total RT reactions
            - ap : phasing information (if allele_resolver is set)
            - TF : total fragments
            - ms : size of the molecule (largest fragment)
        """
        self.is_valid(set_rejection_reasons=True)
        if self.umi is not None:
            self.set_meta('mI', self.umi)
        if self.allele is not None:
            self.set_meta('DA', str(self.allele))

        # Set total amount of associated fragments
        self.set_meta('TF', len(self.fragments) + self.overflow_fragments)
        try:
            self.set_meta('ms',self.estimated_max_length)
        except Exception as e:
            # There is no properly defined aligned length
            pass
        # associatedFragmentCount :
        self.set_meta('af', len(self))
        for rc, frag in enumerate(self):
            frag.set_meta('RC', rc)
            if rc > 0:
                # Set duplicate bit
                for read in frag:
                    if read is not None:
                        read.is_duplicate = True

        # Write RT reaction tags (rt: rt reaction index, rd rt duplicate index)
        # This is only required for fragments which have defined random primers
        rt_reaction_index = None
        for rt_reaction_index, ( (contig, random_primer_start, random_primer_sequence), frags) in enumerate(
                self.get_rt_reactions().items()):

            for rt_duplicate_index, frag in enumerate(frags):
                frag.set_meta('rt', rt_reaction_index)
                frag.set_meta('rd', rt_duplicate_index)
                frag.set_meta('rp', random_primer_start)
        self.set_meta('TR', 0 if (rt_reaction_index is None) else rt_reaction_index + 1)

        if self.allele_resolver is not None:
            self.write_allele_phasing_information_tag()

    def write_tags_to_psuedoreads(self, reads):
        """
        Write molecule information to the supplied reads as BAM tags
        """
        # write methylation tags to new reads if applicable:
        if self.methylation_call_dict is not None:
            self.set_methylation_call_tags(
                self.methylation_call_dict, reads=reads)

        for read in reads:
            read.set_tag('SM', self.sample)
            if hasattr(self, 'get_cut_site'):
                read.set_tag('DS', self.get_cut_site()[1])

            if self.umi is not None:
                read.set_tag('RX', self.umi)
                bc = list(self.get_barcode_sequences())[0]
                read.set_tag('BC', bc)
                read.set_tag('MI', bc + self.umi)

            # Store total amount of RT reactions:
            read.set_tag('TR', len(self.get_rt_reactions()))
            read.set_tag('TF', len(self.fragments) + self.overflow_fragments)

            if self.allele is not None:
                read.set_tag('DA', self.allele)

        if self.allele_resolver is not None:
            self.write_allele_phasing_information_tag(
                self.allele_resolver, reads=reads)

    def deduplicate_to_single(
            self,
            target_bam,
            read_name,
            classifier,
            reference=None):
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

        # We only use the proba:
        base_calling_probs = classifier.predict_proba(features)
        predicted_sequence = ['ACGT'[i] for i in np.argmax(base_calling_probs, 1)]
        phred_scores = np.rint(-10 * np.log10(np.clip(1 - base_calling_probs.max(1), 0.000000001, 0.999999))).astype(
            'B')

        read = self.get_consensus_read(
            read_name=read_name,
            target_file=target_bam,
            consensus=''.join(predicted_sequence),
            phred_scores=phred_scores)
        read.is_read1 = True
        return read

    def deduplicate_to_single_CIGAR_spaced(
            self,
            target_bam,
            read_name,
            classifier=None,
            max_N_span=300,
            reference=None,
            **feature_matrix_args
    ):
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

        if classifier is not None:
            features, reference_bases, CIGAR, alignment_start, alignment_end = self.get_base_calling_feature_matrix_spaced(
                True, reference=reference, **feature_matrix_args)

            base_calling_probs = classifier.predict_proba(features)
            predicted_sequence = ['ACGT'[i] for i in np.argmax(base_calling_probs, 1)]

            reference_sequence = ''.join(
                [base for chrom, pos, base in reference_bases])
            # predicted_sequence[ features[:, [ x*8 for x in range(4) ] ].sum(1)==0 ] ='N'
            predicted_sequence = ''.join(predicted_sequence)

            phred_scores = np.rint(
                -10 * np.log10(np.clip(1 - base_calling_probs.max(1),
                                       0.000000001,
                                       0.999999)
                               )).astype('B')

        reads = []

        query_index_start = 0
        query_index_end = 0
        reference_position = alignment_start  # pointer to current position
        reference_start = alignment_start  # pointer to alignment start of current read
        supplementary = False
        partial_CIGAR = []
        partial_MD = []

        for operation, amount in CIGAR:
            if operation == 'M':  # Consume query and reference
                query_index_end += amount
                reference_position += amount
                partial_CIGAR.append(f'{amount}{operation}')

            if operation == 'N':
                # Consume reference:
                reference_position += amount
                if amount > max_N_span:  # Split up in supplementary alignment
                    # Eject previous
                    # reference_seq =

                    consensus_read = self.get_consensus_read(
                        read_name=read_name,
                        target_file=target_bam,
                        consensus=predicted_sequence[query_index_start:query_index_end],
                        phred_scores=phred_scores[query_index_start:query_index_end],
                        cigarstring=''.join(partial_CIGAR),
                        mdstring=create_MD_tag(
                            reference_sequence[query_index_start:query_index_end],
                            predicted_sequence[query_index_start:query_index_end]
                        ),
                        start=reference_start,
                        supplementary=supplementary
                    )
                    reads.append(consensus_read)
                    if not supplementary:
                        consensus_read.is_read1 = True

                    supplementary = True
                    # Start new:
                    query_index_start = query_index_end
                    reference_start = reference_position
                    partial_CIGAR = []
                else:
                    partial_CIGAR.append(f'{amount}{operation}')

        reads.append(self.get_consensus_read(
            read_name=read_name,
            target_file=target_bam,
            consensus=''.join(predicted_sequence[query_index_start:query_index_end]),
            phred_scores=phred_scores[query_index_start:query_index_end],
            cigarstring=''.join(partial_CIGAR),
            mdstring=create_MD_tag(
                reference_sequence[query_index_start:query_index_end],
                predicted_sequence[query_index_start:query_index_end]

            ),
            start=reference_start,
            supplementary=supplementary
        ))

        # Write last index tag to last read ..
        if supplementary:
            reads[-1].is_read2 = True

        # Write NH tag (the amount of records with the same query read):
        for read in reads:
            read.set_tag('NH', len(reads))

        return reads

    def extract_stretch_from_dict(self, base_call_dict, alignment_start, alignment_end):
        base_calling_probs = np.array(
            [base_call_dict.get((self.chromosome, pos), ('N', 0))[1] for pos in range(alignment_start, alignment_end)])
        predicted_sequence = [base_call_dict.get((self.chromosome, pos), ('N', 0))[0] for pos in
                              range(alignment_start, alignment_end)]
        predicted_sequence = ''.join(predicted_sequence)
        phred_scores = np.rint(
            -10 * np.log10(np.clip(1 - base_calling_probs,
                                   0.000000001,
                                   0.999999999)
                           )).astype('B')
        return predicted_sequence, phred_scores

    def get_base_confidence_dict(self):
        """
        Get dictionary containing base calls per position and the corresponding confidences

        Returns:
            obs (dict) :  (contig (str), position  (int) ) : base (str) : prob correct (list)
        """
        # Convert (contig, position) -> (base_call) into:
        # (contig, position) -> (base_call, confidence)
        obs = defaultdict(lambda: defaultdict(list))
        for read in self.iter_reads():
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                qbase = read.seq[qpos]
                qqual = read.query_qualities[qpos]
                # @ todo reads which span multiple chromosomes
                obs[(self.chromosome, rpos)][qbase].append(1 - np.power(10, -qqual / 10))
        return obs


    def deduplicate_majority(self, target_bam, read_name, max_N_span=None):

        obs = self.get_base_confidence_dict()

        reads = list(self.get_dedup_reads(read_name,
                                     target_bam,
                                     obs={reference_position: phredscores_to_base_call(probs)
                                          for reference_position, probs in obs.items()},
                                     max_N_span=max_N_span))
        self.write_tags_to_psuedoreads([read for read in reads if read is not None])
        return reads

    def generate_partial_reads(self, obs, max_N_span=None):
        CIGAR, alignment_start, alignment_end = self.get_CIGAR()
        query_index_start = 0
        query_index_end = 0
        reference_position = alignment_start  # pointer to current position
        reference_start = alignment_start  # pointer to alignment start of current read
        reference_end = None
        partial_CIGAR = []
        partial_MD = []
        partial_sequence = []
        partial_phred = []

        for operation, amount in CIGAR:
            if operation == 'N':
                if max_N_span is not None and amount > max_N_span:
                    yield reference_start, reference_end, partial_sequence, partial_phred, partial_CIGAR, partial_MD
                    # Clear all
                    partial_CIGAR = []
                    partial_MD = []
                    partial_sequence = []
                    partial_phred = []
                else:
                    # Increment
                    partial_CIGAR.append(f'{amount}{operation}')
                    query_index_start += sum((len(s) for s in partial_sequence))

                reference_position += amount
            elif operation == 'M':  # Consume query and reference

                query_index_end += amount
                if len(partial_CIGAR) == 0:
                    reference_start = reference_position
                start_fetch = reference_position

                reference_position += amount
                reference_end = reference_position
                partial_CIGAR.append(f'{amount}{operation}')

                predicted_sequence, phred_scores = self.extract_stretch_from_dict(obs, start_fetch, reference_end)  # [start .. end)

                partial_sequence.append(predicted_sequence)
                partial_phred.append(phred_scores)

        yield reference_start, reference_end, partial_sequence, partial_phred, partial_CIGAR, partial_MD

    def get_dedup_reads(self, read_name, target_bam, obs, max_N_span=None):
        if self.chromosome is None:
            return None # We cannot perform this action
        for reference_start, reference_end, partial_sequence, partial_phred, partial_CIGAR, partial_MD in self.generate_partial_reads(
                obs, max_N_span=max_N_span):
            consensus_read = self.get_consensus_read(
                read_name=read_name,
                target_file=target_bam,
                consensus=''.join(partial_sequence),
                phred_scores= array('B', np.concatenate(partial_phred)), # Needs to be casted to array
                cigarstring=''.join(partial_CIGAR),
                mdstring=create_MD_tag(
                    self.reference.fetch(self.chromosome, reference_start, reference_end),
                    ''.join(partial_sequence)
                ),
                start=reference_start,
                supplementary=False
            )

            consensus_read.is_reverse = self.strand
            yield consensus_read

    def deduplicate_to_single_CIGAR_spaced_from_dict(
            self,
            target_bam,
            read_name,
            base_call_dict,  # (contig, position) -> (base_call, confidence)
            max_N_span=300,
    ):
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

        CIGAR, alignment_start, alignment_end = self.get_CIGAR()

        reads = []

        query_index_start = 0
        query_index_end = 0
        reference_position = alignment_start  # pointer to current position
        reference_start = alignment_start  # pointer to alignment start of current read
        reference_end = None
        supplementary = False
        partial_CIGAR = []
        partial_MD = []

        partial_sequence = []
        partial_phred = []

        for operation, amount in CIGAR:

            if operation == 'N':
                # Pop the previous read..
                if len(partial_sequence):
                    assert reference_end is not None
                    consensus_read = self.get_consensus_read(
                        read_name=read_name,
                        target_file=target_bam,
                        consensus=''.join(partial_sequence),
                        phred_scores=array('B',np.concatenate(partial_phred)),
                        cigarstring=''.join(partial_CIGAR),
                        mdstring=create_MD_tag(
                            self.reference.fetch(self.chromosome, reference_start, reference_end),
                            ''.join(partial_sequence)
                        ),
                        start=reference_start,
                        supplementary=supplementary
                    )
                    reads.append(consensus_read)
                    if not supplementary:
                        consensus_read.is_read1 = True

                    supplementary = True
                    reference_start = reference_position
                    partial_CIGAR = []
                    partial_phred = []
                    partial_sequence = []

                # Consume reference:
                reference_position += amount
                partial_CIGAR.append(f'{amount}{operation}')

            if operation == 'M':  # Consume query and reference
                query_index_end += amount
                # This should only be reset upon a new read:
                if len(partial_CIGAR) == 0:
                    reference_start = reference_position
                reference_position += amount
                reference_end = reference_position

                partial_CIGAR.append(f'{amount}{operation}')

                predicted_sequence, phred_scores = self.extract_stretch_from_dict(base_call_dict, reference_start,
                                                                                  reference_end)  # [start .. end)

                partial_sequence.append(predicted_sequence)
                partial_phred.append(phred_scores)

        consensus_read = self.get_consensus_read(
            read_name=read_name,
            target_file=target_bam,
            consensus=''.join(partial_sequence),
            phred_scores= array('B',np.concatenate(partial_phred)),
            cigarstring=''.join(partial_CIGAR),
            mdstring=create_MD_tag(
                self.reference.fetch(self.chromosome, reference_start, reference_end),
                ''.join(partial_sequence)
            ),
            start=reference_start,
            supplementary=supplementary
        )
        reads.append(consensus_read)
        if not supplementary:
            consensus_read.is_read1 = True

        supplementary = True
        reference_start = reference_position
        partial_CIGAR = []

        # Write last index tag to last read ..
        if supplementary:
            reads[-1].is_read2 = True
            reads[0].is_read1 = True

        # Write NH tag (the amount of records with the same query read):
        for read in reads:
            read.set_tag('NH', len(reads))

        return reads

    def get_base_calling_feature_matrix(
            self,
            return_ref_info=False,
            start=None,
            end=None,
            reference=None,
            NUC_RADIUS=1,
            USE_RT=True,
            select_read_groups=None):
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
            MATE_INDEX = 3
            CYCLE_INDEX = 4
            MQ_INDEX = 5
            FS_INDEX = 6

            COLUMN_OFFSET = 0
            features_per_block = 8 - (not USE_RT)

            origin_start = start
            origin_end = end

            end += NUC_RADIUS
            start -= NUC_RADIUS

            features = np.zeros(
                (end - start + 1, (features_per_block * BASE_COUNT) + COLUMN_OFFSET))

            if return_ref_info:
                ref_bases = {}

            for rt_id, fragments in self.get_rt_reactions().items():
                # we need to keep track what positions where covered by this RT
                # reaction
                RT_reaction_coverage = set()  # (pos, base_call)
                for fragment in fragments:
                    for read in fragment:
                        if select_read_groups is not None:
                            if not read.has_tag('RG'):
                                raise ValueError(
                                    "Not all reads in the BAM file have a read group defined.")
                            if not read.get_tag('RG') in select_read_groups:
                                continue
                        # Skip reads outside range
                        if read is None or read.reference_start > (
                                end + 1) or read.reference_end < start:
                            continue
                        for cycle, q_pos, ref_pos, ref_base in pysamiterators.ReadCycleIterator(
                                read, matches_only=True, with_seq=True, reference=reference):

                            row_index = ref_pos - start
                            if row_index < 0 or row_index >= features.shape[0]:
                                continue

                            query_base = read.seq[q_pos]
                            # Base index block:
                            block_index = 'ACGTN'.index(query_base)

                            # Update rt_reactions
                            if USE_RT:
                                if not (
                                               ref_pos, query_base) in RT_reaction_coverage:
                                    features[row_index][RT_INDEX +
                                                        COLUMN_OFFSET +
                                                        features_per_block *
                                                        block_index] += 1
                                RT_reaction_coverage.add((ref_pos, query_base))

                            # Update total phred score
                            features[row_index][PHRED_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += read.query_qualities[q_pos]

                            # Update total reads

                            features[row_index][RC_INDEX + COLUMN_OFFSET +
                                                features_per_block * block_index] += 1

                            # Update mate index
                            features[row_index][MATE_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += read.is_read2

                            # Update fragment sizes:
                            features[row_index][FS_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += abs(fragment.span[1] -
                                                                    fragment.span[2])

                            # Update cycle
                            features[row_index][CYCLE_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += cycle

                            # Update MQ:
                            features[row_index][MQ_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += read.mapping_quality

                            # update strand:
                            features[row_index][STRAND_INDEX +
                                                COLUMN_OFFSET +
                                                features_per_block *
                                                block_index] += read.is_reverse

                            if return_ref_info:
                                row_index_in_output = ref_pos - origin_start
                                if row_index_in_output < 0 or row_index_in_output >= origin_end - origin_start + 1:
                                    continue

                                ref_bases[ref_pos] = ref_base.upper()

            # Normalize all and return

            for block_index in range(BASE_COUNT):  # ACGTN
                for index in (
                        PHRED_INDEX,
                        MATE_INDEX,
                        CYCLE_INDEX,
                        MQ_INDEX,
                        FS_INDEX,
                        STRAND_INDEX):
                    features[:, index +
                                COLUMN_OFFSET +
                                features_per_block *
                                block_index] /= features[:, RC_INDEX +
                                                            COLUMN_OFFSET +
                                                            features_per_block *
                                                            block_index]
            # np.nan_to_num( features, nan=-1, copy=False )
            features[np.isnan(features)] = -1

            if NUC_RADIUS > 0:
                # duplicate columns in shifted manner
                x = features
                features = np.zeros(
                    (x.shape[0] - NUC_RADIUS * 2, x.shape[1] * (1 + NUC_RADIUS * 2)))
                for offset in range(0, NUC_RADIUS * 2 + 1):
                    slice_start = offset
                    slice_end = -(NUC_RADIUS * 2) + offset
                    if slice_end == 0:
                        features[:, features_per_block *
                                    BASE_COUNT *
                                    offset:features_per_block *
                                           BASE_COUNT *
                                           (offset +
                                            1)] = x[slice_start:, :]
                    else:
                        features[:, features_per_block *
                                    BASE_COUNT *
                                    offset:features_per_block *
                                           BASE_COUNT *
                                           (offset +
                                            1)] = x[slice_start:slice_end, :]

            if return_ref_info:
                ref_info = [
                    (self.chromosome, ref_pos, ref_bases.get(ref_pos, 'N'))
                    for ref_pos in range(origin_start, origin_end + 1)]
                return features, ref_info
            return features

    def get_CIGAR(self, reference=None):
        """ Get alignment of all associated reads

        Returns:
            y : reference bases
            CIGAR : alignment of feature matrix to reference tuples (operation, count)
            reference(pysam.FastaFile) : reference to fetch reference bases from, if not supplied the MD tag is used
        """

        X = None

        CIGAR = []
        prev_end = None
        alignment_start = None
        alignment_end = None
        for start, end in self.get_aligned_blocks():

            if prev_end is not None:
                CIGAR.append(('N', start - prev_end - 1))
            CIGAR.append(('M', (end - start + 1)))
            prev_end = end

            if alignment_start is None:
                alignment_start = start
                alignment_end = end
            else:
                alignment_start = min(alignment_start, start)
                alignment_end = max(alignment_end, end)

        return CIGAR, alignment_start, alignment_end

    @functools.lru_cache(maxsize=4)
    def get_base_calling_feature_matrix_spaced(
            self,
            return_ref_info=False,
            reference=None,
            **feature_matrix_args):
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
        for start, end in self.get_aligned_blocks():
            if return_ref_info:
                x, y_ = self.get_base_calling_feature_matrix(
                    return_ref_info=return_ref_info, start=start, end=end,
                    reference=reference, **feature_matrix_args
                )
                y += y_
            else:
                x = self.get_base_calling_feature_matrix(
                    return_ref_info=return_ref_info,
                    start=start,
                    end=end,
                    reference=reference,
                    **feature_matrix_args)
            if X is None:
                X = x
            else:
                X = np.append(X, x, axis=0)

            if prev_end is not None:
                CIGAR.append(('N', start - prev_end - 1))
            CIGAR.append(('M', (end - start + 1)))
            prev_end = end

            if alignment_start is None:
                alignment_start = start
                alignment_end = end
            else:
                alignment_start = min(alignment_start, start)
                alignment_end = max(alignment_end, end)

        if return_ref_info:
            return X, y, CIGAR, alignment_start, alignment_end
        else:
            return X, CIGAR, alignment_start, alignment_end

    def get_base_calling_training_data(
            self,
            mask_variants=None,
            might_be_variant_function=None,
            reference=None,
            **feature_matrix_args):
        if mask_variants is not None and might_be_variant_function is None:
            might_be_variant_function = might_be_variant

        features, feature_info, _CIGAR, _alignment_start, _alignment_end = self.get_base_calling_feature_matrix_spaced(
            True, reference=reference, **feature_matrix_args)

        # Edgecase: it can be that not a single base can be used for base calling
        # in that case features will be None
        # when there is no features return None
        if features is None or len(features) == 0:
            return None

        # check which bases should not be used
        use_indices = [
            mask_variants is None or
            not might_be_variant_function(chrom, pos, mask_variants, base)
            for chrom, pos, base in feature_info]

        X_molecule = features[use_indices]
        y_molecule = [
            base for use, (chrom, pos, base) in
            zip(use_indices, feature_info) if use
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

    def set_rejection_reason(self, reason, set_qcfail=False):
        """ Add rejection reason to all fragments associated to this molecule

        Args:
            reason (str) : rejection reason to set

            set_qcfail(bool) : set qcfail bit to True for all associated reads
        """
        for fragment in self:
            fragment.set_rejection_reason(reason, set_qcfail=set_qcfail)

    def is_valid(self, set_rejection_reasons=False):
        """Check if the molecule is valid
        All of the following requirements should be met:
        - no multimapping
        - no low mapping mapping_quality (Change molecule.min_max_mapping_quality to set the threshold)
        - molecule is associated with at least one valid fragment

        Args:
            set_rejection_reasons (bool) : When set to True, all reads get a
            rejection reason (RR tag) written to them if the molecule is rejected.

        Returns:
            is_valid (bool) : True when all requirements are met, False otherwise

        """
        if self.is_multimapped():
            if set_rejection_reasons:
                self.set_rejection_reason('multimapping')
            return False

        if self.min_max_mapping_quality is not None and \
                self.get_max_mapping_qual() < self.min_max_mapping_quality:
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
                 for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False)))))
        )

    def __len__(self):
        """Obtain the amount of fragments associated to the molecule"""
        return len(self.fragments)

    def get_consensus_base_frequencies(self, allow_N=False):
        """Obtain the frequency of bases in the molecule consensus sequence

        Returns:
            base_frequencies (Counter) : Counter containing base frequecies, for example: { 'A':10,'T':3, C:4 }
        """
        return Counter(
            self.get_consensus(
                allow_N=allow_N).values())

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
            reference = self.reference
            if self.reference is None:
                raise ValueError(
                    "refence_backed set to True, but the molecule has no reference assigned. Assing one using pysam.FastaFile()")

        height = max_reads
        chromosome = self.chromosome
        if centroid is None:
            _, centroid, strand = self.get_cut_site()
        span_start = centroid - window_radius
        span_end = centroid + window_radius
        span_len = abs(span_start - span_end)
        base_content_table = np.zeros((height, span_len))
        base_mismatches_table = np.zeros((height, span_len))
        base_indel_table = np.zeros((height, span_len))
        base_qual_table = np.zeros((height, span_len))
        base_clip_table = np.zeros((height, span_len))
        pointer = 0

        mask = None
        if mask_centroid:
            mask = set((chromosome, centroid))

        for _, frags in self.get_rt_reactions().items():
            for frag in frags:
                pointer = frag.write_tensor(
                    chromosome,
                    span_start,
                    span_end,
                    index_start=pointer,
                    base_content_table=base_content_table,
                    base_mismatches_table=base_mismatches_table,
                    base_indel_table=base_indel_table,
                    base_qual_table=base_qual_table,
                    base_clip_table=base_clip_table,
                    height=height,
                    mask_reference_bases=mask,
                    reference=reference,
                    skip_missing_reads=skip_missing_reads)
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
        return (bf['G'] + bf['C']) / sum(bf.values())

    def get_umi_error_rate(self):
        """Obtain fraction of fragments that are associated
        to the molecule with a exact matching UMI vs total amount of associated fragments
        Returns:
            exact_matching_fraction (float)
        """
        mc = 0
        other = 0
        for i, (umi, obs) in enumerate(self.umi_counter.most_common()):
            if i == 0:
                mc = obs
            else:
                other += obs

        return mc / (other + mc)

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

        max_show = 6  # maximum amount of fragments to show
        frag_repr = '\n\t'.join([textwrap.indent(str(fragment), ' ' * 4)
                                 for fragment in self.fragments[:max_show]])

        return f"""{self.__class__.__name__}
        with {len(self.fragments)} assinged fragments
        {"Allele :" + (self.allele if self.allele is not None else "No allele assigned")}
        """ + frag_repr + (
            '' if len(self.fragments) < max_show else f'... {len(self.fragments) - max_show} fragments not shown')

    def update_umi(self):
        """Set UMI
        sets self.umi (str) sets the most common umi associated to the molecule
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
            try:
                site = fragment.get_site_location()
            except AttributeError:
                return None
            if site is not None:
                return tuple((*site, fragment.get_strand()))
        return None

    def get_mean_mapping_qual(self):
        """Get mean mapping quality of the molecule

        Returns:
            mean_mapping_qual (float)
        """
        return np.mean([fragment.mapping_quality for fragment in self])

    def get_max_mapping_qual(self):
        """Get max mapping quality of the molecule
        Returns:
            max_mapping_qual (float)
        """
        return max([fragment.mapping_quality for fragment in self])

    def get_min_mapping_qual(self) -> float:
        """Get min mapping quality of the molecule
        Returns:
            min_mapping_qual (float)
        """
        return min([fragment.mapping_quality for fragment in self])

    def contains_valid_fragment(self):
        """Check if an associated fragment exists which returns True for is_valid()

        Returns:
            contains_valid_fragment (bool) : True when any associated fragment is_valid()
        """
        return any(
            (hasattr(fragment, 'is_valid') and fragment.is_valid()
             for fragment in self.fragments))

    def is_multimapped(self):
        """Check if the molecule is multimapping

        Returns:
            is_multimapped (bool) : True when multimapping
        """
        for fragment in self.fragments:
            if not fragment.is_multimapped:
                return False
        return True

    def add_molecule(self, other):
        """
        Merge other molecule into this molecule.
        Merges by assigning all fragments in other to this molecule.
        """
        for fragment in other:
            self._add_fragment(fragment)


    def get_span_sequence(self, reference=None):
        """Obtain the sequence between the start and end of the molecule
        Args:
            reference(pysam.FastaFile) : reference  to use.
                If not specified `self.reference` is used
        Returns:
            sequence (str)
        """
        if self.chromosome is None:
            return ''

        if reference is None:
            if self.reference is None:
                raise ValueError('Please supply a reference (PySAM.FastaFile)')
            reference = self.reference
        return reference.fetch(
            self.chromosome,
            self.spanStart,
            self.spanEnd).upper()

    def get_fragment_span_sequence(self, reference=None):
        return self.get_span_sequence(reference)

    def _add_fragment(self, fragment):

        # Do not process the fragment when the max_associated_fragments threshold is exceeded
        if self.max_associated_fragments is not None and len(self.fragments) >= (self.max_associated_fragments):
            self.overflow_fragments += 1
            raise OverflowError()

        self.match_hash = fragment.match_hash

        # if we already had a fragment, this fragment is a duplicate:
        if len(self.fragments) > 1:
            fragment.set_duplicate(True)

        self.fragments.append(fragment)

        # Update span:
        add_span = fragment.get_span()

        # It is possible that the span is not defined, then set the respective keys to None
        # This indicates the molecule is qcfail

        #if not self.has_valid_span():
        #    self.spanStart, self.spanEnd, self.chromosome = None,None, None
        #else:
        self.spanStart = add_span[1] if self.spanStart is None else min(
            add_span[1], self.spanStart)
        self.spanEnd = add_span[2] if self.spanEnd is None else max(
            add_span[2], self.spanEnd)
        self.chromosome = add_span[0]

        self.span = (self.chromosome, self.spanStart, self.spanEnd)
        if fragment.strand is not None:
            self.strand = fragment.strand
        self.umi_counter[fragment.umi] += 1
        self.umi_hamming_distance = fragment.umi_hamming_distance
        self.saved_base_obs = None
        self.update_umi()
        return True

    @property
    def aligned_length(self) -> int:
        if self.has_valid_span():
            return self.spanEnd - self.spanStart
        else:
            return None

    @property
    def is_completely_matching(self) -> bool:
        """
        Checks if all associated reads are completely mapped:
        checks if all cigar operations are M,
        Returns True when all cigar operations are M, False otherwise
        """

        return all(
                (
                     all(
                     [ (operation==0)
                        for operation, amount in read.cigartuples] )
                for read in self.iter_reads()))


    @property
    def estimated_max_length(self) -> int:
        """
        Obtain the estimated size of the fragment,
        returns None when estimation is not possible
        Takes into account removed bases (R2)
        Assumes inwards sequencing orientation
        """
        max_size = None
        for frag in self:
            r = frag.estimated_length
            if r is None :
                continue
            if max_size is None:
                max_size = r
            elif r>max_size:
                max_size = r
        return max_size

    def get_safely_aligned_length(self):
        """Get the amount of safely aligned bases (excludes primers)
        Returns:
            aligned_bases (int) : Amount of safely aligned bases
             or None when this cannot be determined because both mates are not mapped
        """
        if self.spanStart is None or self.spanEnd is None:
            return None

        start = None
        end = None
        contig = None
        for fragment in self:
            if not fragment.safe_span:
                continue

            if contig is None:
                contig = fragment.span[0]
            if contig == fragment.span[0]:
                f_start, f_end = fragment.get_safe_span()
                if start is None:
                    start = f_start
                    end = f_end
                else:
                    start = min(f_start, start)
                    end = min(f_end, end)

        if end is None:
            raise ValueError('Not safe')
        return abs(end - start)

    def add_fragment(self, fragment, use_hash=True):
        """Associate a fragment with this Molecule

        Args:
            fragment (singlecellmultiomics.fragment.Fragment) : Fragment to associate
        Returns:
            has_been_added (bool) : Returns False when the fragments which have already been associated to the molecule refuse the fragment

        Raises:
            OverflowError : when too many fragments have been associated already
                            control this with .max_associated_fragments attribute
        """

        if len(self.fragments) == 0:
            self._add_fragment(fragment)
            self.sample = fragment.sample
            return True

        if use_hash:
            if self == fragment:
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
        if chromosome != self.chromosome:
            return True
        return position < (
                self.spanStart -
                self.cache_size *
                0.5) or position > (
                       self.spanEnd +
                       self.cache_size *
                       0.5)

    def get_rt_reactions(self) -> dict:
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

        rt_reactions = molecule_to_random_primer_dict(
            self, max_N_distance=max_N_distance)
        amount_of_rt_reactions = len(rt_reactions)

        # this obtains the maximum fragment size:
        frag_chrom, frag_start, frag_end = pysamiterators.iterators.getListSpanningCoordinates(
            [v for v in itertools.chain.from_iterable(self) if v is not None])

        # Obtain the fragment sizes of all RT reactions:
        rt_sizes = []
        for (rt_contig, rt_end, hexamer), fragments in rt_reactions.items():

            if rt_end is None:
                continue

            rt_chrom, rt_start, rt_end = pysamiterators.iterators.getListSpanningCoordinates(
                itertools.chain.from_iterable(
                    [fragment for fragment in fragments if
                     fragment is not None and fragment.get_random_primer_hash()[0] is not None]))

            rt_sizes.append([rt_end - rt_start])
        return rt_sizes

    def get_mean_rt_fragment_size(self):
        """Obtain the mean RT reaction fragment size

        Returns:
            mean_rt_size (float)
        """
        return np.nanmean(
            self.get_rt_reaction_fragment_sizes()
        )

    def write_pysam(self, target_file, consensus=False, no_source_reads=False, consensus_name=None, consensus_read_callback=None, consensus_read_callback_kwargs=None):
        """Write all associated reads to the target file

        Args:
            target_file (pysam.AlignmentFile) : Target file
            consensus (bool) : write consensus
            no_source_reads (bool) : while in consensus mode, don't write original reads
            consensus_read_callback (function) : this function is called with every consensus read as an arguments
            consensus_read_callback_kwargs (dict) : arguments to pass to the callback function
        """
        if consensus:
            reads = self.deduplicate_majority(target_file,
                                              f'molecule_{uuid4()}' if consensus_name is None else consensus_name)
            if consensus_read_callback is not None:
                if consensus_read_callback_kwargs is not None:
                    consensus_read_callback(reads, **consensus_read_callback_kwargs)
                else:
                    consensus_read_callback(reads)

            for read in reads:
                target_file.write(read)

            if not no_source_reads:
                for read in self.iter_reads():
                    read.is_duplicate=True
                for fragment in self:
                    fragment.write_pysam(target_file)

        else:
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
                                  reads=None
                                  ):
        """Set methylation call tags given a methylation dictionary

        This method sets multiple tags in every read associated to the molecule.
        The tags being set are the bismark_call_tag, every aligned base is annotated
        with a zZxXhH or ".", and a tag for both the total methylated C's and unmethylated C's

        Args:
            call_dict (dict) : Dictionary containing bismark calls (chrom,pos) :
                        {'context':letter,'reference_base': letter   , 'consensus': letter, optiona: 'qual': pred_score (int) }

            bismark_call_tag (str) : tag to write bismark call string

            total_methylated_tag (str) : tag to write total methylated bases

            total_unmethylated_tag (str) : tag to write total unmethylated bases

            reads (iterable) : reads to write the tags to, when not supplied, the tags are written to all associated reads
        Returns:
            can_be_yielded (bool)
        """
        self.methylation_call_dict = call_dict

        # molecule_XM dictionary containing count of contexts
        molecule_XM = Counter(
            list(
                d.get(
                    'context',
                    '.') for d in self.methylation_call_dict.values()))
        # Contruct XM strings
        if reads is None:
            reads = self.iter_reads()
        for read in reads:

            bis_met_call_string = ''.join([
                call_dict.get(
                    (read.reference_name, rpos), {}).get('context', '.')
                # Obtain all aligned positions from the call dict
                # iterate all positions in the alignment
                for qpos, rpos in read.get_aligned_pairs(matches_only=True)
                if qpos is not None and rpos is not None])
            # make sure to ignore non matching positions ? is this neccesary?

            read.set_tag(
                # Write the methylation tag to the read
                bismark_call_tag,
                bis_met_call_string
            )

            # Set total methylated bases
            read.set_tag(
                total_methylated_tag,
                molecule_XM['Z'] + molecule_XM['X'] + molecule_XM['H'])

            # Set total unmethylated bases
            read.set_tag(
                total_unmethylated_tag,
                molecule_XM['z'] + molecule_XM['x'] + molecule_XM['h'])

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


            # Set XR (Read conversion string)
            # @todo: this is TAPS specific, inneficient, ugly
            try:
                fwd = 0
                rev = 0
                for (qpos, rpos, ref_base), call in zip(
                    read.get_aligned_pairs(matches_only=True,with_seq=True),
                    bis_met_call_string):
                    qbase = read.query_sequence[qpos]
                    if call.isupper():
                        if qbase=='A':
                            rev+=1
                        elif qbase=='T':
                            fwd+=1

                # Set XG (genome conversion string)
                if rev>fwd:
                    read.set_tag('XR','CT')
                    read.set_tag('XG','GA')
                else:
                    read.set_tag('XR','CT')
                    read.set_tag('XG','CT')
            except ValueError:
                # Dont set the tag
                pass

    def set_meta(self, tag, value):
        """Set meta information to all fragments

        Args:
            tag (str):
                2 letter tag
            value: value to set

        """
        for f in self:
            f.set_meta(tag, value)

    def __getitem__(self, index):
        """Obtain a fragment belonging to this molecule.

        Args:
            index (int):
                index of the fragment [0 ,1 , 2 ..]

        Returns:
            fragment (singlecellmultiomics.fragment.Fragment)
        """
        return self.fragments[index]

    def get_alignment_stats(self):
        """Get dictionary containing mean amount clip/insert/deletion/matches per base

        Returns:
            cigar_stats(dict): dictionary {
                clips_per_bp(int),
                deletions_per_bp(int),
                matches_per_bp(int),
                inserts_per_bp(int),
                alternative_hits_per_read(int),

                }
        """
        clips = 0
        matches = 0
        inserts = 0
        deletions = 0
        totalbases = 0
        total_reads = 0
        total_alts = 0
        for read in self.iter_reads():
            totalbases += read.query_length
            total_reads += 1
            for operation, amount in read.cigartuples:
                if operation == 4:
                    clips += amount
                elif operation == 2:
                    deletions += amount
                elif operation == 0:
                    matches += amount
                elif operation == 1:
                    inserts += amount
            if read.has_tag('XA'):
                total_alts += len(read.get_tag('XA').split(';'))

        clips_per_bp = clips / totalbases
        inserts_per_bp = inserts / totalbases
        deletions_per_bp = deletions / totalbases
        matches_per_bp = matches / totalbases

        alt_per_read = total_alts / total_reads

        return {
            'clips_per_bp': clips_per_bp,
            'inserts_per_bp': inserts_per_bp,
            'deletions_per_bp': deletions_per_bp,
            'matches_per_bp': matches_per_bp,
            'alt_per_read': alt_per_read,
            'total_bases':totalbases,
            'total_reads':total_reads,
        }

    def get_mean_cycle(
            self,
            chromosome,
            position,
            base=None,
            not_base=None):
        """Get the mean cycle at the supplied coordinate and base-call

        Args:
            chromosome (str)
            position (int)
            base (str) : select only reads with this base
            not_base(str) : select only reads without this base

        Returns:
            mean_cycles (tuple): mean cycle for R1 and R2
        """
        assert (base is not None or not_base is not None), "Supply base or not_base"

        cycles_R1 = []
        cycles_R2 = []
        for read in self.iter_reads():

            if read is None or read.reference_name != chromosome:
                continue


            for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                    read, with_seq=False):

                if query_pos is None or ref_pos != position:
                    continue

                if not_base is not None and read.seq[query_pos] == not_base:
                    continue
                if base is not None and read.seq[query_pos] != base:
                    continue

                if read.is_read2:
                    cycles_R2.append(cycle)
                else:
                    cycles_R1.append(cycle)
        if len(cycles_R2) == 0 and len(cycles_R1)==0:
            raise IndexError(
                "There are no observations if the supplied base/location combination")
        return (np.mean(cycles_R1) if len(cycles_R1) else np.nan),  (np.mean(cycles_R2) if len(cycles_R2) else np.nan)

    def get_mean_base_quality(
                self,
                chromosome,
                position,
                base=None,
                not_base=None):
            """Get the mean phred score at the supplied coordinate and base-call

            Args:
                chromosome (str)
                position (int)
                base (str) : select only reads with this base
                not_base(str) : select only reads without this base

            Returns:
                mean_phred_score (float)
            """
            assert (base is not None or not_base is not None), "Supply base or not_base"

            qualities = []
            for read in self.iter_reads():

                if read is None or read.reference_name != chromosome:
                    continue

                for query_pos, ref_pos in read.get_aligned_pairs(
                        with_seq=False, matches_only=True):

                    if query_pos is None or ref_pos != position:
                        continue

                    if not_base is not None and read.seq[query_pos] == not_base:
                        continue
                    if base is not None and read.seq[query_pos] != base:
                        continue

                    qualities.append(ord(read.qual[query_pos]))
            if len(qualities) == 0:
                raise IndexError(
                    "There are no observations if the supplied base/location combination")
            return np.mean(qualities)

    @cached_property
    def allele_likelihoods(self):
        """
        Per allele likelihood

        Returns:
            likelihoods (dict) : {allele_name : likelihood}

        """
        return self.get_allele_likelihoods()[0]

    @property
    def library(self):
        """
        Associated library

        Returns:
           library (str) : Library associated with the first read of this molecule

        """
        for read in self.iter_reads():
            if read.has_tag('LY'):
                return read.get_tag('LY')

    @cached_property
    def allele_probabilities(self):
        """
        Per allele probability

        Returns:
            likelihoods (dict) : {allele_name : prob}

        """
        return likelihood_to_prob( self.get_allele_likelihoods()[0] )


    @cached_property
    def allele_confidence(self) -> int:
        """
        Returns
            confidence(int) : a phred scalled confidence value for the allele
            assignment, returns zero when no allele is associated to the molecule
        """
        l = self.allele_probabilities
        if l is None or len(l) == 0 :
            return 0
        return int(prob_to_phred( Counter(l).most_common(1)[0][1] ))

    @cached_property
    def base_confidences(self):
        return self.get_base_confidence_dict()

    @cached_property
    def base_likelihoods(self):
        return {(chrom, pos):base_probabilities_to_likelihood(probs) for (chrom, pos),probs in self.base_confidences.items()}

    @cached_property
    def base_probabilities(self):
        # Optimization which is equal to {location:likelihood_to_prob(liks) for location,liks in self.base_likelihoods.items()}
        obs = {}
        for read in self.iter_reads():
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                qbase = read.seq[qpos]
                qqual = read.query_qualities[qpos]
                if qbase=='N':
                    continue
                # @ todo reads which span multiple chromosomes
                k = (self.chromosome, rpos)
                p = 1 - np.power(10, -qqual / 10)

                if not k in obs:
                    obs[k] = {}
                if not qbase in obs[k]:
                    obs[k][qbase] = [p,1] # likelihood, n
                    obs[k]['N'] = [1-p,1] # likelihood, n
                else:
                    obs[k][qbase][0] *= p
                    obs[k][qbase][1] += 1

                    obs[k]['N'][0] *= 1-p # likelihood, n
                    obs[k]['N'][1] += 1 # likelihood, n
        # Perform likelihood conversion and convert to probs
        return { location: likelihood_to_prob({
            base:likelihood/np.power(0.25,n-1)
                    for base,(likelihood,n) in base_likelihoods.items() })
                    for location,base_likelihoods in obs.items()}

    ## This is a duplicate of the above but only calculates for allele informative positions
    @cached_property
    def allele_informative_base_probabilities(self):
        # Optimization which is equal to {location:likelihood_to_prob(liks) for location,liks in self.base_likelihoods.items()}
        obs = {}
        for read in self.iter_reads():
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if not self.allele_resolver.has_location( read.reference_name, rpos ):
                    continue
                qbase = read.seq[qpos]
                qqual = read.query_qualities[qpos]
                if qbase=='N':
                    continue
                # @ todo reads which span multiple chromosomes
                k = (self.chromosome, rpos)
                p = 1 - np.power(10, -qqual / 10)

                if not k in obs:
                    obs[k] = {}
                if not qbase in obs[k]:
                    obs[k][qbase] = [p,1] # likelihood, n
                    obs[k]['N'] = [1-p,1] # likelihood, n
                else:
                    obs[k][qbase][0] *= p
                    obs[k][qbase][1] += 1

                    obs[k]['N'][0] *= 1-p # likelihood, n
                    obs[k]['N'][1] += 1 # likelihood, n
        # Perform likelihood conversion and convert to probs
        return { location: likelihood_to_prob({
            base:likelihood/np.power(0.25,n-1)
                    for base,(likelihood,n) in base_likelihoods.items() })
                    for location,base_likelihoods in obs.items()}



    def calculate_allele_likelihoods(self):
        self.aibd = defaultdict(list)
        self.obtained_allele_likelihoods = Counter()  # Allele -> [prob, prob, prob]

        for (chrom, pos), base_probs in self.allele_informative_base_probabilities.items():

            for base, p in base_probs.items():
                if base == 'N':
                    continue

                assoc_alleles = self.allele_resolver.getAllelesAt(chrom, pos, base)
                if assoc_alleles is not None and len(assoc_alleles) == 1:
                    allele = list(assoc_alleles)[0]
                    self.obtained_allele_likelihoods[allele] += p

                    self.aibd[allele].append((chrom, pos, base, p))



    def get_allele_likelihoods(self,):
        """Obtain the allele(s) this molecule maps to

        Args:
            allele_resolver(singlecellmultiomics.alleleTools.AlleleResolver)  : resolver used
            return_allele_informative_base_dict(bool) : return dictionary containing the bases used for allele determination
            defaultdict(list,
            {'allele1': [
              ('chr18', 410937, 'T'),
              ('chr18', 410943, 'G'),
              ('chr18', 410996, 'G'),
              ('chr18', 411068, 'A')]})

        Returns:
            { 'allele_a': likelihood, 'allele_b':likelihood }
        """
        if self.obtained_allele_likelihoods is None:
            self.calculate_allele_likelihoods()

        return self.obtained_allele_likelihoods, self.aibd



    def get_allele(
            self,
            allele_resolver=None,
            return_allele_informative_base_dict=False):
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
                raise ValueError(
                    "Supply allele resolver or set it to molecule.allele_resolver")

        alleles = set()
        if return_allele_informative_base_dict:
            aibd = defaultdict(list)
        try:
            for (chrom, pos), base in self.get_consensus(
                    base_obs=self.get_base_observation_dict_NOREF()).items():
                c = allele_resolver.getAllelesAt(chrom, pos, base)
                if c is not None and len(c) == 1:
                    alleles.update(c)
                    if return_allele_informative_base_dict:
                        aibd[list(c)[0]].append((chrom, pos, base))

        except Exception as e:
            if return_allele_informative_base_dict:
                return dict()
            else:
                return {}

        if return_allele_informative_base_dict:
            return aibd
        return alleles

    def write_allele_phasing_information_tag(
            self, allele_resolver=None, tag='ap', reads=None):
        """
        Write allele phasing information to ap tag

        For every associated read a tag wil be written containing:
        chromosome,postion,base,allele_name|chromosome,postion,base,allele_name|...
        for all variants found by the AlleleResolver
        """
        if reads is None:
            reads = self.iter_reads()

        use_likelihood = (self.allele_assingment_method==1)

        if not use_likelihood:
            haplotype = self.get_allele(
                return_allele_informative_base_dict=True,
                allele_resolver=allele_resolver)

            phased_locations = [
                (allele, chromosome, position, base)
                for allele, bps in haplotype.items()
                for chromosome, position, base in bps]

            phase_str = '|'.join(
                [
                    f'{chromosome},{position},{base},{allele}' for allele,
                                                                   chromosome,
                                                                   position,
                                                                   base in phased_locations])
        else:

            allele_likelihoods, aibd = self.get_allele_likelihoods()
            allele_likelihoods = likelihood_to_prob(allele_likelihoods)

            phased_locations = [
                (allele, chromosome, position, base, confidence)
                for allele, bps in aibd.items()
                for chromosome, position, base, confidence in bps]

            phase_str = '|'.join(
                [
                    f'{chromosome},{position},{base},{allele},{ prob_to_phred(confidence) }' for allele,
                                                                   chromosome,
                                                                   position,
                                                                   base,
                                                                   confidence in phased_locations])



        if len(phase_str) > 0:
            for read in reads:
                read.set_tag(tag, phase_str)
                if use_likelihood:
                    read.set_tag('al', self.allele_confidence)

    def get_base_observation_dict_NOREF(self, allow_N=False):
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

        base_obs = defaultdict(Counter)

        used = 0  # some alignments yielded valid calls
        ignored = 0
        for fragment in self:
            _, start, end = fragment.span
            for read in fragment:
                if read is None:
                    continue

                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                        read, with_seq=False):

                    if query_pos is None or ref_pos is None or ref_pos < start or ref_pos > end:
                        continue
                    query_base = read.seq[query_pos]
                    if query_base == 'N' and not allow_N:
                        continue
                    base_obs[(read.reference_name, ref_pos)][query_base] += 1

        if used == 0 and ignored > 0:
            raise ValueError('Could not extract any safe data from molecule')

        return base_obs

    def get_base_observation_dict(self, return_refbases=False, allow_N=False,
        allow_unsafe=True, one_call_per_frag=False, min_cycle_r1=None,
         max_cycle_r1=None, min_cycle_r2=None, max_cycle_r2=None, use_cache=True, min_bq=None):
        '''
        Obtain observed bases at reference aligned locations

        Args:
            return_refbases ( bool ):
                return both observed bases and reference bases
            allow_N (bool): Keep N base calls in observations

            min_cycle_r1(int) : Exclude read 1 base calls with a cycle smaller than this value (excludes bases which are trimmed before mapping)

            max_cycle_r1(int) : Exclude read 1 base calls with a cycle larger than this value (excludes bases which are trimmed before mapping)

            min_cycle_r2(int) : Exclude read 2 base calls with a cycle smaller than this value (excludes bases which are trimmed before mapping)

            max_cycle_r2(int) : Exclude read 2 base calls with a cycle larger than this value (excludes bases which are trimmed before mapping)


        Returns:
            { genome_location (tuple) : base (string) : obs (int) }
            and
            { genome_location (tuple) : base (string) if return_refbases is True }
        '''

        # Check if cached is available
        if use_cache:
            if self.saved_base_obs is not None :
                if not return_refbases:
                    return self.saved_base_obs[0]
                else:
                    if self.saved_base_obs[1] is not None:
                        return self.saved_base_obs

        base_obs = defaultdict(Counter)

        ref_bases = {}
        used = 0  # some alignments yielded valid calls
        ignored = 0
        error = None
        for fragment in self:
            _, start, end = fragment.span

            used += 1

            if one_call_per_frag:
                frag_location_obs = set()

            for read in fragment:
                if read is None:
                    continue

                if allow_unsafe:
                    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                        if query_pos is None or ref_pos is None:  # or ref_pos < start or ref_pos > end:
                            continue

                        query_base = read.seq[query_pos]
                        # query_qual = read.qual[query_pos]
                        if min_bq is not None and read.query_qualities[query_pos]<min_bq:
                            continue

                        if query_base == 'N':
                            continue

                        k = (read.reference_name, ref_pos)

                        if one_call_per_frag:
                            if k in frag_location_obs:
                                continue
                            frag_location_obs.add(k)

                        base_obs[k][query_base] += 1

                        if return_refbases:
                            ref_bases[(read.reference_name, ref_pos)
                            ] = ref_base.upper()


                else:
                    for cycle, query_pos, ref_pos, ref_base in pysamiterators.iterators.ReadCycleIterator(
                            read, with_seq=True, reference=self.reference):

                        if query_pos is None or ref_pos is None:  # or ref_pos < start or ref_pos > end:
                            continue

                        # Verify cycle filters:
                        if (not read.is_paired or read.is_read1) and (
                                ( min_cycle_r1 is not None and cycle <  min_cycle_r1 ) or
                                ( max_cycle_r1 is not None and  cycle >  max_cycle_r1 )):
                            continue

                        if (read.is_paired and read.is_read2) and (
                                ( min_cycle_r2 is not None and cycle <  min_cycle_r2 ) or
                                ( max_cycle_r2 is not None and cycle >  max_cycle_r2 )):
                            continue

                        query_base = read.seq[query_pos]
                        # Skip bases with low bq:
                        if min_bq is not None and read.query_qualities[query_pos]<min_bq:
                            continue

                        if query_base == 'N':
                            continue

                        k = (read.reference_name, ref_pos)
                        if one_call_per_frag:
                            if k in frag_location_obs:
                                continue
                            frag_location_obs.add(k)

                        base_obs[(read.reference_name, ref_pos)][query_base] += 1

                        if return_refbases:
                            ref_bases[(read.reference_name, ref_pos)
                            ] = ref_base.upper()

        if used == 0 and ignored > 0:
            raise ValueError(f'Could not extract any safe data from molecule {error}')

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

        base_obs, ref_bases = self.get_base_observation_dict(
            return_refbases=True)
        for location, obs in base_obs.items():
            if ignore_locations is not None and location in ignore_locations:
                continue

            if location in ref_bases:
                ref = ref_bases[location]
                if ref not in 'ACTG':  # don't count weird bases in the reference @warn
                    continue
                matches += obs[ref]
                mismatches += sum((base_obs for base,
                                                base_obs in obs.most_common() if base != ref))

        return matches, mismatches

    def get_consensus(self,
                    dove_safe: bool = False,
                    only_include_refbase: str = None,
                    allow_N=False,
                    with_probs_and_obs=False, **get_consensus_dictionaries_kwargs):
        """
        Obtain consensus base-calls for the molecule
        """

        if allow_N:
            raise NotImplementedError()

        consensii = defaultdict(consensii_default_vector)  # location -> obs (A,C,G,T,N)
        if with_probs_and_obs:
            phred_scores = defaultdict(lambda:defaultdict(list))
        for fragment in self:
            if dove_safe and not fragment.has_R2() or not fragment.has_R1():
                continue

            try:
                for position, (q_base, phred_score) in fragment.get_consensus(
                        dove_safe=dove_safe,only_include_refbase=only_include_refbase,
                        **get_consensus_dictionaries_kwargs).items():

                    if q_base == 'N':
                        continue
                    #    consensii[position][4] += phred_score
                    # else:
                    if with_probs_and_obs:
                        phred_scores[position][q_base].append(phred_score)

                    consensii[position]['ACGTN'.index(q_base)] += 1
            except ValueError as e:
                # For example: ValueError('This method only works for inwards facing reads')
                pass
        if len(consensii)==0:
            if with_probs_and_obs:
                return dict(),None,None
            else:
                return dict()

        locations = np.empty(len(consensii), dtype=object)
        locations[:] = sorted(list(consensii.keys()))

        v = np.vstack([ consensii[location] for location in locations])
        majority_base_indices = np.argmax(v, axis=1)

        # Check if there is ties, this result in multiple hits for argmax (majority_base_indices),
        # such a situtation is of course terrible and should be dropped
        proper = (v == v[np.arange(v.shape[0]), majority_base_indices][:, np.newaxis]).sum(1) == 1

        if with_probs_and_obs:
            return (
                dict(zip(locations[proper], ['ACGTN'[idx] for idx in majority_base_indices[proper]])),
                phred_scores,
                consensii
            )
        else:
            return  dict(zip(locations[proper], ['ACGTN'[idx] for idx in majority_base_indices[proper]]))


    def get_consensus_old(
            self,
            base_obs=None,
            classifier=None,
            store_consensus=True,
            reuse_cached_consensus=True,
            allow_unsafe=False,
            allow_N=False):
        """Get dictionary containing consensus calls in respect to reference.
        By default mayority voting is used to determine the consensus base. If a classifier is supplied the classifier is used to determine the consensus base.

        Args:
            base_obs (defaultdict(Counter)) :
                { genome_location (tuple) : base (string) : obs (int) }

            classifier : fitted classifier to use for consensus calling. When no classifier is provided the consensus is determined by majority voting
            store_consensus (bool) : Store the generated consensus for re-use

        Returns:
            consensus (dict)  :  {location : base}
        """
        consensus = {}  # postion -> base , key is not set when not decided

        if classifier is not None:

            if reuse_cached_consensus and hasattr(
                    self, 'classifier_consensus') and self.classifier_consensus is not None:
                return self.classifier_consensus

            features, reference_bases, CIGAR, alignment_start, alignment_end = self.get_base_calling_feature_matrix_spaced(
                True)

            if features is None:
                # We cannot determine the consensus as there are no features...
                return dict()

            base_calling_probs = classifier.predict_proba(features)
            predicted_sequence = ['ACGT'[i] for i in np.argmax(base_calling_probs, 1)]

            reference_sequence = ''.join(
                [base for chrom, pos, base in reference_bases])

            phred_scores = np.rint(
                -10 * np.log10(np.clip(1 - base_calling_probs.max(1),
                                       0.000000001,
                                       0.999999999)
                               )).astype('B')

            consensus = {(chrom, pos): consensus_base for (
                                                              chrom, pos, ref_base), consensus_base in
                         zip(reference_bases, predicted_sequence)}

            if store_consensus:
                self.classifier_consensus = consensus
                self.classifier_phred_scores = phred_scores
            return consensus

        if base_obs is None:
            try:
                base_obs, ref_bases = self.get_base_observation_dict(
                    return_refbases=True, allow_N=allow_N, allow_unsafe=allow_unsafe)
            except ValueError as e:
                # We cannot determine safe regions
                raise

        for location, obs in base_obs.items():
            votes = obs.most_common()
            if len(votes) == 1 or votes[1][1] < votes[0][1]:
                consensus[location] = votes[0][0]

        if store_consensus:
            self.majority_consensus = consensus

        return consensus

    def get_consensus_base(self, contig, position, classifier=None):
        """Obtain base call at single position of the molecule

        Args:
            contig (str) : contig to extract base call from

            position (int) : position to extract base call from (zero based)

            classifier (obj) : base calling classifier

        Returns:

            base_call (str) : base call, or None when no base call could be made
        """

        try:
            c = self.get_consensus(classifier)
        except ValueError:
            return None
        return c.get((contig, position), None)

    # when enabled other calls (non ref non alt will be set None)
    def check_variants(self, variants, exclude_other_calls=True):
        """Verify variants in molecule

        Args:
            variants (pysam.VariantFile) : Variant file handle to extract variants from

        Returns:
            dict (defaultdict( Counter )) : { (chrom,pos) : ( call (str) ): observations  (int) }
        """
        variant_dict = {}
        for variant in variants.fetch(
                self.chromosome,
                self.spanStart,
                self.spanEnd):
            variant_dict[(variant.chrom, variant.pos - 1)
            ] = (variant.ref, variant.alts)

        variant_calls = defaultdict(Counter)
        for fragment in self:

            _, start, end = fragment.span
            for read in fragment:
                if read is None:
                    continue

                for cycle, query_pos, ref_pos in pysamiterators.iterators.ReadCycleIterator(
                        read):

                    if query_pos is None or ref_pos is None or ref_pos < start or ref_pos > end:
                        continue
                    query_base = read.seq[query_pos]

                    k = (read.reference_name, ref_pos)
                    if k in variant_dict:
                        call = None
                        ref, alts = variant_dict[k]
                        if query_base == ref:
                            call = ('ref', query_base)
                        elif query_base in alts:
                            call = ('alt', query_base)

                        if not exclude_other_calls or call is not None:
                            variant_calls[k][call] += 1

        return variant_calls

    def get_aligned_reference_bases_dict(self):
        """Get dictionary containing all reference bases to which this molecule aligns
        Returns:
            aligned_reference_positions (dict) :  { (chrom,pos) : 'A', (chrom,pos):'T', .. }
        """
        aligned_reference_positions = {}
        for read in self.iter_reads():
            for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                    with_seq=True, matches_only=True):
                aligned_reference_positions[(
                    read.reference_name, ref_pos)] = ref_base.upper()
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

    @property
    def span_len(self):
        return abs(self.spanEnd - self.spanStart)

    def get_methylated_count(self, context=3):
        """Get the total amount of methylated bases

        Args:
            context (int) : 3 or 4 base context

        Returns:
            r (Counter) : sum of methylated bases in contexts
        """

        r = Counter()

    def get_html(
            self,
            reference=None,
            consensus=None,
            show_reference_sequence=True,
            show_consensus_sequence=True,
            reference_bases=None):
        """Get html representation of the molecule
        Returns:
            html_rep(str) : Html representation of the molecule
        """

        html = f"""<h3>{self.chromosome}:{self.spanStart}-{self.spanEnd}
            sample:{self.get_sample()}  {'valid molecule' if self[0].is_valid() else 'Non valid molecule'}</h3>
            <h5>UMI:{self.get_umi()} Mapping qual:{round(self.get_mean_mapping_qual(), 1)} Cut loc: {"%s:%s" % self[0].get_site_location()} </h5>
            <div style="white-space:nowrap; font-family:monospace; color:#888">"""
        # undigested:{self.get_undigested_site_count()}
        consensus = self.get_consensus()

        # Obtain reference bases dictionary:
        if reference_bases is None:
            if reference is None:
                reference_bases = self.get_aligned_reference_bases_dict()

            else:
                # obtain reference_bases from reference file
                raise NotImplementedError()

        for fragment in itertools.chain(*self.get_rt_reactions().values()):
            html += f'<h5>{fragment.get_R1().query_name}</h5>'
            for readid, read in [
                (1, fragment.get_R1()),
                (2, fragment.get_R2())]:  # go over R1 and R2:
                # This is just the sequence:
                if read is None:
                    continue
                html += fragment.get_html(
                    self.chromosome,
                    self.spanStart,
                    self.spanEnd,
                    show_read1=(readid == 1),
                    show_read2=(readid == 2)
                ) + '<br />'

        # Obtain reference sequence and consensus sequence:
        if consensus is None:
            consensus = self.get_consensus()

        span_len = self.spanEnd - self.spanStart
        visualized = ['.'] * span_len
        reference_vis = ['?'] * span_len
        for location, query_base in consensus.items():
            try:
                if reference_bases is None or reference_bases.get(
                        location, '?') == query_base:
                    visualized[location[1] - self.spanStart] = query_base
                    if reference_bases is not None:
                        # or reference_bases.get(location,'?')
                        reference_vis[location[1] -
                                      self.spanStart] = query_base
                else:
                    visualized[location[1] -
                               self.spanStart] = style_str(query_base, color='red', weight=800)
                    if reference_bases is not None:
                        reference_vis[location[1] - self.spanStart] = style_str(
                            reference_bases.get(location, '?'), color='black', weight=800)
            except IndexError as e:
                pass  # Tried to visualize a base outside view

        if show_consensus_sequence:
            html += ''.join(visualized) + '<br />'

        if show_reference_sequence:
            html += ''.join(reference_vis) + '<br />'

        html += "</div>"
        return html

    def get_methylation_dict(self):
        """Obtain methylation dictionary

        Returns:
            methylated_positions (Counter):
                (read.reference_name, rpos) : times seen methylated

            methylated_state (dict):
                {(read.reference_name, rpos) : 1/0/-1 }
                1 for methylated
                0 for unmethylated
                -1 for unknown

        """
        methylated_positions = Counter()  # chrom-pos->count
        methylated_state = dict()  # chrom-pos->1, 0, -1
        for fragment in self:
            for read in fragment:
                if read is None or not read.has_tag('XM'):
                    continue
                methylation_status_string = read.get_tag('XM')
                i = 0
                for qpos, rpos, ref_base in read.get_aligned_pairs(
                        with_seq=True):
                    if qpos is None:
                        continue
                    if ref_base is None:
                        continue
                    if rpos is None:
                        continue
                    methylation_status = methylation_status_string[i]
                    if methylation_status.isupper():
                        methylated_positions[(read.reference_name, rpos)] += 1
                        if methylated_state.get(
                                (read.reference_name, rpos), 1) == 1:
                            methylated_state[(read.reference_name, rpos)] = 1
                        else:
                            methylated_state[(read.reference_name, rpos)] = -1
                    else:
                        if methylation_status == '.':
                            pass
                        else:
                            if methylated_state.get(
                                    (read.reference_name, rpos), 0) == 0:
                                methylated_state[(
                                    read.reference_name, rpos)] = 0
                            else:
                                methylated_state[(
                                    read.reference_name, rpos)] = -1
                    i += 1
        return methylated_positions, methylated_state


    def _get_allele_from_reads(self) -> str:
        """
        Obtain associated allele based on the associated reads of the molecule

        """
        allele = None
        for frag in self:
            for read in frag:
                if read is None or not read.has_tag('DA'):
                    continue
                allele = read.get_tag('DA')
                return allele
        return None
