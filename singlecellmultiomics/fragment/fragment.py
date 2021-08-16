import itertools
from singlecellmultiomics.utils.sequtils import hamming_distance, get_consensus_dictionaries, pick_best_base_call
import pysamiterators.iterators
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
from singlecellmultiomics.utils import style_str
from singlecellmultiomics.bamProcessing import get_read_group_from_read
from singlecellmultiomics.features import FeatureAnnotatedObject

complement = str.maketrans('ATCGN', 'TAGCN')


class Fragment():
    """
    This class holds 1 or more reads which are derived from the same cluster

    Example:
        Generate a Fragment with a single associated read::

            >>> from singlecellmultiomics.molecule import Molecule
            >>> from singlecellmultiomics.fragment import Fragment
            >>> import pysam
            >>> read = pysam.AlignedSegment()
            >>> read.reference_start = 30
            >>> read.query_name = 'R1'
            >>> read.mapping_quality = 30
            >>> read.set_tag('SM','CELL_1') # The sample to which the sample belongs is extracted from the SM tag
            >>> read.set_tag('RX','CAT') # The UMI is extracted from the RX tag
            >>> read.query_sequence = "CATGTATCCGGGCTTAA"
            >>> read.query_qualities = [30] * len(read.query_sequence)
            >>> read.cigarstring = f'{len(read.query_sequence)}M'
            >>> Fragment([read])
            Fragment:
                sample:CELL_1
                umi:CAT
                span:None 30-47
                strand:+
                has R1: yes
                has R2: no
                randomer trimmed: no

    Warning:
        Make sure the RX and SM tags of the read are set! If these are encoded
        in the read name, use singlecellmultiomics.universalBamTagger.customreads
        for conversion.

    """

    def __init__(self, reads,
                 assignment_radius: int = 0,
                 umi_hamming_distance: int = 1,
                 R1_primer_length: int = 0,
                 R2_primer_length: int = 6,
                 tag_definitions: list = None,
                 max_fragment_size: int = None,
                 mapping_dir=(False, True),
                 max_NUC_stretch: int = None,
                 read_group_format: int = 0,  # R1 forward, R2 reverse
                 library_name: str = None, # Overwrites the library name
                 single_end: bool = False
                 ):
        """
        Initialise Fragment

            Args:
                reads( list ):
                    list containing pysam AlignedSegment reads

                assignment_radius( int ):
                    Assignment radius for molecule resolver

                umi_hamming_distance( int ):
                    umi_hamming_distance distance threshold for associating multiple fragments to a molecule

                R1_primer_length(int):
                    length of R1 primer, these bases are not taken into account when calculating a consensus

                R2_primer_length(int):
                    length of R2 primer, these bases are not taken into account when calculating a consensus, this length is auto-detected/ overwritten when the rS tag is set


                tag_definitions(list):
                    sam TagDefinitions

                max_fragment_size (int):
                    When specified, fragments with a fragment size larger than the specified value are flagged as qcfail

                mapping_dir( tuple ):
                    expected mapping direction of mates,
                    True for reverse, False for forward. This parameter is used in dovetail detection

                read_group_format(int) : see get_read_group()
            Returns:
                html(string) : html representation of the fragment
        """
        self.R1_primer_length = R1_primer_length
        self.R2_primer_length = R2_primer_length
        if tag_definitions is None:
            tag_definitions = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TagDefinitions
        self.tag_definitions = tag_definitions
        self.assignment_radius = assignment_radius
        self.umi_hamming_distance = umi_hamming_distance
        self.mapping_dir = mapping_dir
        self.reads = reads
        self.strand = None
        self.meta = {}  # dictionary of meta data
        self.is_mapped = None
        # Check for multimapping
        self.is_multimapped = True
        self.mapping_quality = 0
        self.match_hash = None
        self.safe_span = None  # wether the span of the fragment could be determined
        self.unsafe_trimmed = False  # wether primers have been trimmed off
        self.random_primer_sequence = None
        self.max_fragment_size = max_fragment_size
        self.read_group_format = read_group_format
        self.max_NUC_stretch = max_NUC_stretch
        self.qcfail = False
        self.single_end = single_end

        # Span:\
        self.span = [None, None, None]

        self.umi = None

        # Force R1=read1 R2=read2:

        for i, read in enumerate(self.reads):

            if read is None:
                continue

            if self.max_NUC_stretch is not None and (
                self.max_NUC_stretch*'A' in read.seq or
                self.max_NUC_stretch*'G' in read.seq or
                self.max_NUC_stretch*'T' in read.seq or
                self.max_NUC_stretch*'C' in read.seq):

                self.set_rejection_reason('HomoPolymer', set_qcfail=True)
                self.qcfail=True
                break

            if read.is_qcfail:
                self.qcfail = True

            if read.has_tag('rS'):
                self.random_primer_sequence = read.get_tag('rS')
                self.unsafe_trimmed = True
                self.R2_primer_length = 0 # Set R2 primer length to zero, as the primer has been removed

            if not read.is_unmapped:
                self.is_mapped = True
            if read.mapping_quality > 0:
                self.is_multimapped = False
            self.mapping_quality = max(
                self.mapping_quality, read.mapping_quality)

            if i == 0:
                if read.is_read2:
                    if not read.is_qcfail:
                        raise ValueError('Supply first R1 then R2')
                read.is_read1 = True
                read.is_read2 = False
            elif i == 1:
                read.is_read1 = False
                read.is_read2 = True
        self.set_sample(library_name=library_name)
        self.update_umi()
        if self.qcfail:
            return
        self.set_strand(self.identify_strand())
        self.update_span()


    def write_tensor(self, chromosome=None, span_start=None, span_end=None, height=30, index_start=0,
                     base_content_table=None,
                     base_mismatches_table=None,
                     base_indel_table=None,
                     base_qual_table=None,
                     base_clip_table=None,
                     mask_reference_bases=None,
                     # set with reference bases to mask (converted to N )
                     # (chrom,pos)
                     reference=None,
                     skip_missing_reads=False
                     ):
        """
        Write tensor representation of the fragment to supplied 2d arrays

        Args:
            chromosome( str ):
                chromosome to view

            span_start( int ):
                first base to show

            span_end( int ):
                last base to show

            height (int) :
                height of the tensor (reads)

            index_start(int): start writing tensor at this row

            base_content_table(np.array) : 2d array to write base contents to

            base_indel_table(np.array) : 2d array to write indel information to

            base_content_table(np.array) : 2d array to write base contents to

            base_content_table(np.array) : 2d array to write base contents to

            mask_reference_bases(set): mask reference bases in this set with a N set( (chrom,pos), ... )

            reference(pysam.FastaFile) :  Handle to reference file to use instead of MD tag. If None: MD tag is used.

            skip_missing_reads (bool) : when enabled only existing (non None) reads are added to the tensor. Use this option when mapping single-end

        Returns:
            None

        """
        self.base_mapper = {}
        for i, (strand, base) in enumerate(
                itertools.product((True, False), 'NACTG-')):
            self.base_mapper[(strand, base)] = i + 2

        if chromosome is None and span_start is None and span_end is None:
            chromosome, span_start, span_end = self.get_span()

        span_len = span_end - span_start

        reads = self
        last_ref = None
        used_reads = 0
        for ri, read in enumerate(reads):

            if read is None:
                continue

            alignment_operations = list(itertools.chain(
                *([operation] * l for operation, l in read.cigartuples if operation != 2)))
            print(alignment_operations)
            # Obtain the amount of operations before the alignment start
            operations_before_alignment_start = 0
            for query_pos, ref_pos in read.get_aligned_pairs():
                if ref_pos is None:
                    operations_before_alignment_start += 1
                else:
                    break
            initial_base_offset = read.reference_start - operations_before_alignment_start

            # Obtain row index in the tensor
            if skip_missing_reads:
                row_index = (used_reads + index_start) % height
            else:
                row_index = (ri + index_start) % height

            # Where in the reference are we?
            ref_pointer = read.reference_start

            alignment_started = False
            for operation, (cycle, query_pos, ref_pos, ref_base) in zip(
                alignment_operations,
                pysamiterators.iterators.ReadCycleIterator(
                    read,
                    with_seq=True,
                    reference=reference)):

                # Mask locations from mask_reference_bases
                if mask_reference_bases is not None and ref_pos is not None and (
                        chromosome, ref_pos) in mask_reference_bases:
                    ref_base = 'N'

                if ref_pos is None:
                    if not alignment_started:
                        ref_pointer = initial_base_offset + i
                    else:
                        ref_pointer += 1
                    if operation == 4:
                        try:
                            query_base = read.seq[query_pos]
                            base_clip_table[row_index, ref_pointer -
                                            span_start] = self.base_mapper[(read.is_reverse, query_base)]
                        except IndexError:
                            pass
                    continue
                ref_pointer = ref_pos

                if (ref_pos - span_start) < 0 or (ref_pos -
                                                  span_start) > (span_len - 1):
                    continue

                ref_base = ref_base.upper()
                if query_pos is None:
                    query_base = '-'
                    qual = 0
                    base_indel_table[row_index, ref_pos - span_start] = 1
                else:
                    query_base = read.seq[query_pos]
                    qual = ord(read.qual[query_pos])

                if ref_pos is not None:
                    base_qual_table[row_index, ref_pos - span_start] = qual
                    if query_base != ref_base:
                        base_content_table[row_index, ref_pos -
                                           span_start] = self.base_mapper[(read.is_reverse, query_base)]
                    else:
                        base_mismatches_table[row_index,
                                              ref_pos - span_start] = 0.2
                    base_content_table[row_index, ref_pos -
                                       span_start] = self.base_mapper[(read.is_reverse, query_base)]

            used_reads += 1
        if skip_missing_reads:
            return used_reads + index_start + 1
        else:
            return ri + index_start + 1

    def get_read_group(self, with_attr_dict=False) -> str:

        r = None
        for read in self.reads:
            if read is None:
                continue

            r = get_read_group_from_read(read,
                    format=self.read_group_format,
                    with_attr_dict=with_attr_dict)

            break

        return r


    def set_duplicate(self, is_duplicate):
        """Define this fragment as duplicate, sets the corresponding bam bit flag
        Args:
            value(bool) : is_duplicate
        """
        for read in self.reads:
            if read is not None:
                read.is_duplicate = is_duplicate

    def get_fragment_size(self):
        return abs(self.span[2] - self.span[1])

    def write_tags(self):
        self.set_meta('MQ', self.mapping_quality)
        self.set_meta('MM', self.is_multimapped)
        self.set_meta('RG', self.get_read_group())
        if self.has_valid_span():
            # Write fragment size:
            if self.safe_span:
                self.set_meta('fS', self.get_fragment_size())
                self.set_meta('fe', self.span[1])
                self.set_meta('fs', self.span[2])
            else:
                self.remove_meta('fS')

        else:
            self.set_rejection_reason('FS', set_qcfail=True)

        # Set qcfail bit when the fragment is not valid

        if not self.is_valid():
            for read in self:
                if read is not None:
                    read.is_qcfail = True

    def write_pysam(self, pysam_handle):
        """Write all associated reads to the target file

        Args:
            target_file (pysam.AlignmentFile) : Target file
        """
        self.write_tags()
        for read in self:
            if read is not None:
                pysam_handle.write(read)

    def identify_strand(self):
        # If R2 is rev complement:
        if self.get_R1() is not None:
            # verify if the read has been mapped:
            if self.get_R1().is_unmapped:
                return None
            return self.get_R1().is_reverse
        elif self.get_R2() is not None:
            # verify if the read has been mapped:
            if self.get_R2().is_unmapped:
                return None
            return not self.get_R2().is_reverse
        return None

    def get_random_primer_hash(self):
        """Obtain hash describing the random primer
        this assumes the random primer is on the end of R2 and has a length of self.R2_primer_length
        When the rS tag is set, the value of this tag is used as random primer sequence
        Returns None,None when the random primer cannot be described

        Returns:
            reference_name (str) or None
            reference_start (int) : Int or None
            sequence (str) : Int or None
        """
        R2 = self.get_R2()

        if R2 is None or R2.query_sequence is None:
            return None, None, None

        if self.R2_primer_length == 0 and self.random_primer_sequence is None:
            return None, None, None

        # The read was not mapped
        if R2.is_unmapped:
            # Guess the orientation does not matter
            if self.random_primer_sequence is not None:
                return None, None, self.random_primer_sequence
            return None, None, R2.query_sequence[:self.R2_primer_length]

        if R2.is_reverse:
            global complement
            if self.random_primer_sequence is not None:
                return R2.reference_name, R2.reference_end, self.random_primer_sequence

            return(R2.reference_name, R2.reference_end, R2.query_sequence[-self.R2_primer_length:][::-1].translate(complement))
        else:
            if self.random_primer_sequence is not None:
                return R2.reference_name, R2.reference_start, self.random_primer_sequence
            return(R2.reference_name, R2.reference_start, R2.query_sequence[:self.R2_primer_length])
        raise ValueError()

    def remove_meta(self,key):
        for read in self:
            if read is not None:
                if read.has_tag(key):
                    read.set_tag(key,None)


    def set_meta(self, key, value, as_set=False):
        self.meta[key] = value
        for read in self:
            if read is not None:
                if as_set and read.has_tag(key):

                    data = set(read.get_tag(key).split(','))
                    data.add(value)
                    read.set_tag(key, ','.join(sorted(list(data))))

                else:
                    read.set_tag(key, value)

    def get_R1(self):
        """
        Obtain the AlignedSegment of read 1 of the fragment

        Returns:
            R1 (pysam.AlignedSegment) : Read 1 of the fragment, returns None
                                        when R1 is not mapped
        """
        if len(self.reads) == 0:
            raise IndexError('The fragment has no associated reads')
        return self.reads[0]

    def get_R2(self):
        """
        Obtain the AlignedSegment of read 2 of the fragment

        Returns:
            R1 (pysam.AlignedSegment) : Read 2 of the fragment, returns None
                                        when R2 is not mapped
        """
        if len(self.reads) < 2:
            raise IndexError('The fragment has no associated R2')
        return self.reads[1]

    def has_R1(self):
        """
        Check if the fragment has an associated read 1

        Returns:
            has_r1 (bool)
        """
        if len(self.reads) == 0:
            return False
        return self.reads[0] is not None

    def has_R2(self):
        """
        Check if the fragment has an associated read 2

        Returns:
            has_r2 (bool)
        """
        if len(self.reads) < 2:
            return False
        return self.reads[1] is not None

    @property
    def R1(self):
        return self.reads[0]

    @property
    def R2(self):
        return self.reads[1]

    def get_consensus(self, only_include_refbase: str = None, dove_safe: bool = False, **get_consensus_dictionaries_kwargs) -> dict:
        """
        a dictionary of (reference_pos) : (qbase, quality, reference_base) tuples

        Args:
            only_include_refbase(str) : Only report bases aligned to this reference base, uppercase only
            dove_safe(bool) : Only report bases supported within R1 and R2 start and end coordinates

        Returns:
            consensus(dict) : {reference_position: (qbase, quality)
        """
        r1_consensus, r2_consensus = get_consensus_dictionaries(
                self.R1,
                self.R2,
                only_include_refbase=only_include_refbase,
                dove_safe=dove_safe,
                **get_consensus_dictionaries_kwargs)

        return {
            ref_pos:pick_best_base_call( r1_consensus.get(ref_pos) , r2_consensus.get(ref_pos) )
            for ref_pos in set(r1_consensus.keys()).union(set(r2_consensus.keys()))
        }


    @property
    def estimated_length(self) -> int:
        """
        Obtain the estimated size of the fragment,
        returns None when estimation is not possible
        Takes into account removed bases (R2_primer_length attribute)
        Assumes inwards sequencing orientation, except when self.single_end is set
        """
        if self.single_end:
            if self[0] is None:
                return None
            return self[0].reference_end - self[0].reference_start


        if self.has_R1() and self.has_R2():

            contig = self.R1.reference_name
            if self.R1.is_reverse and not self.R2.is_reverse:
                ##     rp |----- R2 ---->  <--- R1 ----|
                ##    |>----------DISTANCE------------<|
                start, end = self.R2.reference_start - self.R2_primer_length, self.R1.reference_end
                if start<end:
                    return end - start
                else:
                    return None
            elif not self.R1.is_reverse and self.R2.is_reverse:

                ##    |----- R1 ---->  <--- R2 ----|  rp
                ##    |>----------DISTANCE------------<|

                start, end = self.R1.reference_start, self.R2.reference_end + self.R2_primer_length
                if start<end:
                    return end - start
                else:
                    return None
            else:
                return None
        return None


    def update_span(self):
        """
        Update the span (the location the fragment maps to) stored in Fragment

        The span is a single stretch of coordinates on a single contig.
        The contig is determined by the reference_name assocated to the first
        mapped read in self.reads

        This calculation assumes the reads are sequenced inwards and
        dove-tails of the molecule cannot be trusted
        """

        contig = None
        start = None
        end = None


        if self.has_R1() and self.has_R2() and \
           self.R1.reference_start is not None and self.R1.reference_end is not None and \
           self.R2.reference_start is not None and self.R2.reference_end is not None :

            contig = self.R1.reference_name
            if self.R1.is_reverse and not self.R2.is_reverse:
                start, end = self.R2.reference_start, self.R1.reference_end
                self.safe_span = True
            elif not self.R1.is_reverse and self.R2.is_reverse:
                start, end = self.R1.reference_start, self.R2.reference_end
                self.safe_span = True
            else:
                start = min(self.R1.reference_start, self.R2.reference_start)
                end = max(self.R1.reference_start, self.R2.reference_start)
                self.safe_span = False

        elif self.has_R1() and self.R1.reference_start is not None and self.R1.reference_end is not None :
            contig = self.R1.reference_name
            start, end = self.R1.reference_start, self.R1.reference_end
            self.safe_span = False

        elif self.has_R2()  and self.R2.reference_start is not None and self.R2.reference_end is not None :
            contig = self.R2.reference_name
            start, end = self.R2.reference_start, self.R2.reference_end
            self.safe_span = False
        else:

            # Its Sometimes possible that no cigar is set for the alignment, only a start coordinate

            for read in self:
                if read is None:
                    continue
                if len(read.cigar)!=0:
                    raise NotImplementedError("Non implemented span")
                if read.reference_start is not None:
                    start,end = read.reference_start, read.reference_start
                    contig = read.reference_name
                else:
                    raise NotImplementedError("Non implemented span, undefined alignment, and not start coordinate")




        self.span = (contig, start, end)

    @property
    def aligned_length(self):
        return self.span[2]-self.span[1]

    def get_span(self) -> tuple:
        return self.span

    def has_valid_span(self) -> bool:

        # Check if the span could be calculated:
        defined_span = not any([x is None for x in self.get_span()])
        if not defined_span:
            return False

        if self.max_fragment_size is not None and self.get_fragment_size()>self.max_fragment_size:
            return False

        return True

    def set_sample(self, sample: str = None, library_name: str = None):
        """Force sample name or obtain sample name from associated reads"""
        if sample is not None:
            self.sample = sample
        else:
            for read in self.reads:
                if read is not None and read.has_tag('SM'):
                    self.sample = read.get_tag('SM')
                    break

        if library_name is not None:
            sample = self.sample.split('_', 1)[-1]
            self.sample = library_name+'_'+sample
            for read in self.reads:
                if read is not None:
                    read.set_tag('SM', self.sample)
                    read.set_tag('LY', library_name)

    def get_sample(self) -> str:
        """ Obtain the sample name associated with the fragment
        The sample name is extracted from the SM tag of any of the associated reads.

        Returns:
            sample name (str)
        """
        return self.sample

    def __len__(self):
        """Obtain the amount of associated reads to the fragment

        Returns:
            assocoiated reads (int)
        """
        return(sum(read is not None for read in self.reads))

    def set_strand(self, strand: bool):
        """Set mapping strand

        Args:
            strand (bool) : False for Forward, True for reverse
        """
        self.strand = strand

    def get_strand(self):
        """Obtain strand

        Returns:
            strand (bool) : False for Forward, True for reverse
        """
        return self.strand

    def update_umi(self):
        """
        Extract umi from 'RX' tag and store the UMI in the Fragment object
        """
        for read in self.reads:
            if read is not None and read.has_tag('RX'):
                self.umi = read.get_tag('RX')

    def get_umi(self):
        """
        Get the UMI sequence associated to this fragment

        Returns:
            umi(str)
        """
        return self.umi

    def __iter__(self):
        return iter(self.reads)

    def umi_eq(self, other):
        """
        Hamming distance measurement to another Fragment or Molecule,

        Returns :
            is_close (bool) : returns True when the hamming distance between
            the two objects <= umi_hamming_distance

        Example:
            >>> from singlecellmultiomics.molecule import Molecule
            >>> from singlecellmultiomics.fragment import Fragment
            >>> import pysam
            >>> # Create reads (read_A, and read_B), they both belong to the same
            >>> # cell and have 1 hamming distance between their UMI's
            >>> read_A = pysam.AlignedSegment()
            >>> read_A.set_tag('SM','CELL_1') # The sample to which the sample belongs is extracted from the SM tag
            >>> read_A.set_tag('RX','CAT') # The UMI is extracted from the RX tag
            >>> read_B = pysam.AlignedSegment()
            >>> read_B.set_tag('SM','CELL_1') # The sample to which the sample belongs is extracted from the SM tag
            >>> read_B.set_tag('RX','CAG') # The UMI is extracted from the RX tag

            >>> # Create fragment objects for read_A and B:
            >>> frag_A = Fragment([read_A],umi_hamming_distance=0)
            >>> frag_B = Fragment([read_B],umi_hamming_distance=0)
            >>> frag_A.umi_eq(frag_B) # This returns False, the distance is 1, which is higher than 0 (umi_hamming_distance)
            False

            >>> frag_A = Fragment([read_A],umi_hamming_distance=1)
            >>> frag_B = Fragment([read_B],umi_hamming_distance=1)
            >>> frag_A.umi_eq(frag_B) # This returns True, the distance is 1, which is the (umi_hamming_distance)
            True
        """

        if self.umi == other.umi:
            return True
        if self.umi_hamming_distance == 0:
            return False
        else:
            if len(self.umi)!=len(other.umi):
                return False
            return hamming_distance(
                self.umi, other.umi) <= self.umi_hamming_distance

    def __eq__(self, other):  # other can also be a Molecule!
        # Make sure fragments map to the same strand, cheap comparisons
        """
        Check equivalence between two Fragments or Fragment and Molecule.


        Args:
            other (Fragment or Molecule) : object to compare against

        Returns
            is_eq (bool) : True when the other object is (likely) derived from the same molecule

        Example:
            >>> from singlecellmultiomics.molecule import Molecule
            >>> from singlecellmultiomics.fragment import Fragment
            >>> import pysam
            >>> # Create sam file to write some reads to:
            >>> test_sam = pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000])
            >>> read_A = pysam.AlignedSegment(test_sam.header)
            >>> read_A.set_tag('SM','CELL_1') # The sample to which the sample belongs is extracted from the SM tag
            >>> read_A.set_tag('RX','CAT') # The UMI is extracted from the RX tag
            >>> # By default the molecule assignment is done based on the mapping location of read 1:
            >>> read_A.reference_name = 'chr1'
            >>> read_A.reference_start = 100
            >>> read_A.query_sequence = 'ATCGGG'
            >>> read_A.cigarstring = '6M'

            >>> read_B = pysam.AlignedSegment(test_sam.header)
            >>> read_B.set_tag('SM','CELL_1')
            >>> read_B.set_tag('RX','CAT')
            >>> read_B.reference_start = 100
            >>> read_B.query_sequence = 'ATCGG'
            >>> read_B.cigarstring = '5M'

            >>> frag_A = Fragment([read_A],umi_hamming_distance=0)
            >>> frag_A
                Fragment:
                    sample:CELL_1
                    umi:CAT
                    span:chr1 100-106
                    strand:+
                    has R1: yes
                    has R2: no
                    randomer trimmed: no
            >>> frag_B = Fragment([read_B],umi_hamming_distance=0)
            >>> frag_A == frag_B
            True
            # Fragment A and fragment B belong to the same molecule,
            # the UMI is identical, the starting position of R1 is identical and
            # the sample name matches

        When we move one of the reads, the Fragments are not equivalent any more

        Example:
            >>> read_B.reference_start = 150
            >>> frag_B = Fragment([read_B],umi_hamming_distance=0)
            >>> frag_A == frag_B
            False

        Except if the difference <= the assignment_radius

        Example:
            >>> read_B.reference_start = 150
            >>> read_A.reference_start = 100
            >>> frag_B = Fragment([read_B],assignment_radius=300)
            >>> frag_A = Fragment([read_A],assignment_radius=300)
            >>> frag_A == frag_B
            True

         When the UMI's are too far apart, the eq function returns `False`

        Example:
            >>> read_B.reference_start = 100
            >>> read_A.reference_start = 100
            >>> read_A.set_tag('RX','GGG')
            >>> frag_B = Fragment([read_B])
            >>> frag_A = Fragment([read_A])
            >>> frag_A == frag_B
            False

        When the sample of the Fragments are not identical, the eq function returns `False`

        Example:
            >>> read_B.reference_start = 100
            >>> read_A.reference_start = 100
            >>> read_A.set_tag('RX','AAA')
            >>> read_B.set_tag('RX','AAA')
            >>> read_B.set_tag('SM', 'CELL_2' )
            >>> frag_B = Fragment([read_B])
            >>> frag_A = Fragment([read_A])
            >>> frag_A == frag_B
            False

        """
        if self.sample != other.sample:
            return False

        if self.strand != other.strand:
            return False

        if not self.has_valid_span() or not other.has_valid_span():
            return False

        if min(abs(self.span[1] -
                   other.span[1]), abs(self.span[2] -
                                       other.span[2])) > self.assignment_radius:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)

    def __getitem__(self, index):
        """
        Get a read from the fragment

        Args:
            index( int ):
                0 : Read 1
                1: Read 2
                ..
        """
        return self.reads[index]

    def is_valid(self):
        if self.qcfail:
            return False
        return self.has_valid_span()

    def get_strand_repr(self):
        s = self.get_strand()
        if s is None:
            return '?'
        if s:
            return '-'
        else:
            return '+'

    def __repr__(self):
        rep =  f"""Fragment:
        sample:{self.get_sample()}
        umi:{self.get_umi()}
        span:{('%s %s-%s' % self.get_span() if self.has_valid_span() else 'no span')}
        strand:{self.get_strand_repr()}
        has R1: {"yes" if self.has_R1() else "no"}
        has R2: {"yes" if self.has_R2() else "no"}
        randomer trimmed: {"yes" if self.unsafe_trimmed else "no"}
        """

        try:
            rep += '\n\t'.join([f'{key}:{str(value)}' for key, value in self.meta.items()])
        except Exception as e:
            print(self.meta)
            raise
        return rep

    def get_html(
            self,
            chromosome=None,
            span_start=None,
            span_end=None,
            show_read1=None,
            show_read2=None):
        """
        Get HTML representation of the fragment

        Args:
            chromosome( str ):
                chromosome to view

            span_start( int ):
                first base to show

            span_end( int ):
                last base to show

            show_read1(bool):
                show read1

            show_read2(bool):
                show read2

        Returns:
            html(string) : html representation of the fragment

        """

        if chromosome is None and span_start is None and span_end is None:
            chromosome, span_start, span_end = self.get_span()

        span_len = abs(span_end - span_start)
        visualized_frag = ['.'] * span_len

        if show_read1 is None and show_read2 is None:
            reads = self
        else:
            reads = []
            if show_read1:
                reads.append(self.get_R1())
            if show_read2:
                reads.append(self.get_R2())

        for read in reads:
            if read is None:
                continue
            # for cycle, query_pos, ref_pos, ref_base in
            # pysamiterators.iterators.ReadCycleIterator(read,with_seq=True):
            for query_pos, ref_pos, ref_base in read.get_aligned_pairs(
                    with_seq=True, matches_only=True):
                if ref_pos is None or ref_pos < span_start or ref_pos > span_end:
                    continue

                ref_base = ref_base.upper()
                if query_pos is None:
                    query_base = '-'
                else:
                    query_base = read.seq[query_pos]
                if ref_pos is not None:
                    if query_base != ref_base:
                        v = style_str(query_base, color='red', weight=500)
                    else:
                        v = query_base
                    try:
                        visualized_frag[ref_pos - span_start] = v
                    except IndexError:
                        pass  # the alignment is outside the requested view
        return ''.join(visualized_frag)

    def set_rejection_reason(self, reason, set_qcfail=False):
        self.set_meta('RR', reason, as_set=True)
        if set_qcfail:
            for read in self:
                if read is not None:
                    read.is_qcfail = True

    def set_recognized_sequence(self, seq):
        self.set_meta('RZ', seq)



class SingleEndTranscriptFragment(Fragment, FeatureAnnotatedObject):
    def __init__(self, reads, features, assignment_radius=0, stranded=None, capture_locations=False, auto_set_intron_exon_features=True,**kwargs):

        Fragment.__init__(self, reads, assignment_radius=assignment_radius, **kwargs)
        FeatureAnnotatedObject.__init__(self,
                                        features=features,
                                        stranded=stranded,
                                        capture_locations=capture_locations,
                                        auto_set_intron_exon_features=auto_set_intron_exon_features )

        self.match_hash = (self.sample,
                           self.span[0],
                           self.strand)


    def write_tags(self):
        # First decide on what values to write
        FeatureAnnotatedObject.write_tags(self)
        # Additional fragment data
        Fragment.write_tags(self)

    def is_valid(self):
        if len(self.genes)==0:
            return False
        return True

    def annotate(self):
        for read in self:
            if read is None:
                continue

            read_strand = read.is_reverse
            if self.stranded is not None and self.stranded==False:
                read_strand=not read_strand


            strand = ('-' if read_strand else '+')
            read.set_tag('mr',strand)
            for start, end in read.get_blocks():

                for hit in self.features.findFeaturesBetween(
                        chromosome=read.reference_name,
                        sampleStart=start,
                        sampleEnd=end,
                        strand=(None if self.stranded is None else strand)) :

                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    self.hits[hit_ids].add(
                        (read.reference_name, (hit_start, hit_end)))

                    if self.capture_locations:
                        if not hit_id in self.feature_locations:
                            self.feature_locations[hit_id] = []
                        self.feature_locations[hit_id].append( (hit_start, hit_end, hit_strand))


    def get_site_location(self):
        return self.span[0], self.start

    def __eq__(self, other):
        # Check if the other fragment/molecule is aligned to the same gene
        if self.match_hash != other.match_hash:
            return False

        if len(self.genes.intersection(other.genes)):
            return self.umi_eq(other)

        return False




class FeatureCountsSingleEndFragment(Fragment):
    """ Class for fragments annotated with featureCounts

    Extracts annotated gene from the XT tag.
    Deduplicates using XT tag and UMI
    Reads without XT tag are flagged as invalid

    """

    def __init__(self,
                 reads,
                 R1_primer_length=4,
                 R2_primer_length=6,
                 assignment_radius=100_000,
                 umi_hamming_distance=1,
                 invert_strand=False,
                 **kwargs
                 ):

        Fragment.__init__(self,
                          reads,
                          assignment_radius=assignment_radius,
                          R1_primer_length=R1_primer_length,
                          R2_primer_length=R2_primer_length,
                          umi_hamming_distance=umi_hamming_distance,**kwargs)

        self.strand = None
        self.gene = None
        self.identify_gene()

        if self.is_valid():
            self.match_hash = (self.strand, self.gene, self.sample)
        else:
            self.match_hash = None

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash != other.match_hash:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)

    def is_valid(self):
        return self.gene is not None

    def identify_gene(self):
        self.valid = False
        for read in self.reads:
            if read is not None and read.has_tag('XT'):
                self.gene = read.get_tag('XT')
                self.valid = True
                return self.gene



class FeatureCountsFullLengthFragment(FeatureCountsSingleEndFragment):
    """ Class for fragments annotated with featureCounts, with multiple reads covering a gene

    Extracts annotated gene from the XT tag.
    Deduplicates using XT tag and UMI
    Reads without XT tag are flagged as invalid

    """

    def __init__(self,
                 reads,
                 R1_primer_length=4,
                 R2_primer_length=6,
                 assignment_radius=10_000,
                 umi_hamming_distance=1,
                 invert_strand=False,
                 **kwargs
                 ):

        FeatureCountsSingleEndFragment.__init__(self,
                          reads,
                          assignment_radius=assignment_radius,
                          R1_primer_length=R1_primer_length,
                          R2_primer_length=R2_primer_length,
                          umi_hamming_distance=umi_hamming_distance,
                          **kwargs

                          )

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash != other.match_hash:
            return False

        #Calculate distance
        if min(abs(self.span[1] -
                   other.span[1]), abs(self.span[2] -
                                       other.span[2])) > self.assignment_radius:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)


class FragmentWithoutPosition(Fragment):
    """ Fragment without a specific location on a contig"""
    def get_site_location(self):
        if self.has_valid_span():
            return self.span[0], 0
        else:
            return '*', 0

class FragmentStartPosition(Fragment):
    """ Fragment without a specific location on a contig"""
    def get_site_location(self):
        for read in self:
            if read is not None and not read.is_unmapped:
                return read.reference_name, read.reference_start
        return None



class FragmentWithoutUMI(Fragment):
    """
    Use this class when no UMI information is available
    """

    def __init__(self, reads, **kwargs):
        Fragment.__init__(self, reads, **kwargs)

    # remove the set_umi function
    def set_umi(self, **kwargs):
        pass

    # Replace the equality function
    def __eq__(self, other):  # other can also be a Molecule!
        # Make sure fragments map to the same strand, cheap comparisons
        if self.sample != other.sample:
            return False

        if self.strand != other.strand:
            return False

        if min(abs(self.span[1] -
                   other.span[1]), abs(self.span[2] -
                                       other.span[2])) > self.assignment_radius:
            return False

        # Sample matches and starting position is within the defined span
        # radius
        return True
