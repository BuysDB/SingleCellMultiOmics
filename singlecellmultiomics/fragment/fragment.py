from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
from singlecellmultiomics.utils import style_str
complement = str.maketrans('ATCGN', 'TAGCN')

class Fragment():
    def __init__(self, reads, assignment_radius=3, umi_hamming_distance=1,R1_primer_length=0,R2_primer_length=6,
    tag_definitions=None ):
        if tag_definitions is None:
            tag_definitions = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TagDefinitions
        self.tag_definitions = tag_definitions
        self.assignment_radius = assignment_radius
        self.umi_hamming_distance = umi_hamming_distance
        self.reads = reads
        self.strand = None
        self.meta = {} # dictionary of meta data
        self.is_mapped = None
        ### Check for multimapping
        self.is_multimapped = True
        self.mapping_quality = 0
        self.match_hash=None

        # Span:\
        self.span=[None,None,None]

        self.umi=None

        # Force R1=read1 R2=read2:
        for i, read in enumerate(self.reads):

            if read is None:
                continue

            if not read.is_unmapped:
                self.is_mapped = True
            if read.mapping_quality!=0:
                self.is_multimapped = False
            self.mapping_quality = max(self.mapping_quality, read.mapping_quality)
            if i==0:
                if read.is_read2:
                    raise ValueError('Supply first R1 then R2')
                read.is_read1 = True
                read.is_read2 = False
            elif i==1:
                read.is_read1 = False
                read.is_read2 = True

        self.R1_primer_length = R1_primer_length
        self.R2_primer_length = R2_primer_length
        self.set_sample()
        self.set_strand(self.identify_strand())
        self.update_span()
        self.update_umi()

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
        Returns None,None when the random primer cannot be described
        Returns
        -------
        reference_start : int or None
        sequence : str or None
        """
        R2 = self.get_R2()
        if R2 is None or R2.query_sequence is None:
            return None,None
        # The read was not mapped
        if R2.is_unmapped:
            # Guess the orientation does not matter
            return None, R2.query_sequence[:self.R2_primer_length]

        if R2.is_reverse:
            global complement
            return(R2.reference_end, R2.query_sequence[-self.R2_primer_length:][::-1].translate(complement))
        else:
            return(R2.reference_start, R2.query_sequence[:self.R2_primer_length])
        raise ValueError()


    def set_meta(self, key, value):
        self.meta[key] = value

    def get_R1(self):
        if len(self.reads)==0:
            raise IndexError('The fragment has no associated reads')
        return self.reads[0]
    def get_R2(self):
        if len(self.reads)<2:
            raise IndexError('The fragment has no associated R2')
        return self.reads[1]

    def update_span(self):
        surfaceStart = None
        surfaceEnd = None
        contig = None
        for read in self.reads:
            if read is None:
                continue
            if contig is None and read.reference_name is not None:
                contig = read.reference_name
            if read.reference_start is not None and (surfaceStart is None or read.reference_start<surfaceStart):
                surfaceStart=read.reference_start
            if read.reference_end is not None and( surfaceEnd is None or read.reference_end>surfaceEnd):
                surfaceEnd=read.reference_end

        self.span = (contig,surfaceStart, surfaceEnd)

    def get_span(self):
        return self.span

    def has_valid_span(self):
        return not any([x is None for x in self.get_span()])

    def set_sample(self, sample=None):
        """Force sample name or obtain sample name from associated reads"""
        if sample is not None:
            self.sample = sample
        else:
            for read in self.reads:
                if read is not None and read.has_tag('SM'):
                    self.sample = read.get_tag('SM')
                    break


    def get_sample(self):
        return self.sample


    def set_strand(self, strand):
        self.strand = strand

    def get_strand(self):
        return self.strand

    def update_umi(self):
        for read in self.reads:
            if read is not None and read.has_tag('RX'):
                self.umi=read.get_tag('RX')

    def get_umi(self):
        return self.umi

    def __iter__(self):
        return iter(self.reads)

    def umi_eq(self, other):
        if self.umi==other.umi:
            return True
        if self.umi_hamming_distance==0:
            return False
        else:
            return hamming_distance(self.umi,other.umi)<=self.umi_hamming_distance

    def __eq__(self, other): # other can also be a Molecule!
        # Make sure fragments map to the same strand, cheap comparisons
        if self.sample!=other.sample:
            return False

        if self.strand!=other.strand:
            return False

        if min(  abs( self.span[1] - other.span[1] ),  abs( self.span[2] - other.span[2] ) ) >self.assignment_radius:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance calculation
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
        return f"""Fragment:
        sample:{self.get_sample()}
        umi:{self.get_umi()}
        span:{('%s %s-%s' % self.get_span())}
        strand:{self.get_strand_repr()}
        """

    def get_html(self, chromosome=None, span_start=None, span_end=None):
        """
        Get HTML representation of the fragment

        Args:
            chromosome( str ):
                chromosome to view

            span_start( int ):
                first base to show

            span_end( int ):
                last base to show

        Returns:
            html(string) : html representation of the fragment

        """

        if chromosome is None and span_start is None and span_end is None:
            chromosome, span_start, span_end = self.get_span()

        span_len = span_end - span_start
        visualized_frag = ['.']  * span_len
        for read in self:
            if read is None:
                continue
            for cycle, query_pos, ref_pos, ref_base in pysamiterators.iterators.ReadCycleIterator(read,with_seq=True):
                if ref_pos is None:
                    continue
                ref_base = ref_base.upper()
                query_base = read.seq[query_pos]
                if ref_pos is not None:
                    if query_base!=ref_base:
                         v = style_str(query_base,color='red',weight=500)
                    else:
                        v= query_base

                    visualized_frag[ref_pos-span_start] = v
        return ''.join(visualized_frag)

    def set_rejection_reason(self, reason):
        self.set_meta('RR',reason)

    def set_recognized_sequence(self, seq):
        self.set_meta('RZ',seq)


class SingleEndTranscript(Fragment):
    def __init__(self,reads, R2_primer_length=6):
        Fragment.__init__(self, reads, assignment_radius=1_000, umi_hamming_distance=1,R2_primer_length=R2_primer_length )
