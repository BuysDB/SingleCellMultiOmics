from singlecellmultiomics.utils.sequtils import hamming_distance
import pysamiterators.iterators

class Fragment():
    def __init__(self, reads, assignment_radius=3, umi_hamming_distance=1,R1_primer_length=0,R2_primer_length=6 ):
        self.assignment_radius = assignment_radius
        self.umi_hamming_distance = umi_hamming_distance
        self.reads = reads

        self.is_multimapped = True
        # Force R1=read1 R2=read2:
        for i, read in enumerate(self.reads):
            if read is None:
                continue

            if read.mapping_quality!=0:
                self.is_multimapped = False

            if i==0:
                read.is_read1 = True
                read.is_read2 = False
            elif i==1:
                read.is_read1 = False
                read.is_read2 = True

        self.R1_primer_length = R1_primer_length
        self.R2_primer_length = R2_primer_length
        self.set_sample()

    def get_R1(self):
        return self.reads[0]
    def get_R2(self):
        return self.reads[1]

    def get_span(self):
        surfaceStart = None
        surfaceEnd = None
        contig = None
        for read in self.reads:
            if read is None:
                continue
            if contig is None and read.reference_name is not None:
                contig = read.reference_name
            if surfaceStart is None or read.reference_start<surfaceStart:
                surfaceStart=read.reference_start
            if surfaceEnd is None or read.reference_end>surfaceEnd:
                surfaceEnd=read.reference_end

        return contig,surfaceStart,surfaceEnd

    def set_sample(self):
        for read in self.reads:
            if read is not None and read.has_tag('SM'):
                self.sample = read.get_tag('SM')
        return None

    def get_sample(self):
        return self.sample

    def get_strand(self):
        for read in self.reads:
            if read is not None and  read.has_tag('RS'):
                return read.get_tag('RS')
        return None

    def get_umi(self):
        for read in self.reads:
            if read is not None and read.has_tag('RX'):
                return read.get_tag('RX')
        return None

    def __iter__(self):
        return iter(self.reads)

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.get_sample()!=other.get_sample():
            return False

        if self.get_strand()!=other.get_strand():
            return False

        spanSelf =  self.get_span()
        spanOther =  other.get_span()

        if min(  abs( spanSelf[1] - spanOther[1] ),  abs( spanSelf[2] - spanOther[2] ) ) >self.assignment_radius:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance calculation
        if self.umi_hamming_distance==0:
            return self.get_umi()==other.get_umi()
        else:
            return hamming_distance(self.get_umi(),other.get_umi())<=self.umi_hamming_distance

    def __repr__(self):
        return f"""Fragment:
        sample:{self.get_sample()}
        umi:{self.get_umi()}
        span:{('%s %s-%s' % self.get_span())}
        strand:{self.get_strand()}
        """

class SingleEndTranscript(Fragment):
    def __init__(self,reads, random_primer_length=6):
        Fragment.__init__(self, reads, assignment_radius=1_000, umi_hamming_distance=1 )
