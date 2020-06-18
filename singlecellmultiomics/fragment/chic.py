from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.sequtils import hamming_distance


class CHICFragment(Fragment):
    def __init__(self,
                 reads,
                 R1_primer_length=4,
                 R2_primer_length=6,
                 assignment_radius=1_000,
                 umi_hamming_distance=1,
                 invert_strand=False,
                 no_umi_cigar_processing=False,
                 **kwargs
                 ):
        self.invert_strand = invert_strand
        Fragment.__init__(self,
                          reads,
                          assignment_radius=assignment_radius,
                          R1_primer_length=R1_primer_length,
                          R2_primer_length=R2_primer_length,
                          umi_hamming_distance=umi_hamming_distance,
                          max_NUC_stretch = 18,
                          **kwargs

                )

        # set CHIC cut site given reads
        self.no_umi_cigar_processing = no_umi_cigar_processing
        self.strand = None
        self.ligation_motif = None
        self.site_location = None
        self.cut_site_strand = None
        self.identify_site()
        if self.is_valid():
            self.match_hash = (
                self.strand,
                self.cut_site_strand,
                self.site_location[0],
                self.site_location[1],
                self.sample)
        else:
            self.match_hash = None

    def set_site(
            self,
            site_chrom,
            site_pos,
            site_strand=None,
            is_trimmed=False):
        self.set_meta('DS', site_pos)
        if site_strand is not None:
            self.set_meta('RS', site_strand)
        self.set_strand(site_strand)
        self.site_location = (site_chrom, site_pos)
        self.cut_site_strand = site_strand

    def identify_site(self):

        R1 = self.get_R1()

        if R1 is None:
            self.set_rejection_reason("R1_undefined")
            return None

        if R1.has_tag('lh'):
            self.ligation_motif = R1.get_tag('lh')

        """ Valid configs:
        Observed read pair:
        T######## R1 ########## ^ ########## R2 ##########
        Real molecule:
        A###########################################################

        Real molecule: and adapter:
        ADAPTER: 5' NNNNT 3'
                        A{A:80%,N:20%}NNN [CHIC MOLECULE]
                                      ^ real cut site
        """

        is_trimmed = (R1.has_tag('MX') and R1.get_tag(
            'MX').startswith('scCHIC'))

        if R1.is_unmapped:
            self.set_rejection_reason("R1_unmapped")
            return(None)

        # Identify the start coordinate of Read 1 by reading the amount of softclips on the start of the read
        r1_start =(R1.reference_end if R1.is_reverse else R1.reference_start)
        if not self.no_umi_cigar_processing:
            if R1.is_reverse:
                if R1.cigartuples[-1][0]==4: # softclipped at end
                    r1_start+=R1.cigartuples[-1][1]
            else:
                if R1.cigartuples[0][0]==4: # softclipped at start
                    r1_start-=R1.cigartuples[0][1]

        if is_trimmed:
            # The first base of the read has been taken off and the lh tag is
            # already set, this can be copied to RZ

            self.set_site(
                site_strand=R1.is_reverse if not self.invert_strand else not R1.is_reverse,
                # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                site_chrom=R1.reference_name,
                site_pos=r1_start,
                is_trimmed=True
            )

        else:

            self.set_site(
                site_strand=R1.is_reverse if not self.invert_strand else not R1.is_reverse,
                # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                site_chrom=R1.reference_name,
                site_pos=(r1_start - 1 if R1.is_reverse else r1_start + 1),
                is_trimmed=False)

    def is_valid(self):
        if self.qcfail:
            return False

        if self.max_fragment_size is not None:
            try:
                size = self.get_fragment_size()
                if size>self.max_fragment_size:
                    return False
            except Exception:
                pass
        return self.site_location is not None

    def get_site_location(self):
        return self.site_location

    def __repr__(self):
        return Fragment.__repr__(
            self) + f'\n\tMNase cut site:{self.get_site_location()}'

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash != other.match_hash:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)
