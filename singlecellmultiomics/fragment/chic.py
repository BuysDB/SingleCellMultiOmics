from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.sequtils import hamming_distance
class CHICFragment(Fragment):
    def __init__(self,
        reads,
        R1_primer_length=4,
        R2_primer_length=6,
        assignment_radius=1_000,
        umi_hamming_distance=1
        ):
        Fragment.__init__(self,
            reads,
            assignment_radius=assignment_radius,
            R1_primer_length=R1_primer_length,
            R2_primer_length=R2_primer_length,
            umi_hamming_distance=umi_hamming_distance )
        # set CHIC cut site given reads
        self.strand = None
        self.site_location = None
        self.cut_site_strand = None
        self.identify_site()
        if self.is_valid():
            self.match_hash = (self.strand, self.cut_site_strand,
             self.site_location[0],self.site_location[1], self.sample)
        else:
            self.match_hash = None

    def set_site(self, site_chrom, site_pos, site_strand=None ):
        self.set_meta('DS',site_pos)
        if site_strand is not None:
            self.set_meta('RS',site_strand)
        self.set_strand(site_strand)
        self.site_location = (site_chrom, site_pos)
        self.cut_site_strand = site_strand

    def identify_site(self):

        if R1 is None or R2 is None:
            return None
        self.addAlleleInfo([read for read in reads if read is not None])

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

        is_trimmed = (R1.has_tag('MX') and R1.get_tag('MX').startswith('scCHIC'))

        if R1.is_unmapped or R2.is_unmapped:
            return(None)
        try:
            start, end = tagtools.getPairGenomicLocations(R1,R2, R1PrimerLength=1 - int(is_trimmed), R2PrimerLength=6)
            self.setFragmentSize(R1, end-start)
            self.setFragmentSize(R2, end-start)
            self.setFragmentTrust(R1, start, end)
            self.setFragmentTrust(R2, start, end)

        except Exception as e:
            self.setFragmentSize(R1, 'unknown')
            self.setFragmentSize(R2, 'unknown')


        #if R1.seq[0]=='T': # Site on the start of R1, R2 should map behind
        if is_trimmed:
            # The first base of the read has been taken off and the lh tag is already set, this can be copied to RZ

            self.set_site(
                site_strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                site_chrom=R1.reference_name,
                site_pos=(R1.reference_end if R1.is_reverse else R1.reference_start),
                is_trimmed=True

                  )
                return(site_pos)

        else:

            self.set_site(
                site_strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                site_chrom=R1.reference_name,
                site_pos=(R1.reference_end-1 if R1.is_reverse else R1.reference_start+1),
                is_trimmed=False)
                return(site_pos)

    def is_valid(self):
        return self.site_location is not None

    def get_site_location(self):
        return self.site_location

    def __repr__(self):
        return Fragment.__repr__(self)+f'\n\tRestriction site:{self.get_site_location()}'

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash!=other.match_hash:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance calculation
        return self.umi_eq(other)
