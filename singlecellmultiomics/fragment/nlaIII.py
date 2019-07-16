from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.sequtils import hamming_distance
class NLAIIIFragment(Fragment):
    def __init__(self,reads, R1_primer_length=4,R2_primer_length=6,assignment_radius=1_000, umi_hamming_distance=1  ):
        Fragment.__init__(self, reads, assignment_radius=assignment_radius,R1_primer_length=R1_primer_length,R2_primer_length=R2_primer_length, umi_hamming_distance=umi_hamming_distance )
        # set NLAIII cut site given reads
        self.strand = None
        self.site_location = None
        self.cut_site_strand = None
        self.identify_site()

    def set_site(self, site_chrom, site_pos, site_strand=None ):

        self.set_strand(site_strand)
        self.site_location = (site_chrom, site_pos)
        self.cut_site_strand = site_strand

    def identify_site(self):

        R1,R2 = self.reads

        """ Valid configs:
        CATG######## R1 ########## ^ ########## R2 ##########
        ############ R2 ########## ^ ########### R1 #####CATG  reverse case
        !BWA inverts the query sequence if it maps to the negative strand!

        or R2.is_unmapped:
            if R1.is_unmapped and R2.is_unmapped:
                self.set_rejection_reason(R1, 'unmapped R1;R2')
            elif R1.is_unmapped:
                self.set_rejection_reason(R1, 'unmapped R1')
                self.set_rejection_reason(R2, 'unmapped R1')
            else:
                self.set_rejection_reason(R1, 'unmapped R2')
                self.set_rejection_reason(R2, 'unmapped R2')
            return(None)
        """
        if R1 is None or R1.is_unmapped :
            self.set_rejection_reason('unmapped R1')
            return None

        if R1.seq[:4]=='CATG' and not R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_start)
            self.set_site(  site_strand=0, site_chrom=rpos[0], site_pos=rpos[1] )
            self.set_recognized_sequence( 'CATG')
            return(rpos)
        elif R1.seq[-4:]=='CATG' and R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_end-4)
            self.set_site(  site_strand=1, site_chrom=rpos[0], site_pos=rpos[1] )
            self.set_recognized_sequence( 'CATG')
            return(rpos)

        # Sometimes the cycle is off
        elif R1.seq[:3]=='ATG' and not R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_start-1)
            self.set_site(  site_strand=0, site_chrom=rpos[0], site_pos=rpos[1] )
            self.set_recognized_sequence( 'ATG')
            return(rpos)
        elif R1.seq[-3:]=='CAT' and R1.is_reverse: # First base was trimmed or lost
            rpos = ( R1.reference_name, R1.reference_end-3)
            self.set_site(   site_strand=1, site_chrom=rpos[0], site_pos=rpos[1] )
            self.set_recognized_sequence( 'CAT')
            return(rpos)

        else:
            if R1.seq[:4]=='CATG' and R1.is_reverse:
                self.set_rejection_reason('found CATG R1 REV exp FWD')

            elif R1.seq[-4:]=='CATG' and not R1.is_reverse:
                self.set_rejection_reason( 'found CATG R1 FWD exp REV')
            else:
                self.set_rejection_reason( 'no CATG')
            return None


    def is_valid(self):
        return self.site_location is not None

    def get_site_location(self):
        return self.site_location

    def __repr__(self):
        return Fragment.__repr__(self)+f'\n\tRestriction site:{self.get_site_location()}'

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.get_sample()!=other.get_sample():
            return False

        if self.get_strand()!=other.get_strand():
            return False

        if self.get_site_location()!=other.get_site_location():
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance calculation
        if self.umi_hamming_distance==0:
            return self.get_umi()==other.get_umi()
        else:
            return hamming_distance(self.get_umi(),other.get_umi())<=self.umi_hamming_distance
