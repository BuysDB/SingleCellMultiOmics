from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.sequtils import hamming_distance


class NLAIIIFragment(Fragment):
    def __init__(self,
                 reads,
                 R1_primer_length=4,
                 R2_primer_length=6,
                 assignment_radius=1_000,
                 umi_hamming_distance=1,
                 invert_strand=False,
                 no_overhang =False, # CATG is present OUTSIDE the fragment
                 reference=None #Reference is required when no_overhang=True
                 ):
        self.invert_strand = invert_strand
        self.no_overhang = no_overhang
        self.reference = reference
        if self.no_overhang and reference  is None:
            raise ValueError('Supply a reference handle when no_overhang=True')

        Fragment.__init__(self,
                          reads,
                          assignment_radius=assignment_radius,
                          R1_primer_length=R1_primer_length,
                          R2_primer_length=R2_primer_length,
                          umi_hamming_distance=umi_hamming_distance)
        # set NLAIII cut site given reads
        self.strand = None
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

    def set_site(self, site_chrom, site_pos, site_strand=None):
        self.set_meta('DS', site_pos)
        if site_strand is not None:
            if self.invert_strand:
                self.set_meta('RS', not site_strand)
            else:
                self.set_meta('RS', site_strand)
        self.set_strand(site_strand)
        self.site_location = (site_chrom, site_pos)
        self.cut_site_strand = site_strand

    def identify_site(self):

        R1, R2 = self.reads
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
        if R1 is None or R1.is_unmapped:
            self.set_rejection_reason('unmapped R1')
            return None

        if self.no_overhang:
            forward_motif = self.reference.fetch(R1.reference_name,  R1.reference_start-4,  R1.reference_start)
            rev_motif = self.reference.fetch(R1.reference_name, R1.reference_end, R1.reference_end+4)
        else:
            forward_motif = R1.seq[:4]
            rev_motif = R1.seq[-4:]


        if forward_motif == 'CATG' and not R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_start)
            self.set_site(site_strand=1, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('CATG')
            return(rpos)
        elif rev_motif == 'CATG' and R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_end - 4)
            self.set_site(site_strand=0, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('CATG')
            return(rpos)

        # Sometimes the cycle is off, this is dangerous though because the cell barcode and UMI might be shifted too!
        elif forward_motif.startswith( 'ATG') and not R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_start - 1)
            self.set_site(site_strand=1, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('ATG')
            return(rpos)
        elif rev_motif.startswith( 'ATG') and R1.is_reverse:  # First base was trimmed or lost
            rpos = (R1.reference_name, R1.reference_end - 3)
            self.set_site(site_strand=0, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('CAT')
            return(rpos)

        else:
            if forward_motif == 'CATG' and R1.is_reverse:
                self.set_rejection_reason('found CATG R1 REV exp FWD')

            elif rev_motif == 'CATG' and not R1.is_reverse:
                self.set_rejection_reason('found CATG R1 FWD exp REV')
            else:
                self.set_rejection_reason('no CATG')
            return None

    def get_undigested_site_count(self, reference_handle):
        """
        Obtain the amount of undigested sites in the span of the fragment

        Parameters
        ----------
        reference : pysam.FastaFile or similiar

        Returns
        -------
        undigested_site_count : int
            amount of undigested cut sites in the mapping span of the fragment

        Raises:
        -------
        ValueError : when the span of the molecule is not properly defined
        """
        if any(e is None for e in self.span):
            raise ValueError('Span of the fragment is not well defined')

        total = reference_handle.fetch(*self.span).count('CATG')
        # ignore 1 CATG if it was Detected:
        if self.meta.get('RZ') == 'CATG':
            total -= 1
        return total

    def is_valid(self):
        return self.site_location is not None

    def get_site_location(self):
        return self.site_location

    def __repr__(self):
        return Fragment.__repr__(
            self) + f'\n\tRestriction site:{self.get_site_location()}'

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash != other.match_hash:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)
