from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.sequtils import hamming_distance


class NlaIIIFragment(Fragment):
    def __init__(self,
                 reads,
                 R1_primer_length=4,
                 R2_primer_length=6,
                 assignment_radius=1_000,
                 umi_hamming_distance=1,
                 invert_strand=False,
                 check_motif=True,
                 no_overhang =False, # CATG is present OUTSIDE the fragment
                 cut_location_offset=-4,
                 reference=None, #Reference is required when no_overhang=True
                 allow_cycle_shift=False,
                 use_allele_tag = False, # Use existing DA tag for molecule deduplication
                 no_umi_cigar_processing=False, **kwargs
                 ):
        self.invert_strand = invert_strand
        self.no_overhang = no_overhang
        self.reference = reference
        self.allow_cycle_shift = allow_cycle_shift
        self.cut_location_offset = cut_location_offset
        self.no_umi_cigar_processing= no_umi_cigar_processing
        self.use_allele_tag = use_allele_tag
        self.check_motif=check_motif

        if self.no_overhang and reference  is None:
            raise ValueError('Supply a reference handle when no_overhang=True')
        if self.no_overhang and not self.check_motif:
            raise ValueError('no_overhang=True is not compatible with check_motif=False, as there is no way to deduplicate. Consider using method "chic"')

        Fragment.__init__(self,
                          reads,
                          assignment_radius=assignment_radius,
                          R1_primer_length=R1_primer_length,
                          R2_primer_length=R2_primer_length,
                          umi_hamming_distance=umi_hamming_distance, **kwargs)
        # set NLAIII cut site given reads
        self.strand = None
        self.site_location = None
        self.cut_site_strand = None
        if self.identify_site():

            self.found_valid_site = True
        else:
            self.found_valid_site = False

        if self.is_valid():

            if self.use_allele_tag:

                allele = 'n'
                for read in self.reads:
                    if read is not None and read.has_tag('DA'):
                        allele = read.get_tag('DA')

                self.match_hash = (
                    allele,
                    self.strand,
                    self.cut_site_strand,
                    self.site_location[0],
                    self.site_location[1],
                    self.sample)

            else:
                self.match_hash = (
                    self.strand,
                    self.cut_site_strand,
                    self.site_location[0],
                    self.site_location[1],
                    self.sample)
        else:
            self.match_hash = None

    def set_site(self, site_chrom, site_pos, site_strand=None, valid=True):
        if not valid :
            self.found_valid_site = False
        else:
            self.found_valid_site = True
            self.set_meta('DS', site_pos)

        if site_strand is not None:
            if self.invert_strand:
                self.set_meta('RS', not site_strand)
            else:
                self.set_meta('RS', site_strand)
        self.set_strand(site_strand)
        self.site_location = (site_chrom, site_pos)
        self.cut_site_strand = site_strand

    def get_safe_span(self, allow_unsafe=True):

    # For a mapped read pair it can be important to figure out which bases are actual genomic signal
    # Al sequence behind and aligned to the random primer(s) cannot be trusted and should be masked out
    # secondly all signal before the starting location of R1 cannot be trusted
    # This function returns a lower and higher bound of the locations within the fragment that can be trusted
    # ASCII art: (H is primer sequence)
    #           R1 H------------------------>
    #   <------------------HH R2

    # Output: (E is emitted)
    #           R1 HEEEEEEE----------------->
    #   <------------------HH R2


        if self.R1_primer_length==0 and self.R2_primer_length==0:
            starts = tuple( read.reference_start for read in self if read is not None and not read.is_unmapped )
            ends = tuple( read.reference_end for read in self if read is not None and not read.is_unmapped )
            return min( min(starts), min(ends) ), max( max(starts), max(ends) )

        R1, R2 = self.reads

        if (R1 is None or R1.is_unmapped) and allow_unsafe:
            return R2.reference_start,R2.reference_end

        if (R2 is None or R2.is_unmapped) and allow_unsafe:
            return R1.reference_start,R1.reference_end


        if  (R1 is None and not allow_unsafe) or R2 is None or R1.is_unmapped:
            # This is an annoying situation, we cannot determine what bases can be trusted
            raise ValueError('Genomic locations cannot be determined')

        if R2.is_unmapped:
            # This is an annoying situation, we cannot determine what bases can be trusted
            raise ValueError('Genomic locations cannot be determined')

        if R1.is_reverse==R2.is_reverse:
            raise ValueError('Fragment incorrectly mapped') # The fragment is not correctly mapped

        if not R1.is_reverse: # R1 is forward, R2 is reverse
            #           R1 H------------------------>
            #   <------------------HH R2
            start = R1.reference_start+ self.R1_primer_length
            end = R2.reference_end -  self.R2_primer_length
        else:
            #           R2 HH------------------------>
            #   <------------------HH R1
            start = R2.reference_start+ self.R2_primer_length
            end = R1.reference_end -  self.R1_primer_length

        if start>=end:
            raise ValueError('Fragment has no size')

        return start,end


    def identify_site(self):
        self.found_valid_site = False
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

            # scan 3 bp of sequence for CATG
            scan_extra_bp = 3
            site_coordinate = None
            if R1.is_reverse:
                motif = self.reference.fetch(R1.reference_name, R1.reference_end, R1.reference_end-self.cut_location_offset+scan_extra_bp)
                if 'CATG' in motif:
                    offset = motif.find('CATG')
                    site_coordinate = R1.reference_end + offset
            else:
                motif = self.reference.fetch(R1.reference_name,  R1.reference_start+self.cut_location_offset-scan_extra_bp,  R1.reference_start)
                if 'CATG' in motif:
                    offset = motif[::-1].find('GTAC')
                    site_coordinate = R1.reference_start - offset - 4

            if site_coordinate is None:
                self.set_rejection_reason('no_CATG_in_ref', set_qcfail=True)
                return None

            self.set_site(site_strand=R1.is_reverse, site_chrom=R1.reference_name, site_pos=site_coordinate)
            self.set_recognized_sequence(motif)
            return site_coordinate

        else:
            forward_motif = R1.seq[:4]
            rev_motif = R1.seq[-4:]


        r1_start =(R1.reference_end if R1.is_reverse else R1.reference_start)
        if not self.no_umi_cigar_processing:
            if R1.is_reverse:
                if R1.cigartuples[-1][0]==4: # softclipped at end
                    r1_start+=R1.cigartuples[-1][1]
            else:
                if R1.cigartuples[0][0]==4: # softclipped at start
                    r1_start-=R1.cigartuples[0][1]

        if (not self.check_motif or forward_motif == 'CATG') and not R1.is_reverse: # Not reverse = forward
            rpos = (R1.reference_name, r1_start)
            self.set_site(site_strand=R1.is_reverse, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence(forward_motif)
            return(rpos)
        elif (not self.check_motif or rev_motif == 'CATG') and R1.is_reverse:
            rpos = (R1.reference_name, r1_start - 4)
            self.set_site(site_strand=R1.is_reverse, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence(rev_motif)
            return(rpos)

        # Sometimes the cycle is off, this is dangerous though because the cell barcode and UMI might be shifted too!
        elif self.allow_cycle_shift and  forward_motif.startswith( 'ATG') and not R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_start - 1)
            self.set_site(site_strand=R1.is_reverse, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('ATG')
            return(rpos)
        elif self.allow_cycle_shift and rev_motif.startswith( 'ATG') and R1.is_reverse:  # First base was trimmed or lost
            rpos = (R1.reference_name, R1.reference_end - 3)
            self.set_site(site_strand=R1.is_reverse, site_chrom=rpos[0], site_pos=rpos[1])
            self.set_recognized_sequence('CAT')
            return(rpos)

        else:

            if forward_motif == 'CATG' and R1.is_reverse:
                self.set_rejection_reason('found CATG R1 REV exp FWD', set_qcfail=True)

            elif rev_motif == 'CATG' and not R1.is_reverse:
                self.set_rejection_reason('found CATG R1 FWD exp REV', set_qcfail=True)
            else:
                self.set_rejection_reason('no CATG', set_qcfail=True)

            # Every fragment needs to have a site. Otherwise it will get lost. Use R1 start location as anchor
            if R1.is_reverse:
                rpos = (R1.reference_name, r1_start - 4)
            else:
                rpos = (R1.reference_name, r1_start)

            self.set_site(site_strand=R1.is_reverse, site_chrom=rpos[0], site_pos=rpos[1], valid=False)
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
        if self.qcfail:
            return False

        try:
            if self.max_fragment_size is not None and self.get_fragment_size()>self.max_fragment_size:
                self.set_rejection_reason('FS', set_qcfail=True)
                return False
        except TypeError:
            pass

        return self.found_valid_site

    def get_site_location(self):
        if self.site_location is not None:
            return self.site_location
        else:
            # We need some kind of coordinate...
            for read in self:
                if read is not None and read.reference_name is not None and read.reference_start is not None:
                    return read.reference_name, read.reference_start

    def __repr__(self):
        site_loc = self.get_site_location()

        if site_loc is None or len(site_loc)==0:
            return Fragment.__repr__(
            self) + f'\n\tNo restriction site found'
        else:
            site_def_str = f'\n\tRestriction site:{site_loc[0]}:{site_loc[1]}'
            try:
                return Fragment.__repr__(
                    self) + site_def_str
            except Exception as e:
                print(self.get_site_location())
                raise

    def __eq__(self, other):
        # Make sure fragments map to the same strand, cheap comparisons
        if self.match_hash != other.match_hash:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance
        # calculation
        return self.umi_eq(other)
