from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.tagtools import tagtools
complement = str.maketrans('ATGC', 'TACG')
class ChicSeqFlagger( DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def addSite(self, reads, strand, restrictionChrom, restrictionPos,is_trimmed=False ):

        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(  sample,  restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            if read is None:
                continue
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, restrictionPos)

            if is_trimmed:
                self.setRecognizedSequence(read,read.get_tag('lh'))
            else:
                if reads[0].is_reverse: #
                    self.setRecognizedSequence(read, reads[0].seq[-2:][::-1].translate(complement)) # the first two base, this should be A{A:80%, N:20%}, we take the complement because the reads captures the complement strand
                else:
                    self.setRecognizedSequence(read, reads[0].seq[:2]) # the last two bases
            self.setSource(read, 'CHIC')
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, '+' if strand==1 else ('-' if strand==0 else '?')) #Note: strand 1 is +

    def digest(self, reads):
        if len(reads)!=2:
            return None # Only made for mate pair

        R1,R2 = reads
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

            self.addSite( [R1,R2],
                strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                restrictionChrom=R1.reference_name,
                restrictionPos=(R1.reference_end if R1.is_reverse else R1.reference_start),
                is_trimmed=True

                  )

        else:

            self.addSite( [R1,R2],
                strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                restrictionChrom=R1.reference_name,
                restrictionPos=(R1.reference_end-1 if R1.is_reverse else R1.reference_start+1),
                is_trimmed=False)
        return(( R1.reference_name, R1.reference_start))
