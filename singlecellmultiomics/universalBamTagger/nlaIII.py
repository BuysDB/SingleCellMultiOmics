from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.tagtools import tagtools


class NlaIIIFlagger(DigestFlagger):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs)

    def addSite(self, reads, strand, restrictionChrom, restrictionPos):

        if not reads[0].has_tag(
                self.sampleTag) or not reads[0].has_tag(
                self.umiTag):
            return

        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(
            self.alleleTag) else reads[0].get_tag(
            self.alleleTag)
        siteInfo = tuple([x for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(
            sample, restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            if read is None:
                continue
            self.setSiteOversequencing(read, moleculeId)
            self.setSiteCoordinate(read, restrictionPos)
            self.setSource(read, 'NLA'), {}
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, '+' if strand ==
                           1 else ('-' if strand == 0 else '?'))

    def digest(self, reads):
        if len(reads) != 2:
            if len(reads) == 1:
                self.setRejectionReason(reads[0], 'unmapped mate')
            else:
                self.setRejectionReason(reads[0], 'nopair')
            return None  # Only made for mate pair
        R1, R2 = reads

        self.addAlleleInfo([read for read in reads if read is not None])

        """ Valid configs:
        CATG######## R1 ########## ^ ########## R2 ##########
        ############ R2 ########## ^ ########### R1 #####CATG  reverse case
        !BWA inverts the query sequence if it maps to the negative strand!

        or R2.is_unmapped:
            if R1.is_unmapped and R2.is_unmapped:
                self.setRejectionReason(R1, 'unmapped R1;R2')
            elif R1.is_unmapped:
                self.setRejectionReason(R1, 'unmapped R1')
                self.setRejectionReason(R2, 'unmapped R1')
            else:
                self.setRejectionReason(R1, 'unmapped R2')
                self.setRejectionReason(R2, 'unmapped R2')
            return(None)
        """

        # Obtain RT hexamer:
        if R2 is not None:
            hstart, hseq = tagtools.getRandomPrimerHash(
                R2, onStart=True, primerLength=6)
            self.setRandomPrimer(R1, R2, hstart, hseq)

        if R1 is None or R1.is_unmapped:
            self.setRejectionReason(R1, 'unmapped R1')
            self.setRejectionReason(R2, 'unmapped R1')
            return None

        if R1.seq[:4] == 'CATG' and not R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_start)
            self.addSite([R1, R2], strand=0,
                         restrictionChrom=rpos[0], restrictionPos=rpos[1])
            self.setRecognizedSequence(R1, 'CATG')
            self.setRecognizedSequence(R2, 'CATG')
            return(rpos)
        elif R1.seq[-4:] == 'CATG' and R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_end - 4)
            self.addSite([R1, R2], strand=1,
                         restrictionChrom=rpos[0], restrictionPos=rpos[1])
            self.setRecognizedSequence(R1, 'CATG')
            self.setRecognizedSequence(R2, 'CATG')
            return(rpos)

        # Sometimes the cycle is off
        elif R1.seq[:3] == 'ATG' and not R1.is_reverse:
            rpos = (R1.reference_name, R1.reference_start - 1)
            self.addSite([R1, R2], strand=0,
                         restrictionChrom=rpos[0], restrictionPos=rpos[1])
            self.setRecognizedSequence(R1, 'ATG')
            self.setRecognizedSequence(R2, 'ATG')
            return(rpos)
        elif R1.seq[-3:] == 'CAT' and R1.is_reverse:  # First base was trimmed or lost
            rpos = (R1.reference_name, R1.reference_end - 3)
            self.addSite([R1, R2], strand=1,
                         restrictionChrom=rpos[0], restrictionPos=rpos[1])
            self.setRecognizedSequence(R1, 'CAT')
            self.setRecognizedSequence(R2, 'CAT')
            return(rpos)

        else:
            if R1.seq[:4] == 'CATG' and R1.is_reverse:
                self.setRejectionReason(R1, 'found CATG R1 REV exp FWD')
                self.setRejectionReason(R2, 'found CATG R1 REV exp FWD')

            elif R1.seq[-4:] == 'CATG' and not R1.is_reverse:
                self.setRejectionReason(R1, 'found CATG R1 FWD exp REV')
                self.setRejectionReason(R2, 'found CATG R1 FWD exp REV')
            else:
                self.setRejectionReason(R1, 'no CATG')
                self.setRejectionReason(R2, 'no CATG')
            return None

        try:
            start, end = tagtools.getPairGenomicLocations(
                R1, R2, R1PrimerLength=4, R2PrimerLength=6)
            self.setFragmentSize(R1, end - start)
            self.setFragmentSize(R2, end - start)
            self.setFragmentTrust(R1, start, end)
            self.setFragmentTrust(R2, start, end)

        except Exception as e:
            self.setFragmentSize(R1, 'unknown')
            self.setFragmentSize(R2, 'unknown')

        """
        if R1.seq[:4]=='CATG' and R1.reference_start<=R2.reference_start: # Site on the start of R1, R2 should map behind
            self.addSite( [R1,R2],  strand=0, restrictionChrom=R1.reference_name, restrictionPos=R1.reference_start )
            return(( R1.reference_name, R1.reference_start))

        if R1.seq[-4:]=='CATG' and R1.reference_start>=R2.reference_start: # Site on the end of R1, R2 should map before
            self.addSite( [R1,R2],  strand=1, restrictionChrom=R1.reference_name, restrictionPos=R1.reference_end-4 )
            return( (R1.reference_name, R1.reference_end-4))
        """
