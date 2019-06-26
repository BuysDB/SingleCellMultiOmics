from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.tagtools import tagtools

class MSPJIFlagger(DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )
        if self.reference is None:
            raise ValueError('The MSPJI digest flagger requires a reference genome file (FASTA)')

    def addSite(self, reads, strand, context, restrictionChrom, restrictionPos, ligationSite ):
        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(  sample,  restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, restrictionPos)
            self.setRecognizedSequence(read, context)
            self.setSource(read, 'MSPJI')
            self.setLigationSite(read, ligationSite)
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, strand)

    def digest(self, reads):
        if len(reads)!=2:
            return None # Only made for mate pair
        R1,R2 = reads
        if R1 is None or R2 is None:
            return None
        window = 20
        allele = self.addAlleleInfo(reads)
        drop=False
        for operation, length in R1.cigartuples:
            if operation!=0:
                drop=True
                break

        methylatedSite = None

        # CNNR NNNN NNNN | NNNN
        # CNNY NNNN NNNN  NNNNN
        if not R1.is_reverse:
            adapter = R1.seq[:4]
            seq =  self.reference.fetch( R1.reference_name, R1.reference_start-window, R1.reference_start+window+1 )
            detectedForwardFlank = seq[window-13]=='C'
            detectedForwardFlankMotif = seq[window-13:window-11]
            detectedForwardFragmentMotif = seq[window+15:window+17].translate(complement)[::-1]
            detectedForwardFragment = detectedForwardFragmentMotif[0]=='C'

            if detectedForwardFlank and not detectedForwardFragment:
                patternType = f'FWD13_{detectedForwardFlankMotif}'
                methylatedSite = R1.reference_start-13
                strand = '+'
            elif detectedForwardFragment and not detectedForwardFlank:
                patternType = f'REV16_{detectedForwardFragmentMotif}'
                methylatedSite = R1.reference_start-16
                strand = '-'

            elif detectedForwardFlank and detectedForwardFragment:
                patternType = 'FWDB'
                methylatedSite = None
            else:
                patternType = 'FWDO'
                methylatedSite = None
        else:
            adapter = R1.seq[-4:].translate(complement)[::-1]
            seq =  self.reference.fetch( R1.reference_name, R1.reference_end-window-1, R1.reference_end+window )
            detectedReverseFlank = seq[window-16]=='C'
            detectedReverseFlankMotif = seq[window-16:window-14]
            detectedReverseFragmentMotif = seq[window+12:window+14].translate(complement)[::-1]
            detectedReverseFragment = detectedReverseFragmentMotif[0]=='C'
            if detectedReverseFlank and not detectedReverseFragment:
                patternType = f'FWD16_{detectedReverseFlankMotif}'
                methylatedSite = R1.reference_start-16
                strand = '+'
            elif detectedReverseFragment and not detectedReverseFlank:
                patternType = f'REV12_{detectedReverseFragmentMotif}'
                methylatedSite = R1.reference_start-12
                strand = '-'
            elif detectedReverseFragment and detectedReverseFlank:
                methylatedSite = None
                patternType = f'REVB'
            else:
                patternType = f'REVO'
                methylatedSite = None

        #Check if we already saw this site:
        if methylatedSite is None:
            return None


        self.addSite( [R1,R2],
            strand=strand,
            restrictionChrom=R1.reference_name,
            restrictionPos=methylatedSite,
            context=patternType.split('_')[-1],
            ligationSite=adapter
         )
