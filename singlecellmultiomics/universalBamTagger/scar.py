from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.tagtools import tagtools
from singlecellmultiomics.utils.sequtils import phred_to_prob
import numpy as np


class ScarFlagger(DigestFlagger):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs)

    def addSite(self, reads, scarChromosome, scarPrimerStart):

        R1_primer_length = 20
        R2_primer_length = 18

        sample = reads[0].get_tag(self.sampleTag)
        allele = None if not reads[0].has_tag(
            self.alleleTag) else reads[0].get_tag(
            self.alleleTag)

        R1, R2 = reads
        if R1 is None or R2 is None:
            for read in reads:
                if read is not None:
                    self.setRejectionReason(read, 'NotPaired')

            return None
        # find all deletions:
        scarDescription = set()

        qualities = []
        for read in reads:

            firstCigarOperation, firstCigarLen = read.cigartuples[0]
            insertPos = 0
            lastNonNoneRefPos = read.reference_start if firstCigarOperation != 4 else read.reference_start - firstCigarLen

            expandedCigar = []
            for cigarOperation, cigarLen in read.cigartuples:
                expandedCigar += [cigarOperation] * cigarLen

            for queryPos, referencePos in read.get_aligned_pairs(
                    matches_only=False):
                if queryPos is not None:
                    qualities.append(phred_to_prob(read.qual[queryPos]))

                if queryPos is None and referencePos is None:
                    continue

                if referencePos is not None:
                    lastNonNoneRefPos = referencePos
                    insertPos = 0
                    if queryPos is not None:  # If both the reference and query match, there is not scar information
                        continue

                if queryPos is None:
                    scarDescription.add(f'{referencePos}.D')
                elif referencePos is None:  # insert or clip:
                    operation = expandedCigar[queryPos]
                    if operation == 1:
                        if lastNonNoneRefPos is None:
                            raise ValueError('Unsolvable :(')
                        queryBase = read.seq[queryPos]
                        scarDescription.add(
                            f'{queryBase}.{insertPos+lastNonNoneRefPos}.I')
                        insertPos += 1

        scarDescription = ','.join(sorted(list(scarDescription)))

        siteInfo = tuple(
            [x for x in [allele, scarDescription] if x is not None])

        moleculeId = self.increaseAndRecordOversequencing(
            sample, scarChromosome, scarPrimerStart, siteInfo=siteInfo)

        # Add average base calling quality excluding primers:
        meanQ = np.mean(qualities)
        for read in reads:

            read.set_tag('SQ', 1 - meanQ)

            self.setSiteOversequencing(read, moleculeId)
            if len(scarDescription) == 0:
                scarDescription = 'WT'
            read.set_tag('SD', scarDescription)
            self.setSiteCoordinate(read, R1.reference_start)
            self.setRecognizedSequence(read, 'SCAR')
            self.setSource(read, 'SCAR')
            if allele is not None:
                self.setAllele(read, allele)

    def digest(self, reads):
        if len(reads) != 2:
            return None  # Only made for mate pair
        R1, R2 = reads

        self.addAlleleInfo(reads)
        if R1 is None or R1.is_unmapped or R2 is None or R2.is_unmapped:
            return(None)

        self.addSite(reads, scarChromosome=R1.reference_name,
                     scarPrimerStart=R1.reference_start)
