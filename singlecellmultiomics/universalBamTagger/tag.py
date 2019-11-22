from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.tagtools import tagtools


class TagFlagger(DigestFlagger):

    def __init__(self, tag=None, **kwargs):
        DigestFlagger.__init__(self, **kwargs)
        if tag is None or (tag != 'chrom' and len(tag) != 2):
            raise ValueError(f'Invalid tag:{tag}')
        self.tag = tag

    def addSite(self, reads, strand, siteDef):
        sample = reads[0].get_tag(self.sampleTag)
        if reads[0].has_tag(self.umiTag):
            umi = reads[0].get_tag(self.umiTag)
        else:
            umi = 'N'
        allele = None if not reads[0].has_tag(
            self.alleleTag) else reads[0].get_tag(
            self.alleleTag)
        siteInfo = tuple(
            [x for x in [strand, allele, umi, siteDef] if x is not None])

        moleculeId = self.increaseAndRecordOversequencing(
            sample,
            reads[0].reference_name,
            0, siteInfo=siteInfo)
        # print(siteDef,siteInfo,moleculeId)

        for read in reads:
            self.setSiteOversequencing(read, moleculeId)
            self.setSiteCoordinate(read, 0)
            self.setRecognizedSequence(read, f'TAG:{self.tag}')
            self.setSource(read, 'TAG')
            if allele is not None:
                self.setAllele(read, allele)

            # Note: strand 1 is +
            self.setStrand(read, '+' if strand ==
                           1 else ('-' if strand == 0 else '?'))

    def digest(self, reads):
        self.addAlleleInfo([read for read in reads if read is not None])
        R1 = reads[0]
        if len(reads) == 2:
            R2 = reads[1]
        else:
            R2 = None
        for read in reads:
            if self.tag == 'chrom':
                if R1.reference_name is not None:
                    siteDef = str(R1.reference_name)
                if R2 is not None and R2.reference_name is not None:
                    siteDef = str(R2.reference_name)
            else:
                if R1.has_tag(self.tag):
                    siteDef = R1.get_tag(self.tag)
                elif R2 is not None and R2.has_tag(self.tag):
                    siteDef = R2.get_tag(self.tag)
                else:
                    return None

            if reads[0].is_read1:
                strand = reads[0].is_reverse
            else:
                strand = not reads[0].is_reverse
            self.addSite(reads,
                         strand=int(strand),
                         siteDef=siteDef
                         )
