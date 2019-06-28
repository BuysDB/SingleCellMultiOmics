import collections

class RangeCache():

    def __init__(self, maxRange=1200):
        self.d = dict() #Chrom -> pos -> value
        self.current = 0
        self.maxRange = maxRange
        self.freedMemory = 0

    # Obtain record within a radius (center, Radius) , returns first (closest) instance it finds
    def getWithinRange(self, chrom, center, radius):
        if not chrom in self.d:
            return None
        for i in range(0, radius+1):
            r = self.get(chrom, i+center)
            if r is not None:
                return r
            r = self.get(chrom, center-i)
            if r is not None:
                return r

    def get(self, chrom, pos ):
        if not chrom in self.d:
            return None

        if not pos in self.d[chrom]:
            return None
        return self.d[chrom][pos]

    def set(self, chrom, pos, data):
        self.purge(chrom, pos)
        if not chrom in self.d:
            self.d[chrom] = {}
        self.d[chrom][pos] = data

    def purge( self, chrom, pos ):
        if chrom in self.d:
            drop = []
            for position in self.d[chrom]:
                if abs(position-pos)>self.maxRange:
                    drop.append(position)
            for d in drop:
                if d in self.d[chrom]:
                    del self.d[chrom][d]
                    self.freedMemory+=1
        # Remove all data on other chromosomes:
        drop = [x for x in self.d if x!=chrom]
        for d in drop:
            del self.d[d]
            self.freedMemory+=1

    def __len__(self):
        return sum( [len(self.d[chrom]) for chrom in self.d])

# A digest function should fullfill the following conditions:
# Accept : a read pair of R1, and R2 optional
# a reference genome handle (optional)
# an allele tool handle (optional)

class DigestFlagger():

    def __init__(self, reference=None, alleleResolver=None, moleculeRadius=0, verbose=False, **kwargs):
        self.reference = reference
        self.alleleResolver = alleleResolver
        self.verbose = verbose
        self.siteCoordinateTag = 'DS'
        self.oversequencedMoleculeTag = 'RC'
        self.recognitionSequenceTag = 'RZ'
        self.sourceTypeTag = 'DT' # DataType
        self.alleleTag = 'DA'
        self.alleleSetTag = 'AA'
        self.strandTag = 'RS'
        self.umiTag = 'RX'
        self.sampleTag = 'SM'
        self.ligationTag = 'LI'
        self.rejectionTag = 'RR'
        self.randomPrimerPosTag = 'rP'
        self.randomPrimerPosSeq = 'rS'

        self.fragmentSizeTag = 'fS'
        self.fragmentEndTag = 'fe'
        self.fragmentStartTag = 'fs'

        self.cacheWrites = 0

        self.cachePurgeEvery = 10_000 # Clean the cache every N writes

        self.moleculeRadius = moleculeRadius
        # This hold which sites we saw to count oversequencing
        self.observedSites = collections.defaultdict( RangeCache ) # sample -> (chrom,allele,strand,..) -> pos -> seen

    def __repr__(self):
        return f'Tagger with {self.cacheWrites} cache writes, total consumption: {self.getTotalConsumption()} '

    def getTotalConsumption(self):
        return sum( [len(self.observedSites[sample]) for sample in self.observedSites])


    def setRejectionReason(self, read, reason):
        self.appendTag(read, self.rejectionTag, reason)

    def setLigationSite(self, read, site):
        read.set_tag( self.ligationTag, site  )

    def setSiteCoordinate( self, read, coordinate ):
        self.appendTag(read, self.siteCoordinateTag, coordinate )

    def setSiteOversequencing( self, read, moleculeIndex=1 ): # 1 if first seen 2, second, -1 if None
        read.set_tag( self.oversequencedMoleculeTag, -1 if moleculeIndex is None else moleculeIndex  )
        # Decribe as string and set tag:

        if moleculeIndex>1:
            read.is_duplicate = True

    def setFragmentSize(self, read, size):
        read.set_tag( self.fragmentSizeTag,size)

    def setFragmentTrust(self, read, start, end):
        read.set_tag( self.fragmentStartTag,start)
        read.set_tag( self.fragmentEndTag,end)


    def setAllele( self, read, allele): # 1 if first seen 2, second, -1 if None
        read.set_tag( self.alleleTag, allele)

    def setRecognizedSequence( self, read, sequence ):
        self.appendTag(read, self.recognitionSequenceTag, sequence )

    def setSource( self, read, source ):
        self.appendTag(read, self.sourceTypeTag, source )

    def setStrand( self, read, strand ):
        read.set_tag( self.strandTag, strand )


    def setRandomPrimer(self, R1,R2, hstart, hseq ):
        if R1 is not None:
            R1.set_tag(self.randomPrimerPosSeq, hseq)
            R1.set_tag(self.randomPrimerPosTag, hstart)
        if R2 is not None:
            R2.set_tag(self.randomPrimerPosSeq, hseq)
            R2.set_tag(self.randomPrimerPosTag, hstart)

    def appendTag(self, read, tag, value):
        if read is None:
            return

        if not read.has_tag(tag):
            read.set_tag(tag, value)
        else:
            read.set_tag(tag, f'{read.get_tag(tag)},{value}' )

    def addAlleleInfo(self, reads):
        allele = None if self.alleleResolver is None else self.alleleResolver.getAllele(reads)
        if self.verbose:
            print(allele,reads)
        if allele is not None and len(allele)>0:
            allele = ','.join(sorted(list(allele)))

            for read in reads:
                if read is not None:
                    self.setAllele(read,allele)
        return allele

    # siteInfo describes the site as tuple: (allele, recognizedStrand, umiSequence, ..)
    # the site info can have arbitrary length, with as only requirement that it should be hashable
    def increaseAndRecordOversequencing( self, sample, chrom, pos, siteInfo=() ):
        cutSite = (chrom, pos)

        self.cacheWrites += 1
        if self.cacheWrites>self.cachePurgeEvery:
            for sample in self.observedSites:
                self.observedSites[sample].purge(*cutSite)
            self.cacheWrites=0

        if self.moleculeRadius!=0:
            current = self.observedSites[sample].getWithinRange(*cutSite, radius=self.moleculeRadius)
        else:
            current = self.observedSites[sample].get( *cutSite )

        # There are no matching molecules nearby:
        if current is None:
            current = collections.Counter()
            self.observedSites[sample].set( *cutSite, current )

        UMI = siteInfo
        current[UMI]+=1
        #print(current)
        return  current[UMI]

    ### ADD THIS FUNCTION:
    #def digest( self, reads )
    #
