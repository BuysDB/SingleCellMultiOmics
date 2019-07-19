#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import singlecellmultiomics.fastqProcessing.fastqIterator as fastqIterator
import string


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

class SamTag:

    def __init__(self, tag, humanName, isPhred=False, doNotWrite=False):
        self.tag = tag
        self.humanName = humanName
        self.isPhred=isPhred
        self.doNotWrite = doNotWrite
        if len(tag)!=2:
            raise ValueError('Invalid tag length')

tags = [
    SamTag('SM', 'sampleName'),
    SamTag('LY', 'sequencingLibrary'),
    SamTag('EX', 'experiment'),
    SamTag('MX', 'demultiplexingStrategy'),
    SamTag('Fc', 'flowCellId'),
    SamTag('La', 'lane'),
    SamTag('Ti', 'tile'),
    SamTag('Is', 'instrument'),
    SamTag('BK', 'isBulk'),
    SamTag('XC', 'plateCoord'),
    SamTag('CX', 'clusterXpos'),
    SamTag('CY', 'clusterYpos'),
    SamTag('RP', 'readPairNumber', doNotWrite=True),
    SamTag('RN', 'runNumber'),
    SamTag('Fi', 'filterFlag'),
    SamTag('CN', 'controlNumber'),
    SamTag('aA', 'correctedIlluminaIndexSequence'),
    SamTag('aI', 'indexSequence'),
    SamTag('aa', 'rawIlluminaIndexSequence'),
    SamTag('ah', 'hammingDistanceToIndex'),
    SamTag('BC', 'barcode'),
    SamTag('CB', 'cellBarcode'), # Same as bc
    SamTag('QX', 'barcodeQual', isPhred=True),
    SamTag('bc', 'rawBarcode'),
    SamTag('hd', 'hammingDistanceRawBcAssignedBc'),
    SamTag('BI', 'barcodeIndex'),
    SamTag('QT', 'sampleBarcodeQuality', isPhred=True),
    SamTag('RX', 'umi'),
    SamTag('uL', 'umiLength'),
    SamTag('RQ', 'umiQual', isPhred=True),
    SamTag('BX', 'umiCorrected'),
    SamTag('BZ', 'umiCorrectedQuality',isPhred=True),
    SamTag('MI', 'molecularIdentifier'),
    SamTag('QM', 'molecularIdentifierQuality',isPhred=True),
    SamTag('DS', 'siteCoordinate'),
    SamTag('DA', 'allele'),
    SamTag('RZ', 'recognizedSequence'),
    SamTag('RS', 'recognizedStrand'),
    SamTag('LI', 'ligationMotif'),
    SamTag('RC', 'moleculeOverseqCountIndex'),
    SamTag('RR', 'rejectionReason'),
    SamTag('DT', 'sourceType'),
    SamTag('EX', 'exons'),
    SamTag('IN', 'introns'),
    SamTag('GN', 'genes'),
    SamTag('nM', 'editDistanceToReference'),
    SamTag('NM', 'editDistanceToReference'),
    SamTag('AS', 'alignmentScore'),
    SamTag('NH', 'occurencesOfFragment'),
    SamTag('HI', 'queryHitIndex'),
    SamTag('XS', 'assignmentStatus'),
    SamTag('SQ', 'meanBaseQuality'),
    SamTag('SD', 'scarDescription'),
    SamTag('RG', 'readGroup'),
    SamTag('XA', 'alternativeAlignmentHits'),
    SamTag('MD', 'alignment'),


    SamTag('fS', 'fragmentSizeTag'),
    SamTag('fe', 'fragmentEndTag'),
    SamTag('fs', 'fragmentStartTag'),

    SamTag('lh', 'ligationOverhangSequence'),
    SamTag('lq', 'ligationOverhangQuality',isPhred = True),

    SamTag('a1', 'Read1Clipped3primeBasesKmer'),
    SamTag('A2', 'Read2Clipped3primeBasesKmer'),
    SamTag('aQ', 'Read1Clipped3primeBasesKmerQuality',isPhred = True),
    SamTag('AQ', 'Read2Clipped3primeBasesKmerQuality',isPhred = True),

    SamTag('e1', 'Read1Clipped5primeBasesKmer'),
    SamTag('eQ', 'Read2Clipped5primeBasesKmer'),
    SamTag('eB', 'Read1Clipped5primeCycle'),
    SamTag('EB', 'Read2Clipped5primeCycle'),
    SamTag('E2', 'Read1Clipped5primeBasesKmerQuality',isPhred = True),
    SamTag('EQ', 'Read2Clipped5primeBasesKmerQuality',isPhred = True),

    SamTag('ES', 'EnzymeSequence'),  # Christoph
    SamTag('eq', 'EnzymeSequenceQualities', isPhred = True),
    SamTag('IS', 'ispcrSequence'),
    SamTag('is', 'ispcrSequenceQualities', isPhred = True),

    SamTag('H0', 'HexamerSequenceR1'),  # anna
    SamTag('H1', 'HexamerPhredScoreR1', isPhred = True),
    SamTag('H2', 'HexamerSequenceR2'),
    SamTag('H3', 'HexamerPhredScoreR2', isPhred = True),

    SamTag('XM', 'perBaseMethylationStatus'),
    SamTag('Cm', 'modified3BPContexts'),
    SamTag('Qm', 'modified4BPContexts'),
    SamTag('Cu', 'unmodified3BPContexts'),
    SamTag('Qu', 'unmodified4BPContexts'),
    SamTag('MC', 'methylatedBaseCountForFragment'),

    SamTag('rS', 'randomPrimerStart'),
    SamTag('rP', 'randomPrimerSequence')

]


TagDefinitions =  { tag.tag: tag for tag in tags }

## Obtain metadata from read

def metaFromRead(read, tag):

    if tag=='chrom':
        return read.reference_name

    if read.has_tag(tag):
        return read.get_tag(tag)

    try:
        return getattr(read, tag)
    except Exception as e:
        pass
        #print(e)

    return None


# Clean a string to be able to be used in a fastq file header
fastqCleanerRegex = re.compile('[^a-zA-Z0-9-_]', re.UNICODE)
def fqSafe(string):
    global fastqCleanerRegex
    return( fastqCleanerRegex.sub('', string) )


illuminaHeaderSplitRegex = re.compile(':| ', re.UNICODE)

class TaggedRecord():
    def __init__(self, tagDefinitions, rawRecord=False, library=None, **kwargs):
        self.tags = {} # 2 character Key -> value
        self.tagDefinitions = tagDefinitions
        if rawRecord is not False:
            try:
                self.fromRawFastq(rawRecord, **kwargs)
            except NonMultiplexable:
                raise
        if library is not None:
            self.addTagByTag( 'LY',library)

        if type(rawRecord) is fastqIterator.FastqRecord:
            self.sequence = rawRecord.sequence
            self.qualities = rawRecord.qual
            self.plus = rawRecord.plus

    def addTagByTag( self,tagName,  value, isPhred=None, decodePhred=False, cast_type=str):

        if isPhred is None:
            isPhred = self.tagDefinitions[tagName].isPhred
        if type(value)!=cast_type:
            value = cast_type(value)
        if isPhred:

            if decodePhred:
                # Convert encoded phred scores back to original ascii
                self.tags[tagName] = fastqHeaderSafeQualitiesToPhred( value, method=3 )
            else:
                self.tags[tagName] = phredToFastqHeaderSafeQualities( value, method=3 )
        else:
            self.tags[tagName] = fqSafe(value)

    def __repr__(self):
        return self.asFastq()

    def asFastq(self, sequence=None, dirAtt=None, baseQualities=None, format='illumina'):
        if sequence is None:
            if self.sequence is None:
                raise ValueError()
            sequence = self.sequence
        if dirAtt is None:
            if self.plus is None:
                raise ValueError()
            dirAtt = self.plus
        if baseQualities is None:
            if self.qualities is None:
                raise ValueError()
            baseQualities = self.qualities

        header = ";".join([ f"{attribute}:{value}" for attribute,value in self.tags.items() if not self.tagDefinitions[attribute].doNotWrite ])
        if len(header)>255: # the header length is stored as uint_8 and includes a null character. The maximum length is thus 255
            raise ValueError("The length of the demultiplexed header is longer than 255 characters. Reduce the length of the header. For example by using -merge _ which will not put the flow cell in the sample name")
        return "@%s\n%s\n%s\n%s\n" % (
            header,
            sequence,
            dirAtt,
            baseQualities
        )

    def has_tag(self, tag):
        return tag in self.tags

    def asIlluminaHeader(self):
        return '{Is}:{RN}:{Fc}:{La}:{Ti}:{CX}:{CY}'.format( **self.tags )


    def fromRawFastq(self, fastqRecord, indexFileParser=None, indexFileAlias=None):
        global illuminaHeaderSplitRegex
        try:
            instrument, runNumber, flowCellId, lane, tile, clusterXpos, clusterYpos, readPairNumber, isFiltered, controlNumber, indexSequence = illuminaHeaderSplitRegex.split(fastqRecord.header.strip())
        except:
            instrument, runNumber, flowCellId, lane, tile, clusterXpos, clusterYpos, readPairNumber, isFiltered, controlNumber = illuminaHeaderSplitRegex.split(fastqRecord.header.strip().replace('::',''))
            indexSequence="N"
            #NS500413:32:H14TKBGXX:2:11101:16448:1664 1:N:0::
        self.addTagByTag( 'Is',instrument)
        self.addTagByTag('RN',runNumber)
        self.addTagByTag('Fc',flowCellId)
        self.addTagByTag('La',lane)
        self.addTagByTag('Ti',tile)
        self.addTagByTag('CX',clusterXpos)
        self.addTagByTag('CY',clusterYpos)
        self.addTagByTag('RP',readPairNumber)
        self.addTagByTag('Fi',isFiltered)
        self.addTagByTag('CN',controlNumber)


        if indexFileParser is not None and indexFileAlias is not None:
            indexIdentifier, correctedIndex, hammingDistance = indexFileParser.getIndexCorrectedBarcodeAndHammingDistance(alias=indexFileAlias, barcode=indexSequence)


            self.addTagByTag('aa',indexSequence)
            if correctedIndex is not None:
                #
                self.addTagByTag('aA',correctedIndex)
                self.addTagByTag('aI',indexIdentifier)
            else:
                raise NonMultiplexable('Could not obtain index for %s  %s %s' % ( indexSequence, correctedIndex, indexIdentifier))
                #self.addTagByTag('aA',"None")
                #self.addTagByTag('aI',-1)
                #self.addTagByTag('ah',hammingDistance)

        else:
            self.addTagByTag('aA',indexSequence)

    def tagPysamRead(self, read):

        moleculeIdentifier = ""
        moleculeQuality = ""
        moleculeIdentifiyingTags = [
            ('BC','QT', False),
            ('RX','RQ',False),
            ('aA',None, True)
        ]
        try:
            for tag,qualityTag, required in moleculeIdentifiyingTags:
                if self.has_tag(tag) and not self.tags.get(tag) is None:
                    moleculeIdentifier += self.tags[tag]
                    if qualityTag is None: # Padding:
                        moleculeQuality+= ( 'o'* len(self.tags[tag]))
                    else:
                        moleculeQuality+=self.tags[qualityTag]

                if required and not self.has_tag(tag):
                    raise NonMultiplexable('Tag was defined to be required')

            if len(moleculeQuality)!=len(moleculeIdentifier):
                raise ValueError('Could not reconstruct molecule identifier')

            self.addTagByTag('MI',moleculeIdentifier, isPhred=False)
            self.addTagByTag('QM',moleculeQuality, isPhred=False)

        except NonMultiplexable:

            # Its bulk
            self.tags['BK'] = True

        # Add sample tag: (If possible)
        # If the BI tag is present it means we know the index of the cell
        # if no BI tag is present, assume the sample is bulk
        if 'BI' in self.tags:
            self.addTagByTag('SM', f'{self.tags["LY"]}_{self.tags["BI"]}')
        else:
            self.addTagByTag('SM', f'{self.tags["LY"]}_BULK')

        # Now we defined the desired values of the tags. Write them to the record:
        for tag,value in self.tags.items():
            #print(tag,value)
            if self.tagDefinitions[tag].isPhred:
                value = fastqHeaderSafeQualitiesToPhred(value, method=3)
            read.set_tag(tag, value)
        if read.has_tag('QM') and len(read.get_tag('QM'))!=len(read.get_tag('MI')):
            raise ValueError('QM and MI tag length not matching')

    def fromTaggedFastq(self, fastqRecord ):
        for keyValue in fastqRecord.header.replace('@','').strip().split(';'):
            key, value = keyValue.split(':')
            self.addTagByTag(key, value, decodePhred=True)

    def fromTaggedBamRecord(self, pysamRecord):
        for keyValue in pysamRecord.query_name.strip().split(';'):
            key, value = keyValue.split(':')
            self.addTagByTag(key, value, isPhred=False)


def reverseComplement(seq):
    global complement
    return( "".join(complement.get(base, base) for base in reversed(seq)) )

def phredToFastqHeaderSafeQualities( asciiEncodedPhredScores, method=3 ):
    """ Convert ASCII encoded pred string to fastq safe string.
    numeric encoded string (method 0),
    or 65 shifted (method 1) which is safe to use in the fastq header"""
    if method==0:
        return( ",".join([str(ord(phred)-33) for phred in asciiEncodedPhredScores]) )
    elif method==1:
        return( "".join([chr(ord(phred)+32) for phred in asciiEncodedPhredScores]) )
    else:
        return( "".join([string.ascii_letters[min(max(0,ord(phred)-33), len(string.ascii_letters))] for phred in asciiEncodedPhredScores]) )

def fastqHeaderSafeQualitiesToPhred(phred, method=3):
    #print(phred)
    return "".join(( chr(string.ascii_letters.index(v)+33) for v in phred ))


class NonMultiplexable(Exception):
    pass

# The demultplexing strategy converts a tuple if fastq record(s) into a demultiplexed record

# Every demultiplexing strategy has :
# a full name: the name shown in the tool
# shortName : a short name, will be put into EVERY fastq record, so keep it really short
# autoDetectable: a flag indicating if this method should be auto detected;
# for example for whole genome sequencing reads, we cannot tell from the reads that it is this data, and the flag should be False

# Method demultiplex( *records (R1, R2, R3 ... )
# Raises NonMultiplexable Exception if the records do not yield a valid result (which is used to determine if the demultplexing method is valid to use)

# Upon initialisation the strategies recieve a dictionary containing the barcodes loaded by the barcode file parser
# (barcodeMapping)

class DemultiplexingStrategy(object):

    def __init__(self):
        self.shortName = 'place holder demultiplexing method'
        self.longName = 'placeHolder'
        self.autoDetectable = False
        self.description = 'inherit this class to build your own demultipexing strategy'
        
        self.indexSummary = ''
        self.barcodeSummary = ''

    def demultiplex(self, records, **kwargs):
        raise NotImplementedError()

    def __repr__(self):
        return f'{self.longName} {self.shortName} {self.description} DemultiplexingStrategy'

    def getParserSummary(self):
        return(' '+self.indexSummary + '\n barcodes:' +self.barcodeSummary )



class IlluminaBaseDemultiplexer(DemultiplexingStrategy):

    def __init__(self, indexFileParser,  indexFileAlias='illumina_merged_ThruPlex48S_RP', **kwargs):
        DemultiplexingStrategy.__init__(self)
        self.indexFileParser = indexFileParser
        self.illuminaIndicesAlias = indexFileAlias
        self.shortName = 'ILLU'
        self.longName = 'IlluminaDemux'
        self.description = 'Demultiplex as a bulk sample'
        self.indexSummary = f'sequencing indices: {indexFileAlias}'

    def demultiplex(self, records, inherited=False, library=None, **kwargs):
        global TagDefinitions

        try:
            if inherited:
                return [ TaggedRecord(rawRecord=record,tagDefinitions=TagDefinitions, indexFileParser=self.indexFileParser, indexFileAlias=self.illuminaIndicesAlias, library=library) for record in records ]
            else:
                return [TaggedRecord(rawRecord=record,tagDefinitions=TagDefinitions,indexFileParser=self.indexFileParser, indexFileAlias=self.illuminaIndicesAlias, library=library).asFastq(record.sequence, record.plus, record.qual) for record in records]
        except NonMultiplexable:
            raise

# Base strategy for read pairs which have both an umi and sample barcode
class UmiBarcodeDemuxMethod(IlluminaBaseDemultiplexer):

    def __init__(self,
        umiRead=0, umiStart = 0, umiLength=6,
        barcodeRead=0, barcodeStart = 6, barcodeLength=8,
         barcodeFileParser=None, barcodeFileAlias=None, indexFileParser=None,
         indexFileAlias = 'illumina_merged_ThruPlex48S_RP',
         **kwargs ):
        self.description=''
        self.barcodeFileAlias = barcodeFileAlias
        self.barcodeFileParser = barcodeFileParser
        IlluminaBaseDemultiplexer.__init__(self, indexFileParser=indexFileParser,indexFileAlias=indexFileAlias)
        self.barcodeSummary = self.barcodeFileAlias
        self.umiRead = umiRead # 0:Read 1, 1: Read 2 etc
        self.umiStart = umiStart # First base
        self.umiLength = umiLength

        self.barcodeRead = barcodeRead
        self.barcodeStart = barcodeStart
        self.barcodeLength = barcodeLength
        self.autoDetectable = False

        self.sequenceCapture = [slice(None) , slice(None) ] # ranges
        if umiLength==0:
            # Barcode only
            if barcodeStart!=0:
                raise NotImplementedError('Complicated slice where we need to capture around a region')
            self.sequenceCapture[barcodeRead] =     slice( barcodeLength, None)
        else:
            if umiRead!=barcodeRead:
                raise NotImplementedError()
            if not( umiStart==0 or barcodeStart==0 ):
                raise NotImplementedError('Complicated slice where we need to capture around a region')
            self.sequenceCapture[barcodeRead] =     slice( barcodeLength+umiLength, None)

    def __repr__(self):
        return f'{self.longName} bc: {self.barcodeStart}:{self.barcodeLength}, umi: {self.umiStart}:{self.umiLength} {self.description}'

    def demultiplex(self, records, **kwargs):

        # Check if the supplied reads are mate-pair:
        if len(records)!=2:
            raise NonMultiplexable('Not mate pair')

        # Perform first pass demultiplexing of the illumina fragments:
        try:
            taggedRecords = IlluminaBaseDemultiplexer.demultiplex(self, records, inherited=True,  **kwargs)
        except NonMultiplexable:
            raise

        rawBarcode = records[self.barcodeRead].sequence[self.barcodeStart:self.barcodeStart+self.barcodeLength]
        barcodeQual =  records[self.barcodeRead].qual[self.barcodeStart:self.barcodeStart+self.barcodeLength]

        barcodeIdentifier, barcode, hammingDistance = self.barcodeFileParser.getIndexCorrectedBarcodeAndHammingDistance(alias=self.barcodeFileAlias, barcode=rawBarcode)
        #print(self.barcodeFileParser,self.barcodeFileAlias,rawBarcode,barcodeIdentifier, barcode, hammingDistance)
        if barcodeIdentifier is None:
            raise NonMultiplexable('barcode not set')

        if self.umiLength!=0:
            umi = records[self.umiRead].sequence[self.umiStart:self.umiStart+self.umiLength]
            umiQual = records[self.umiRead].qual[self.umiStart:self.umiStart+self.umiLength]


        for tr in  taggedRecords:
            #tr.addTagByTag('uL', self.umiLength, isPhred=False)
            if self.umiLength==0:
                #tr.addTagByTag('MI', barcode, isPhred=False)
                #tr.addTagByTag('QM', barcodeQual, isPhred=True)
                pass
            else:
                tr.addTagByTag('RX', umi, isPhred=False)
                tr.addTagByTag('RQ', umiQual, isPhred=True)
                #tr.addTagByTag('MI', barcode+umi, isPhred=False)
                #tr.addTagByTag('QM', barcodeQual+umiQual, isPhred=True)

            tr.addTagByTag('BI', barcodeIdentifier, isPhred=False)
            tr.addTagByTag('bc', rawBarcode, isPhred=False)
            #tr.addTagByTag('hd', hammingDistance, isPhred=False)

            tr.addTagByTag('BC', barcode, isPhred=False )
            tr.addTagByTag('QT', barcodeQual, isPhred=True)

            if len(barcode)!=len(barcodeQual):
                raise ValueError()


            tr.addTagByTag('MX', self.shortName)

        for rid,(record, taggedRecord) in enumerate( zip(records, taggedRecords)):
            taggedRecord.sequence = record.sequence[self.sequenceCapture[rid]]
            taggedRecord.qualities =  record.qual[self.sequenceCapture[rid]]
            taggedRecord.plus = record.plus

        return taggedRecords
        #return [ tr.asFastq(record.sequence[self.sequenceCapture[rid]], record.plus, record.qual[self.sequenceCapture[rid]]) for rid,(tr,record) in enumerate(zip(taggedRecords, records))]
        # Add information and rebuild header
        #header = f'@UMI:{umi};UMIQ:{umiQual};CBI:{barcodeIdentifier};CB:{barcode};CBQ:{barcodeQual};'


        #return fastqIterator.FastqRecord(header, records[1].sequence, records[1].plus,  records[1].qual )
