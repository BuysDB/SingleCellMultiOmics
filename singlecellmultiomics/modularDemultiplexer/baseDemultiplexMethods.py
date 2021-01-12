#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import singlecellmultiomics.fastqProcessing.fastqIterator as fastqIterator
import string
from singlecellmultiomics.utils.sequtils import hamming_distance
from singlecellmultiomics.tags import *

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


TagDefinitions = {tag.tag: tag for tag in tags}

# Obtain metadata from read


def metaFromRead(read, tag):

    if tag == 'chrom':
        return read.reference_name

    if read.has_tag(tag):
        return read.get_tag(tag)

    # Backwards compatibility with BI > bi tag:
    if tag=='BI' and read.has_tag('bi'):
        return read.get_tag('bi')

    # Forwards compatibility:
    if tag=='bi' and read.has_tag('BI'):
        return read.get_tag('BI')

    try:
        return getattr(read, tag)
    except Exception as e:
        pass
        # print(e)

    return None


# Clean a string to be able to be used in a fastq file header
fastqCleanerRegex = re.compile('[^a-zA-Z0-9-_]', re.UNICODE)

def fqSafe(string) -> str:
    """
    Convert input string into a representation which can be stored in a fastq header

    Input:
        string(str) : string to clean

    Returns:
        cleaned(str)
    """
    global fastqCleanerRegex
    return(fastqCleanerRegex.sub('', string))


illuminaHeaderSplitRegex = re.compile(':| ', re.UNICODE)


class TaggedRecord():
    def __init__(
            self,
            tagDefinitions,
            rawRecord=False,
            library=None,
            reason=None,
            **kwargs):
        self.tags = {}  # 2 character Key -> value
        self.tagDefinitions = tagDefinitions
        if rawRecord is not False:
            try:
                self.fromRawFastq(rawRecord, **kwargs)
            except NonMultiplexable:
                raise
        if library is not None and not 'LY' in self.tags:
            self.addTagByTag('LY', library, isPhred=False)
        if reason is not None:
            self.tags['RR'] = reason

        if isinstance(rawRecord, fastqIterator.FastqRecord):
            self.sequence = rawRecord.sequence
            self.qualities = rawRecord.qual
            self.plus = rawRecord.plus

    def addTagByTag(
            self,
            tagName,
            value,
            isPhred=None,
            decodePhred=False,
            cast_type=str,
            make_safe=True):

        if isPhred is None:
            isPhred = self.tagDefinitions[tagName].isPhred
        if not isinstance(value, cast_type):
            value = cast_type(value)
        if isPhred:
            if decodePhred:
                # Convert encoded phred scores back to original ascii
                self.tags[tagName] = fastqHeaderSafeQualitiesToPhred(
                    value, method=3)
            else:
                self.tags[tagName] = phredToFastqHeaderSafeQualities(
                    value, method=3)
        else:
            if cast_type is str:
                if make_safe:
                    self.tags[tagName] = fqSafe(value)
                else:
                    self.tags[tagName] = value
            else:
                self.tags[tagName] = value

    def __repr__(self):
        return self.asFastq()

    def asFastq(
            self,
            sequence=None,
            dirAtt=None,
            baseQualities=None,
            format='illumina'):
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

        header = ";".join([f"{attribute}:{value}" for attribute, value in self.tags.items(
        ) if not self.tagDefinitions[attribute].doNotWrite])
        if len(header) > 255:  # the header length is stored as uint_8 and includes a null character. The maximum length is thus 255
            raise ValueError(
                f"The length of the demultiplexed header is longer than 255 characters. Try to keep your library name below 60 characters. Reduce the length of the header. For example by using -merge _ which will not put the flow cell in the sample name. The header looks like this: {header}")
        return f'@{header}\n{sequence}\n{dirAtt}\n{baseQualities}\n'

    def has_tag(self, tag):
        return tag in self.tags

    def asIlluminaHeader(self):
        return '{Is}:{RN}:{Fc}:{La}:{Ti}:{CX}:{CY}'.format(**self.tags)




    def parse_3dec_header(self,fastqRecord, indexFileParser,  indexFileAlias):

        instrument = 'UNK'
        runNumber = 'UNK'
        flowCellId = 'UNK'
        indexSequence = 'N'
        lane = 'UNK'
        tile = 'UNK'
        clusterXpos = '-1'
        clusterYpos = '-1'
        readPairNumber = '0'
        isFiltered = '0'
        controlNumber = '0'

        # 3-DEC: @Cluster_s_1_1101_2
        if fastqRecord.header.count('_') == 4:
            _cluster_, _s_, lane, tile, readPairNumber = fastqRecord.header.split(
                '_')
            # check  that this s thingy is at the right place
            assert(_s_ == 's')
        else:
            raise

        self.tags.update({
            'Is': instrument,
            'RN': runNumber,
            'Fc': flowCellId,
            'La': lane,
            'Ti': tile,
            'CX': clusterXpos,
            'CY': clusterYpos,
            'RP': readPairNumber,
            'Fi': isFiltered,
            'CN': controlNumber
        })

    def parse_illumina_header(self,fastqRecord, indexFileParser,  indexFileAlias):
        try:
            instrument, runNumber, flowCellId, lane, tile, clusterXpos, clusterYpos, readPairNumber, isFiltered, controlNumber, indexSequence = illuminaHeaderSplitRegex.split(
                fastqRecord.header.strip())
        except BaseException:

            try:
                instrument, runNumber, flowCellId, lane, tile, clusterXpos, clusterYpos, readPairNumber, isFiltered, controlNumber = illuminaHeaderSplitRegex.split(
                    fastqRecord.header.strip().replace('::', ''))
                indexSequence = "N"
            except BaseException:
                raise

        self.tags.update({
            'Is': instrument,
            'RN': runNumber,
            'Fc': flowCellId,
            'La': lane,
            'Ti': tile,
            'CX': clusterXpos,
            'CY': clusterYpos,
            'RP': readPairNumber,
            'Fi': isFiltered,
            'CN': controlNumber
        })

        if indexFileParser is not None and indexFileAlias is not None:
            # Check if the index is an integer:
            try:
                indexInteger = int(indexSequence)
                indexIdentifier, correctedIndex, hammingDistance = indexSequence, indexSequence, 0
            except Exception:
                indexIdentifier, correctedIndex, hammingDistance = indexFileParser.getIndexCorrectedBarcodeAndHammingDistance(
                    alias=indexFileAlias, barcode=indexSequence)

            self.addTagByTag('aa', indexSequence, isPhred=False)
            if correctedIndex is not None:
                #
                #self.addTagByTag('aA',correctedIndex, isPhred=False)
                #self.addTagByTag('aI',indexIdentifier, isPhred=False)
                self.tags.update({'aA': correctedIndex, 'aI': indexIdentifier})
            else:
                raise NonMultiplexable(
                    'Could not obtain index for %s  %s %s' %
                    (indexSequence, correctedIndex, indexIdentifier))
                # self.addTagByTag('aA',"None")
                # self.addTagByTag('aI',-1)
                # self.addTagByTag('ah',hammingDistance)

        else:
            #self.addTagByTag('aA',indexSequence, isPhred=False)
            self.tags['aa'] = indexSequence

    def parse_scmo_header(self, fastqRecord, indexFileParser,  indexFileAlias):
        self.tags.update( dict( kv.split(':') for kv in fastqRecord.header.strip()[1:].split(';') ) )

    def fromRawFastq(
            self,
            fastqRecord,
            indexFileParser=None,
            indexFileAlias=None):

        try:
            self.parse_illumina_header(fastqRecord, indexFileParser,  indexFileAlias)
        except BaseException:
            if fastqRecord.header.startswith('@Is'):
                self.parse_scmo_header(fastqRecord, indexFileParser,  indexFileAlias)
            else:
                self.parse_3dec_header(fastqRecord, indexFileParser,  indexFileAlias)


            # NS500413:32:H14TKBGXX:2:11101:16448:1664 1:N:0::
        """ This is the nice and safe way:
        self.addTagByTag( 'Is',instrument, isPhred=False)
        self.addTagByTag('RN',runNumber, isPhred=False)
        self.addTagByTag('Fc',flowCellId, isPhred=False)
        self.addTagByTag('La',lane, isPhred=False)
        self.addTagByTag('Ti',tile, isPhred=False)
        self.addTagByTag('CX',clusterXpos, isPhred=False)
        self.addTagByTag('CY',clusterYpos, isPhred=False)
        self.addTagByTag('RP',readPairNumber, isPhred=False)
        self.addTagByTag('Fi',isFiltered, isPhred=False)
        self.addTagByTag('CN',controlNumber, isPhred=False)
        """

    def tagPysamRead(self, read):

        moleculeIdentifier = ""
        moleculeQuality = ""
        moleculeIdentifiyingTags = [
            ('BC', 'QT', False),
            ('RX', 'RQ', False),
            ('aA', None, True)
        ]
        try:
            QT_missing=False
            for tag, qualityTag, required in moleculeIdentifiyingTags:
                if self.has_tag(tag) and not self.tags.get(tag) is None:
                    moleculeIdentifier += self.tags[tag]
                    if qualityTag is None:  # Padding:
                        moleculeQuality += ('o' * len(self.tags[tag]))
                    else:
                        if qualityTag in self.tags:
                            moleculeQuality += self.tags[qualityTag]
                        elif qualityTag=='QT':
                            QT_missing = True

                if required and not self.has_tag(tag):
                    raise NonMultiplexable('Tag was defined to be required')

            # if len(moleculeQuality)!=len(moleculeIdentifier):
            # raise ValueError('Could not reconstruct molecule identifier')
            # @todo set this back when we recover QT

            correctedIndex = None if not self.has_tag(
                'aA') else self.tags['aA']
            indexSequence = None if not self.has_tag('aa') else self.tags['aa']

            if correctedIndex is not None and indexSequence is not None:
                hd = hamming_distance(indexSequence, correctedIndex)
                if hd is None:
                    raise ValueError(
                        "Could not resolve hamming distance between {correctedIndex} and {indexSequence}")
                self.addTagByTag('ah', hd, isPhred=False, cast_type=int)

            self.addTagByTag('MI', moleculeIdentifier, isPhred=False)
            if not QT_missing:
                self.addTagByTag('QM', moleculeQuality, isPhred=False)

        except NonMultiplexable:

            # Its bulk
            self.tags['BK'] = True

        # Add sample tag: (If possible)
        # If the bi tag is present it means we know the index of the cell
        # if no bi tag is present, assume the sample is bulk
        if 'bi' in self.tags:
            self.addTagByTag(
                'SM',
                f'{self.tags["LY"]}_{self.tags["bi"]}',
                isPhred=False)
        elif 'BI' in self.tags:
            self.addTagByTag(
                'SM',
                f'{self.tags["LY"]}_{self.tags["BI"]}',
                isPhred=False)

            # Remove BI tag
            self.tags['bi'] = self.tags["BI"]
            del self.tags['BI']
        else:
            self.addTagByTag('SM', f'{self.tags["LY"]}_BULK', isPhred=False)




        # Now we defined the desired values of the tags. Write them to the
        # record:
        for tag, value in self.tags.items():
            # print(tag,value)
            if self.tagDefinitions[tag].isPhred:
                value = fastqHeaderSafeQualitiesToPhred(value, method=3)
            read.set_tag(tag, value)

        if not QT_missing and read.has_tag('QM') and len(
                read.get_tag('QM')) != len(
                read.get_tag('MI')):
            raise ValueError('QM and MI tag length not matching')

    def fromTaggedFastq(self, fastqRecord):
        for keyValue in fastqRecord.header.replace('@', '').strip().split(';'):
            key, value = keyValue.split(':')
            self.addTagByTag(key, value, decodePhred=True)

    def fromTaggedBamRecord(self, pysamRecord):
        for keyValue in pysamRecord.query_name.strip().split(';'):
            key, value = keyValue.split(':')
            self.addTagByTag(key, value, isPhred=False)


def reverseComplement(seq):
    global complement
    return("".join(complement.get(base, base) for base in reversed(seq)))


def phredToFastqHeaderSafeQualities(asciiEncodedPhredScores, method=3):
    """ Convert ASCII encoded pred string to fastq safe string.
    numeric encoded string (method 0),
    or 65 shifted (method 1) which is safe to use in the fastq header"""
    if method == 0:
        return(",".join([str(ord(phred) - 33) for phred in asciiEncodedPhredScores]))
    elif method == 1:
        return("".join([chr(ord(phred) + 32) for phred in asciiEncodedPhredScores]))
    else:
        return("".join([string.ascii_letters[min(max(0, ord(phred) - 33), len(string.ascii_letters))] for phred in asciiEncodedPhredScores]))


def fastqHeaderSafeQualitiesToPhred(phred, method=3):
    # print(phred)
    return "".join((chr(string.ascii_letters.index(v) + 33) for v in phred))


class NonMultiplexable(Exception):
    pass

# The demultplexing strategy converts a tuple if fastq record(s) into a
# demultiplexed record

# Every demultiplexing strategy has :
# a full name: the name shown in the tool
# shortName : a short name, will be put into EVERY fastq record, so keep it really short
# autoDetectable: a flag indicating if this method should be auto detected;
# for example for whole genome sequencing reads, we cannot tell from the
# reads that it is this data, and the flag should be False

# Method demultiplex( *records (R1, R2, R3 ... )
# Raises NonMultiplexable Exception if the records do not yield a valid
# result (which is used to determine if the demultplexing method is valid
# to use)

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
        return(' ' + self.indexSummary + '\n barcodes:' + self.barcodeSummary)


class IlluminaBaseDemultiplexer(DemultiplexingStrategy):

    def __init__(
            self,
            indexFileParser,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',
            **kwargs):
        DemultiplexingStrategy.__init__(self)
        self.indexFileParser = indexFileParser
        self.illuminaIndicesAlias = indexFileAlias
        self.shortName = 'ILLU'
        self.longName = 'IlluminaDemux'
        self.description = 'Demultiplex as a bulk sample'
        self.indexSummary = f'sequencing indices: {indexFileAlias}'

    def demultiplex(
            self,
            records,
            inherited=False,
            library=None,
            reason=None,
            **kwargs):
        global TagDefinitions

        try:
            if inherited:
                return [
                    TaggedRecord(
                        rawRecord=record,
                        tagDefinitions=TagDefinitions,
                        indexFileParser=self.indexFileParser,
                        indexFileAlias=self.illuminaIndicesAlias,
                        library=library,
                        reason=reason) for record in records]
            else:
                return [
                    TaggedRecord(
                        rawRecord=record,
                        tagDefinitions=TagDefinitions,
                        indexFileParser=self.indexFileParser,
                        indexFileAlias=self.illuminaIndicesAlias,
                        library=library,
                        reason=reason).asFastq(
                        record.sequence,
                        record.plus,
                        record.qual) for record in records]
        except NonMultiplexable:
            raise

# Base strategy for read pairs which have both an umi and sample barcode


class UmiBarcodeDemuxMethod(IlluminaBaseDemultiplexer):

    def __init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=6,
            barcodeRead=0,
            barcodeStart=6,
            barcodeLength=8,
            barcodeFileParser=None,
            barcodeFileAlias=None,
            indexFileParser=None,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',
            random_primer_read=None,
            random_primer_length=6,
            random_primer_end=False, # True for end, False for start
            **kwargs):
        self.description = ''
        self.barcodeFileAlias = barcodeFileAlias
        self.barcodeFileParser = barcodeFileParser
        IlluminaBaseDemultiplexer.__init__(
            self,
            indexFileParser=indexFileParser,
            indexFileAlias=indexFileAlias)
        self.barcodeSummary = self.barcodeFileAlias
        self.umiRead = umiRead  # 0:Read 1, 1: Read 2 etc
        self.umiStart = umiStart  # First base
        self.umiLength = umiLength

        self.barcodeRead = barcodeRead
        self.barcodeStart = barcodeStart
        self.barcodeLength = barcodeLength
        self.autoDetectable = False

        self.random_primer_read = random_primer_read
        self.random_primer_length = random_primer_length
        self.random_primer_end = random_primer_end

        # ranges to capture for read 1 and read 2
        self.sequenceCapture = [slice(None), slice(None)]
        if umiLength == 0:
            # Barcode only
            if barcodeStart != 0:
                raise NotImplementedError(
                    'Complicated slice where we need to capture around a region')
            self.sequenceCapture[barcodeRead] = slice(barcodeLength, None)
        else:
            if umiRead != barcodeRead:
                raise NotImplementedError()
            if not(umiStart == 0 or barcodeStart == 0):
                raise NotImplementedError(
                    'Complicated slice where we need to capture around a region')
            self.sequenceCapture[barcodeRead] = slice(
                barcodeLength + umiLength, None)

        if random_primer_read is not None:

            if self.sequenceCapture[random_primer_read].stop is not None:
                raise NotImplementedError()
            if random_primer_end:

                self.sequenceCapture[random_primer_read] = slice(
                    self.sequenceCapture[random_primer_read].start,
                    -random_primer_length,
                    self.sequenceCapture[random_primer_read].step
                )
                self.random_primer_slice = slice(-random_primer_length, None, None)
            else:
                self.sequenceCapture[random_primer_read] = slice(
                    random_primer_length,
                    None,
                    self.sequenceCapture[random_primer_read].step
                )
                self.random_primer_slice = slice(0, random_primer_length, None)

    def __repr__(self):
        return f'{self.longName} bc: {self.barcodeStart}:{self.barcodeLength}, umi: {self.umiStart}:{self.umiLength} {self.description}'

    def demultiplex(self, records, **kwargs):

        # Check if the supplied reads are mate-pair or single end
        if len(records) not in (1, 2):
            raise NonMultiplexable('Not mate pair or single end')


        # Perform first pass demultiplexing of the illumina fragments:
        try:
            taggedRecords = IlluminaBaseDemultiplexer.demultiplex(
                self, records, inherited=True, **kwargs)
        except NonMultiplexable:
            raise

        rawBarcode = records[self.barcodeRead].sequence[self.barcodeStart:
                                                        self.barcodeStart + self.barcodeLength]
        barcodeQual = records[self.barcodeRead].qual[self.barcodeStart:
                                                     self.barcodeStart + self.barcodeLength]

        barcodeIdentifier, barcode, hammingDistance = self.barcodeFileParser.getIndexCorrectedBarcodeAndHammingDistance(
            alias=self.barcodeFileAlias, barcode=rawBarcode)
        #print(self.barcodeFileParser,self.barcodeFileAlias,rawBarcode,barcodeIdentifier, barcode, hammingDistance)
        if barcodeIdentifier is None:
            raise NonMultiplexable(
                f'bc:{rawBarcode}_not_matching_{self.barcodeFileAlias}')

        random_primer = None
        if self.random_primer_read is not None:
            random_primer = records[self.random_primer_read].sequence[self.random_primer_slice]
        if self.umiLength != 0:
            umi = records[self.umiRead].sequence[self.umiStart:self.umiStart + self.umiLength]
            umiQual = records[self.umiRead].qual[self.umiStart:self.umiStart + self.umiLength]

        for tr in taggedRecords:
            #tr.addTagByTag('uL', self.umiLength, isPhred=False)
            if self.umiLength == 0:
                #tr.addTagByTag('MI', barcode, isPhred=False)
                #tr.addTagByTag('QM', barcodeQual, isPhred=True)
                pass
            else:
                tr.addTagByTag('RX', umi, isPhred=False, make_safe=False)
                tr.addTagByTag('RQ', umiQual, isPhred=True)
                #tr.addTagByTag('MI', barcode+umi, isPhred=False)
                #tr.addTagByTag('QM', barcodeQual+umiQual, isPhred=True)

            """ These can be updated at once
            tr.addTagByTag('bi', barcodeIdentifier, isPhred=False)
            tr.addTagByTag('bc', rawBarcode, isPhred=False)
            tr.addTagByTag('MX', self.shortName, isPhred=False)
            tr.addTagByTag('BC', barcode, isPhred=False )
            """
            tr.tags.update({
                'bi': barcodeIdentifier,
                'bc': rawBarcode,
                'MX': self.shortName,
                'BC': barcode
            })
            #tr.addTagByTag('hd', hammingDistance, isPhred=False)
            if random_primer is not None:
                tr.addTagByTag('rS',
                               random_primer,
                               isPhred=False,
                               make_safe=False)

            #tr.addTagByTag('QT', barcodeQual, isPhred=True)

            if len(barcode) != len(barcodeQual):
                raise ValueError()

        for rid, (record, taggedRecord) in enumerate(
                zip(records, taggedRecords)):
            taggedRecord.sequence = record.sequence[self.sequenceCapture[rid]]
            taggedRecord.qualities = record.qual[self.sequenceCapture[rid]]
            taggedRecord.plus = record.plus

        return taggedRecords
        # return [ tr.asFastq(record.sequence[self.sequenceCapture[rid]], record.plus, record.qual[self.sequenceCapture[rid]]) for rid,(tr,record) in enumerate(zip(taggedRecords, records))]
        # Add information and rebuild header
        #header = f'@UMI:{umi};UMIQ:{umiQual};CBI:{barcodeIdentifier};CB:{barcode};CBQ:{barcodeQual};'

        # return fastqIterator.FastqRecord(header, records[1].sequence,
        # records[1].plus,  records[1].qual )
