#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import singlecellmultiomics.fastqProcessing.fastqIterator as fastqIterator
import string
from singlecellmultiomics.utils.sequtils import hamming_distance
from singlecellmultiomics.tags import *
from singlecellmultiomics.modularDemultiplexer.headerCodecs import HeaderCodecRegistry
from singlecellmultiomics.tags import TagDefinitions, NonFastqWritableTags

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


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
fastqCleanerRegex = re.compile('[^a-zA-Z0-9-_:.]', re.UNICODE)

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
            tagDefinitions=None,
            rawRecord=False,
            library=None,
            reason=None,
            parse_func=None,
            **kwargs):
        self.tags = {}  # 2 character key -> value
        self.tagDefinitions = tagDefinitions if tagDefinitions is not None else TagDefinitions

        if rawRecord is not False:
            try:
                self.fromRawFastq(rawRecord, parse_func=parse_func, **kwargs)
            except NonMultiplexable:
                raise
        if library is not None and not 'LY' in self.tags:
            self.tags['LY'] = library
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
            make_safe=True): # Make safe: make the tag value safe for fastq headers, only applicable for str type
        
        if isPhred is None:
            isPhred = self.tagDefinitions[tagName].isPhred
        if cast_type and not isinstance(value, cast_type):
            value = cast_type(value)
        if isPhred:
            if decodePhred:
                # Convert encoded phred scores back to original ascii
                self.tags[tagName] = fastqHeaderSafeQualitiesToPhred(value)
            else:
                self.tags[tagName] = phredToFastqHeaderSafeQualities(value)
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

        header = ";".join( [f"{k}:{self.tags[k]}" for k in self.tags if k not in NonFastqWritableTags])
        
        if len(header) > 255:  # the header length is stored as uint_8 and includes a null character. The maximum length is thus 255
            raise ValueError(
                f"The length of the demultiplexed header is longer than 255 characters. Try to keep your library name below 60 characters. Reduce the length of the header. For example by using -merge _ which will not put the flow cell in the sample name. The header looks like this: {header}")
        return f'@{header}\n{sequence}\n{dirAtt}\n{baseQualities}\n'

    def has_tag(self, tag):
        return tag in self.tags

    def asIlluminaHeader(self):
        try:
            return '{Is}:{RN}:{Fc}:{La}:{Ti}:{CX}:{CY}'.format(**self.tags)
        except KeyError as e:
            if 'oh' in self.tags:
                return self.tags['oh']
            raise e

    def fromRawFastq(
            self,
            fastqRecord,
            indexFileParser=None,
            indexFileAlias=None,
            parse_func=None):
        """Parse FASTQ record using header codec.
        
        Args:
            fastqRecord: FastqRecord object to parse
            indexFileParser: Parser for barcode index detection (BarcodeFileParser)
            indexFileAlias: Alias for index file, see modularDemultiplexer/indices (exclude .bc)
            parse_func: Header parsing function (Auto-detects the header format and function if None)
        """
        # Autodetect a header parse function when this is not set
        if parse_func is None:
            registry = HeaderCodecRegistry()
            registry.detect(fastqRecord.header)
            parse_func = registry._parse_func
        
        self.tags.update(parse_func(
            fastqRecord.header,
            indexFileParser,
            indexFileAlias
        ))


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

            self.addTagByTag('MI', moleculeIdentifier, isPhred=False, make_safe=False)
            if not QT_missing:
                self.addTagByTag('QM', moleculeQuality, isPhred=False, make_safe=False)

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
                isPhred=False, make_safe=False)
        elif 'BI' in self.tags:
            self.addTagByTag(
                'SM',
                f'{self.tags["LY"]}_{self.tags["BI"]}',
                isPhred=False, make_safe=False)

            # Remove BI tag
            self.tags['bi'] = self.tags["BI"]
            del self.tags['BI']
        elif "LY" in self.tags:
            self.addTagByTag('SM', f'{self.tags["LY"]}_BULK', isPhred=False, make_safe=False)




        # Now we defined the desired values of the tags. Write them to the
        # record:
        for tag, value in self.tags.items():
            # print(tag,value)
            if tag in self.tagDefinitions and self.tagDefinitions[tag].isPhred:
                value = fastqHeaderSafeQualitiesToPhred(value)
            read.set_tag(tag, value)

        if not QT_missing and read.has_tag('QM') and len(
                read.get_tag('QM')) != len(
                read.get_tag('MI')):
            raise ValueError('QM and MI tag length not matching')

    def fromTaggedFastq(self, fastqRecord):
        for keyValue in fastqRecord.header.replace('@', '').strip().split(';'):
            if ':' not in keyValue:
                continue # This is not a key value apparently...
            # maxsplit=1: tag values may contain colons, e.g. the 'oh' fallback tag
            # (original read header) or rejection reasons written to reject files
            key, value = keyValue.split(':', 1)
            if key in self.tagDefinitions:
                isPhred = self.tagDefinitions[key].isPhred
            else:
                # Unregistered keys are not phred encoded
                isPhred = False
            self.addTagByTag(key, value, isPhred=isPhred, decodePhred=True)

    def fromTaggedBamRecord(self, pysamRecord):
        """Parse header from tagged BAM record (already processed, usually in SCMO key-value format)."""
        try:
            for keyValue in pysamRecord.query_name.strip().split(';'):
                if ':' in keyValue:
                    key, value = keyValue.split(':', 1)
                    value = value.strip()
                    if value.startswith('@'): # Ugly slow backwards compatibility to malformed headers @TODO: remove this at some point
                        value = value[1:]
                    if key.startswith('@'): # Ugly slow backwards compatibility to malformed headers @TODO: remove this at some point
                        key = value[1:]
                    self.addTagByTag(key, value, isPhred=False, make_safe=False) # Slow @TODO
        except ValueError as e:
            # Try to parse "Single Cell Discoveries" header
            # Format: illumina_header;key:value;key:value...
            query_name = pysamRecord.query_name.strip()
            if ';' in query_name:
                illumina_header, attributes = query_name.split(';', 1)
                
                # Parse illumina part using codec @TODO: this is obviously too slow. Need to pass a appropriate codec
                registry = HeaderCodecRegistry()
                if registry.detect('@' + illumina_header): 
                    tags = registry.parse('@' + illumina_header)
                    self.tags.update(tags)
                
                # Parse attributes
                for keyValue in attributes.split(';'):
                    if ':' in keyValue:
                        key, value = keyValue.split(':', 1)
                        value = value.strip()
                        if value.startswith('@'): # Backwards compatibility with incorrectly encoded fastq headers (with @Is)
                            value = value[1:]
                        self.addTagByTag(key, value, isPhred=False, make_safe=False)


def reverseComplement(seq):
    global complement
    return("".join(complement.get(base, base) for base in reversed(seq)))


def phredToFastqHeaderSafeQualities(asciiEncodedPhredScores):
    """ Convert ASCII encoded phred string to a fastq header safe string.

    Quality values are mapped onto string.ascii_letters (Q0 -> 'a' ... Q51 -> 'Z'),
    higher qualities are clamped to the last letter."""
    # Clamp to the last valid index: quality chars above 'T' would otherwise index out of range:
    return("".join([string.ascii_letters[min(max(0, ord(phred) - 33), len(string.ascii_letters) - 1)] for phred in asciiEncodedPhredScores]))


def fastqHeaderSafeQualitiesToPhred(phred):
    """ Decode a header-safe quality string (see phredToFastqHeaderSafeQualities)
    back to ASCII encoded phred scores."""
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
            parse_func=None,
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
                        reason=reason,
                        parse_func=parse_func) for record in records]
            else:
                return [
                    TaggedRecord(
                        rawRecord=record,
                        tagDefinitions=TagDefinitions,
                        indexFileParser=self.indexFileParser,
                        indexFileAlias=self.illuminaIndicesAlias,
                        library=library,
                        reason=reason,
                        parse_func=parse_func).asFastq(
                        record.sequence,
                        record.plus,
                        record.qual) for record in records]
        except NonMultiplexable:
            raise

# Base strategy for read pairs which have both an umi and sample barcode, the barcode and umi are next to another and continuous 
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
        try:
            return f'{self.longName} bc: {self.barcodeStart}:{self.barcodeLength}, umi: {self.umiStart}:{self.umiLength} {self.description}'
        except AttributeError: # Happens when this class is inherited and some of the attributes might have not been set
            return f'{self.longName}, {self.description}'
    

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
                tr.tags['RX'] =  umi
                tr.addTagByTag('RQ', umiQual, isPhred=True,cast_type=None, make_safe=False)
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


# Base strategy for read pairs which have both an umi and sample barcode, but are not contiguous, 
# for example UMI-BCA-UMI_BCB

def apply_slices_seq(record, slices):
    if len(slices)==0 or record is None:
        return ''
    return ''.join( (record.sequence[sc] for sc in slices))

def apply_slices_qual(record, slices):
    if len(slices)==0 or record is None:
        return ''
    return ''.join( (record.qual[sc] for sc in slices))


class ScatteredUmiBarcodeDemuxMethod(IlluminaBaseDemultiplexer):

    def __init__(
            self,
            barcode_slices = None, # 2-len Tuple of List of slices where the cell barcode is present in the read (slice(), slice(), ...], [ slice())
            umi_slices = None, # 2-len List of slices where the UMI is present in the read [slice(), slice(), ...], [slice()] # 
            capture_slices = None, # Slices which indicate what region of read 1 and read2 to store in the final read (required) (slice, slice)
            # Note that there is only one slice per read for the capture
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
        self.barcode_slices = barcode_slices
        self.umi_slices = umi_slices
        self.capture_slices = capture_slices
        self.random_primer_read = random_primer_read
        self.random_primer_length = random_primer_length
        self.random_primer_end = random_primer_end

        

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
        

        # Extract the UMI barcode_slices
        umi, umi_qual = ''.join( (apply_slices_seq(record, slicer) for record, slicer in zip(records, self.umi_slices)) ), ''.join( (apply_slices_qual(record, slicer) for record, slicer in zip(records, self.umi_slices)) )
        raw_barcode, raw_barcode_qual = ''.join( (apply_slices_seq(record, slicer) for record, slicer in zip(records, self.barcode_slices)) ), ''.join( (apply_slices_qual(record, slicer) for record, slicer in zip(records, self.barcode_slices)) )

        barcodeIdentifier, barcode, hammingDistance = self.barcodeFileParser.getIndexCorrectedBarcodeAndHammingDistance(
            alias=self.barcodeFileAlias, barcode=raw_barcode)
        
        if barcodeIdentifier is None:
            raise NonMultiplexable(
                f'bc:{raw_barcode}_not_matching_{self.barcodeFileAlias}')

        random_primer = None
        if self.random_primer_read is not None:
            random_primer = records[self.random_primer_read].sequence[self.random_primer_slice]
        

        for tr in taggedRecords:
            #tr.addTagByTag('uL', self.umiLength, isPhred=False)
            if len(umi)>0:
                tr.tags['RX'] =  umi
                tr.addTagByTag('RQ', umi_qual, isPhred=True,cast_type=None)
               
            tr.tags.update({
                'bi': barcodeIdentifier,
                'bc': raw_barcode,
                'MX': self.shortName,
                'BC': barcode
            })
            if random_primer is not None:
                tr.addTagByTag('rS',
                               random_primer,
                               isPhred=False,
                               make_safe=False)

        for rid, (record, taggedRecord) in enumerate(
                zip(records, taggedRecords)):
            taggedRecord.sequence = record.sequence[self.capture_slices[rid]]
            taggedRecord.qualities = record.qual[self.capture_slices[rid]]
            taggedRecord.plus = record.plus

        return taggedRecords