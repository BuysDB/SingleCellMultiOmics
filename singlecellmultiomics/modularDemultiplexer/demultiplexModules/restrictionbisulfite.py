import sys
# sys.path.append("/hpc/hub_oudenaarden/group_scripts/")
# sys.path.append("/hpc/hub_oudenaarden/group_scripts/modularDemultiplexer/")

from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod, IlluminaBaseDemultiplexer

from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import IlluminaBaseDemultiplexer
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable

# TO do, currently not implemented:
# enzyme ID correction based on Hamming Distance
# same for ISPCR sequence


class Base_RestrictionBisulfiteDemuxMethod(UmiBarcodeDemuxMethod):

    def __init__(self,
                 umiRead=0, umiStart=0, umiLength=8,  # default settings UMI
                 barcodeRead=0, barcodeStart=8, barcodeLength=8,  # default settings Barcode
                 enzymeRead=0, enzymeStart=16, enzymeLength=3,  # default settings Enzyme ID
                 ispcrRead=0, ispcrStart=19, ispcrLength=15,  # default settings ISPCR
                 ispcrSeq="CAGTGGTATCAGAGT",
                 barcodeFileParser=None,  # compatible, no need to change
                 barcodeFileAlias=None,  # passed from lower-level Classes, e.g. "reBS_nla384w"
                 indexFileParser=None,  # compatible, no need to change
                 **kwargs):  # additional arguments
        self.description = 'base class for restriction bisulfite'
        self.barcodeFileAlias = barcodeFileAlias  # description , e.g. "maya_384NLA"
        self.barcodeFileParser = barcodeFileParser  # Namespace for barcode file parse
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=umiRead,
            umiStart=umiStart,
            umiLength=umiLength,
            barcodeRead=barcodeRead,
            barcodeStart=barcodeStart,
            barcodeLength=barcodeLength,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)

        self.barcodeSummary = self.barcodeFileAlias
        self.umiRead = umiRead  # 0:Read 1, 1: Read 2 etc
        self.umiStart = umiStart  # First base
        self.umiLength = umiLength
        self.shortName = 'RB'
        self.longName = 'base class for restriction bisulfite'
        self.illumina_mux = IlluminaBaseDemultiplexer(
            indexFileParser=indexFileParser,
            indexFileAlias='illumina_merged_ThruPlex48S_RP')

        self.barcodeRead = barcodeRead
        self.barcodeStart = barcodeStart
        self.barcodeLength = barcodeLength

        self.enzymeRead = enzymeRead
        self.enzymeStart = enzymeStart
        self.enzymeLength = enzymeLength

        self.ispcrRead = ispcrRead
        self.ispcrStart = ispcrStart
        self.ispcrLength = ispcrLength

        self.autoDetectable = False

        self.sequenceCapture = [slice(None), slice(None)]  # ranges
        # TAKE OUT IF STATEMENT
        if umiLength == 0:
            # if there is a barcode only
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
                barcodeLength + umiLength + enzymeLength + ispcrLength, None)

    def __repr__(self):
        return f'{self.longName} {self.description}'

    def demultiplex(self, records, **kwargs):

        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')

        # Perform first pass demultiplexing of the illumina fragments:
        try:
            taggedRecords = self.illumina_mux.demultiplex(
                records, inherited=True, **kwargs)
        except NonMultiplexable:
            raise
        rawBarcode = records[self.barcodeRead].sequence[self.barcodeStart:
                                                        self.barcodeStart + self.barcodeLength]
        barcodeQual = records[self.barcodeRead].qual[self.barcodeStart:
                                                     self.barcodeStart + self.barcodeLength]

        barcodeIdentifier, barcode, hammingDistance = self.barcodeFileParser.getIndexCorrectedBarcodeAndHammingDistance(
            alias=self.barcodeFileAlias, barcode=rawBarcode)
        #print(barcodeIdentifier, barcode, hammingDistance)
        if barcodeIdentifier is None:
            raise NonMultiplexable('Illumina barcode not set')

        if self.umiLength != 0:
            umi = records[self.umiRead].sequence[self.umiStart:self.umiStart + self.umiLength]
            umiQual = records[self.umiRead].qual[self.umiStart:self.umiStart + self.umiLength]

        if self.enzymeLength != 0:
            enz = records[self.enzymeRead].sequence[self.enzymeStart:
                                                    self.enzymeStart + self.enzymeLength]
            enzQual = records[self.enzymeRead].qual[self.enzymeStart:
                                                    self.enzymeStart + self.enzymeLength]

        if self.ispcrLength != 0:
            ispcr = records[self.ispcrRead].sequence[self.ispcrStart:
                                                     self.ispcrStart + self.ispcrLength]
            ispcrQual = records[self.ispcrRead].qual[self.ispcrStart:
                                                     self.ispcrStart + self.ispcrLength]

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

            tr.addTagByTag(
                'bi',
                barcodeIdentifier,
                isPhred=False,
                make_safe=False)
            tr.addTagByTag('bc', rawBarcode, isPhred=False, make_safe=False)
            #tr.addTagByTag('hd', hammingDistance, isPhred=False)

            tr.addTagByTag('BC', barcode, isPhred=False, make_safe=False)
            tr.addTagByTag('QT', barcodeQual, isPhred=True)

            if len(barcode) != len(barcodeQual):
                raise ValueError()

         # Enzyme sequence
            tr.addTagByTag(
                'ES',
                enz,
                isPhred=False,
                make_safe=False)  # Enzyme sequence
         # Enzyme quality
            tr.addTagByTag('eq', enzQual, isPhred=True)

            # ISPCR sequence
            tr.addTagByTag('IS', ispcr, isPhred=False, make_safe=False)
            # ISPCR quality
            #tr.addTagByTag('is', ispcrQual, isPhred=True)

            tr.addTagByTag(
                'MX',
                self.shortName,
                make_safe=False,
                isPhred=False)

        for rid, (record, taggedRecord) in enumerate(
                zip(records, taggedRecords)):
            taggedRecord.sequence = record.sequence[self.sequenceCapture[rid]]
            taggedRecord.qualities = record.qual[self.sequenceCapture[rid]]
            taggedRecord.plus = record.plus

        return taggedRecords
        # Add information and rebuild header
        #header = f'@UMI:{umi};UMIQ:{umiQual};CBI:{barcodeIdentifier};CB:{barcode};CBQ:{barcodeQual};'

        # return fastqIterator.FastqRecord(header, records[1].sequence,
        # records[1].plus,  records[1].qual )


# NLAIII, 384 well format with 3bp UMI
class Nla_384w_u8_c8_ad3_is15(Base_RestrictionBisulfiteDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'nla_bisulfite'
        Base_RestrictionBisulfiteDemuxMethod.__init__(
            self,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'RBSN' #'ReBsNla384C8U8E3I15'
        self.longName = 'Restriction-BS, NlaIII; 384w; UMI: 8 bp, CB: 8bp, Enz. ID: 3bp, ISPCR: 15 bp'
        self.autoDetectable = True
        self.description = 'Bisulfite-compatible barcoded Nla Adapters (384w). R1 contains UMI, BC, Enzyme ID and ISPCR.'
