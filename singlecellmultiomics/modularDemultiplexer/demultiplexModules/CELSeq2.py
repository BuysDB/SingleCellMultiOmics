from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

# Cell seq 2 with 6bp UMI


class CELSeq2_c8_u6(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=6,
            barcodeRead=0,
            barcodeStart=6,
            barcodeLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U6'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 6bp'
        self.autoDetectable = True
        self.description = 'R1 starts with a 6bp UMI  followed by a 8bp cell barcode. R2 ends with a 6bp random primer'


class CELSeq2_c8_u6_NH(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=6,
            barcodeRead=0,
            barcodeStart=6,
            barcodeLength=8,
            random_primer_read=None,
            random_primer_length=None,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U6NH'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 6bp, NO random primer'
        self.autoDetectable = False
        self.description = 'R1 starts with a 6bp UMI  followed by a 8bp cell barcode. R2 has no random primer. Use this demultiplexing method for VASA'

# Reversed case:
class CELSeq2_c8_u6_swapped_reads(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=1,
            umiStart=0,
            umiLength=6,
            barcodeRead=1,
            barcodeStart=6,
            barcodeLength=8,
            random_primer_read=0,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U6S'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 6bp, RP: 6bp'
        self.autoDetectable = True
        self.description = 'R2 starts with a 6bp UMI  followed by a 8bp cell barcode. R1 ends with a 6bp random primer'


# Cell seq 2 with 8bp UMI
class CELSeq2_c8_u8(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeRead=0,
            barcodeStart=8,
            barcodeLength=8,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U8'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 8bp'
        self.autoDetectable = True
        self.description = 'R1 starts with a longer 8bp UMI  followed by a 8bp cell barcode. R2 ends with a 6bp random primer'

# Reversed case:


class CELSeq2_c8_u8(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=1,
            umiStart=0,
            umiLength=8,
            barcodeRead=1,
            barcodeStart=8,
            barcodeLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U8S'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 8bp'
        self.autoDetectable = True
        self.description = 'R2 starts with a longer 8bp UMI  followed by a 8bp cell barcode. R1 ends with a 6bp primer'


class CELSeq2_c8_u8_NNLAIII(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq2_noNla'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=8,
            barcodeRead=0,
            barcodeStart=8,
            barcodeLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS2C8U8NNLA'
        self.longName = 'CELSeq 2, CB: 8bp, UMI: 8bp, NLAIII free'
        self.autoDetectable = False
        self.description = 'CEL-Seq2 without NLAIII digestable barcodes '
