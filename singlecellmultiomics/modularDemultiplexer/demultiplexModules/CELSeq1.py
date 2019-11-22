from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod


# Cell seq 1 with 6bp UMI
class CELSeq1_c8_u4(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'celseq1'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=8,
            umiLength=4,
            barcodeRead=0,
            barcodeStart=0,
            barcodeLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CS1C8U4'
        self.longName = 'CELSeq 1, CB: 8bp, UMI: 4bp'
        self.autoDetectable = True
        self.description = 'R1 starts with a 8bp cell barcode followed by a 4bp UMI. R2 ends with a 6bp random primer'
