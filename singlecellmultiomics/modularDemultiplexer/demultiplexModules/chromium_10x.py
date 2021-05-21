from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

class chrom10x_c16_u12(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = '10x_3M-february-2018'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=16,
            umiLength=12,
            barcodeRead=0,
            barcodeStart=0,
            barcodeLength=16,
            random_primer_read=None,
            random_primer_length=None,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'CHROMC16U12'
        self.longName = 'Chromium 10x, CB: 16bp, UMI: 12bp'
        self.autoDetectable = False
        self.description = 'R1 starts with a 16bp cell barcode followed by a 12bp UMI.'
