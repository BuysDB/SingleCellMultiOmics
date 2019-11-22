from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

# ScarTrace


class ScartraceR1(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'scartrace'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=0,
            barcodeRead=0,
            barcodeStart=0,
            barcodeLength=8,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'SCARC8R1'
        self.longName = 'Scartrace, CB: 8bp'
        self.description = '384 well format. Scar amplicon demultiplexing, cell barcode in read 1'
        self.autoDetectable = True


class ScartraceR2(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'scartrace'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=0,
            barcodeRead=1,
            barcodeStart=0,
            barcodeLength=8,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'SCARC8R2'
        self.longName = 'Scartrace, CB: 8bp'
        self.description = '384 well format. Scar amplicon demultiplexing, cell barcode in read 2'
        self.autoDetectable = True
