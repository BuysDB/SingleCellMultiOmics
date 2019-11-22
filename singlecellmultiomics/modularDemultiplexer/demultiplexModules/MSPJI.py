from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

# MSPJI with 3bp UMI


class MSPJI_c8_u3(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'maya_mspj1'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'MSPJIC8U3'
        self.longName = 'MSPJI, CB: 8bp UMI: 3bp'
        self.autoDetectable = True
        self.description = 'MSPJI barcoded fragments. 3bp umi followed by 8bp cell barcode.'
