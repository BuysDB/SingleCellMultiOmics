from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod, NonMultiplexable

# NLAIII, 96 well format with 3bp UMI


class NLAIII_96w_c8_u3(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'lennart96NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            random_primer_read=1,
            random_primer_length=6,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'NLAIII96C8U3'
        self.longName = 'NLAIII, 96well CB: 8bp UMI: 3bp RP: 6bp'
        self.autoDetectable = True
        self.description = '96 well format. 3bp umi followed by 8bp barcode. R2 ends with a 6bp random primer'


# NLAIII, 384 well format with 3bp UMI
class NLAIII_384w_c8_u3(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=1,
            random_primer_length=6,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'NLAIII384C8U3'
        self.longName = 'NLAIII, 384well CB: 8bp UMI: 3bp RP:6bp'
        self.autoDetectable = True
        self.description = '384 well format. 3bp umi followed by 8bp barcode. R2 ends with a 6bp random primer'

    def demultiplex(self, records, **kwargs):
        if kwargs.get('probe') and records[0].sequence[self.barcodeLength + \
                      self.umiLength: self.barcodeLength + self.umiLength + 4] != 'CATG':
            raise NonMultiplexable

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, **kwargs)
        return taggedRecords
