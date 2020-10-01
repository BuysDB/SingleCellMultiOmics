from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod, NonMultiplexable

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

    def demultiplex(self, records, **kwargs):

        if kwargs.get(
                'probe') and not records[0].sequence.startswith('CCTTGAACTTCTGGTTGTAG'):
            raise NonMultiplexable

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
                self, records, **kwargs)
        return taggedRecords



class ScartraceR2RP4(UmiBarcodeDemuxMethod):
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

            random_primer_end=False,
            random_primer_read=0,
            random_primer_length=4,

            **kwargs)
        self.shortName = 'SCARC8R2R4'
        self.longName = 'Scartrace, CB: 8bp, with 4bp random sequence in read 1'
        self.description = '384 well format. Scar amplicon demultiplexing, cell barcode in read , 4bp random sequence in R1'
        self.autoDetectable = True

    def demultiplex(self, records, **kwargs):

        if kwargs.get(
                'probe') and not records[0].sequence[4:].startswith('CCTTGAACTTCTGGTTGTAG'):
            raise NonMultiplexable

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
                self, records, **kwargs)
        return taggedRecords
