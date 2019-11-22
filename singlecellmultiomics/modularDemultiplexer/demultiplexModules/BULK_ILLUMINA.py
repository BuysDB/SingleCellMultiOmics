from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import DemultiplexingStrategy
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TaggedRecord
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TagDefinitions
import singlecellmultiomics.barcodeFileParser.barcodeFileParser as barcodeFileParser


class IlluminaBaseDemultiplexer(DemultiplexingStrategy):

    def __init__(
            self,
            indexFileParser,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',
            **kwargs):
        self.barcodeFileParser = None
        DemultiplexingStrategy.__init__(self)
        self.indexFileParser = indexFileParser

        self.illuminaIndicesAlias = indexFileAlias
        self.shortName = 'ILLU'
        self.longName = 'IlluminaDemux'
        self.autoDetectable = False
        self.description = 'Demultiplex as a bulk sample'
        self.barcodeSummary = 'Bulk, no cell barcodes'
        self.indexSummary = f'sequencing indices: {indexFileAlias}'

    def demultiplex(self, records, library=None, reason=None, **kwargs):
        global TagDefinitions

        try:
            return [
                TaggedRecord(
                    rawRecord=record,
                    tagDefinitions=TagDefinitions,
                    indexFileParser=self.indexFileParser,
                    indexFileAlias=self.illuminaIndicesAlias,
                    library=library,
                    reason=reason) for record in records]
            #[TaggedRecord(rawRecord=record,tagDefinitions=TagDefinitions,indexFileParser=self.indexFileParser, indexFileAlias=self.illuminaIndicesAlias, library=library).asFastq(record.sequence, record.plus, record.qual) for record in records]
        except NonMultiplexable:
            raise
        except Exception as e:
            raise
