from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import DemultiplexingStrategy
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TaggedRecord
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TagDefinitions


# ask buys or annaa about it

class HexamerBaseDemultiplexer(DemultiplexingStrategy):

    def __init__(
            self,
            indexFileParser,
            illuminaIndicesAlias='illumina_merged_iPCR_RP',
            **kwargs):
        self.barcodeFileParser = None
        DemultiplexingStrategy.__init__(self)
        self.indexFileParser = indexFileParser
        self.illuminaIndicesAlias = illuminaIndicesAlias
        self.shortName = 'HEX'
        self.longName = 'HexamersDemux'
        self.autoDetectable = False
        self.hexLength = 6
        self.description = 'Demultiplex as a bulk sample and add hexamer H0..H3 tags'
        self.barcodeSummary = 'Bulk, no cell barcodes'
        self.indexSummary = f'sequencing indices: {illuminaIndicesAlias}'

    def demultiplex(self, records, library=None, **kwargs):
        global TagDefinitions

        try:
            h0 = records[0].sequence[:self.hexLength]
            h1 = records[0].qual[:self.hexLength]
            if len(records) > 1:
                h2 = records[1].sequence[:self.hexLength]
                h3 = records[1].qual[:self.hexLength]
            trs = []
            for record in records:
                t = TaggedRecord(
                    rawRecord=record,
                    tagDefinitions=TagDefinitions,
                    indexFileParser=self.indexFileParser,
                    indexFileAlias=self.illuminaIndicesAlias,
                    library=library)
                t.addTagByTag('H0', h0, isPhred=False)
                t.addTagByTag('H1', h1, isPhred=True)
                if len(records) > 1:
                    t.addTagByTag('H2', h2, isPhred=False)
                    t.addTagByTag('H3', h3, isPhred=True)
                trs.append(
                    t.asFastq(
                        record.sequence,
                        record.plus,
                        record.qual))

            return trs

        except NonMultiplexable:
            raise
        except Exception as e:
            print(e)
