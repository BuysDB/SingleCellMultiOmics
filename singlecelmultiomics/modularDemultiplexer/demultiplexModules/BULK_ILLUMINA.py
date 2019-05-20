from baseDemultiplexMethods import UmiBarcodeDemuxMethod
from baseDemultiplexMethods import DemultiplexingStrategy
from baseDemultiplexMethods import NonMultiplexable
from baseDemultiplexMethods import TaggedRecord
from baseDemultiplexMethods import TagDefinitions
import barcodeFileParser


class IlluminaBaseDemultiplexer(DemultiplexingStrategy):

	def __init__(self, indexFileParser,  illuminaIndicesAlias='illumina_merged_iPCR_RP', **kwargs):
		self.barcodeFileParser = None
		DemultiplexingStrategy.__init__(self)
		self.indexFileParser = indexFileParser

		self.illuminaIndicesAlias = illuminaIndicesAlias
		self.shortName = 'ILLU'
		self.longName = 'IlluminaDemux'
		self.autoDetectable = False
		self.description = 'Demultiplex as a bulk sample'
		self.barcodeSummary='Bulk, no cell barcodes'
		self.indexSummary = f'sequencing indices: {illuminaIndicesAlias}'

	def demultiplex(self, records, library=None):
		global TagDefinitions

		try:
			return [TaggedRecord(rawRecord=record,tagDefinitions=TagDefinitions,indexFileParser=self.indexFileParser, indexFileAlias=self.illuminaIndicesAlias, library=library).asFastq(record.sequence, record.plus, record.qual) for record in records]
		except NonMultiplexable:
			raise
		except Exception as e:
			print( e )
