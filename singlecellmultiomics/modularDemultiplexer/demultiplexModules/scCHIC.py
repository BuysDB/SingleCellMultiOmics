from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod,NonMultiplexable

# SCCHIC using NLAIII adapter, 384 well format with 3bp UMI followed by "A" base

class SCCHIC_384w_c8_u3(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs ):
		self.barcodeFileAlias = 'maya_384NLA'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=3,
		barcodeRead=0, barcodeStart = 3, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs  )
		self.shortName = 'SCHIC384C8U3'
		self.longName = 'Single cell CHIC, 384well CB: 8bp UMI: 3bp'
		self.autoDetectable = True
		self.description = '384 well format. 3bp umi followed by 8bp barcode and a single A'

		self.sequenceCapture[0] = slice( self.barcodeLength+ self.umiLength + 1, None) # dont capture the first base


	def demultiplex(self, records, **kwargs):

		if kwargs.get('probe') and records[0].sequence[self.barcodeLength+ self.umiLength+1]!='A':
			raise NonMultiplexable

		taggedRecords = UmiBarcodeDemuxMethod.demultiplex(self,records, **kwargs)
		# add first 2 bases as ligation tag:
		taggedRecords[0].addTagByTag('lh', records[0].sequence[self.barcodeLength+ self.umiLength+2], isPhred=False)
		taggedRecords[0].addTagByTag('lq', records[0].qual[self.barcodeLength+ self.umiLength+2], isPhred=True)
		taggedRecords[1].addTagByTag('lh', records[0].sequence[self.barcodeLength+ self.umiLength+2], isPhred=False)
		taggedRecords[1].addTagByTag('lq', records[0].qual[self.barcodeLength+ self.umiLength+2], isPhred=True)
		taggedRecords[0].sequence = taggedRecords[0].sequence[1:]
		taggedRecords[0].qualities = taggedRecords[0].qualities[1:]
		return taggedRecords
