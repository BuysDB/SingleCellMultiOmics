from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

# NLAIII, 96 well format with 3bp UMI
class NLAIII_96w_c8_u3(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs ):
		self.barcodeFileAlias = 'lennart96NLA'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=3,
		barcodeRead=0, barcodeStart = 3, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs  )
		self.shortName = 'NLAIII96C8U3'
		self.longName = 'NLAIII, 96well CB: 8bp UMI: 3bp'
		self.autoDetectable = True
		self.description = '96 well format. 3bp umi followed by 8bp barcode'

# NLAIII, 384 well format with 3bp UMI
class NLAIII_384w_c8_u3(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs ):
		self.barcodeFileAlias = 'maya_384NLA'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=3,
		barcodeRead=0, barcodeStart = 3, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs  )
		self.shortName = 'NLAIII384C8U3'
		self.longName = 'NLAIII, 384well CB: 8bp UMI: 3bp'
		self.autoDetectable = True
		self.description = '384 well format. 3bp umi followed by 8bp barcode'
