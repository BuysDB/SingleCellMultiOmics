from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod

# SCCHIC using NLAIII adapter, 384 well format with 3bp UMI followed by "A" base
"""
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
"""
