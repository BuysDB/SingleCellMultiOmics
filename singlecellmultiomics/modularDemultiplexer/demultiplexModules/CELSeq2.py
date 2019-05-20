from baseDemultiplexMethods import UmiBarcodeDemuxMethod
import barcodeFileParser

# Cell seq 2 with 6bp UMI
class CELSeq2_c8_u6(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs  ):
		self.barcodeFileAlias = 'celseq2'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=6,
		barcodeRead=0, barcodeStart = 6, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs   )
		self.shortName = 'CS2C8U6'
		self.longName = 'CELSeq 2, CB: 8bp, UMI: 6bp'
		self.autoDetectable = True
		self.description = 'R1 starts with a 6bp UMI  followed by a 8bp cell barcode'


# Cell seq 2 with 8bp UMI
class CELSeq2_c8_u8(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs  ):
		self.barcodeFileAlias = 'celseq2'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=8,
		barcodeRead=0, barcodeStart = 8, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs   )
		self.shortName = 'CS2C8U8'
		self.longName = 'CELSeq 2, CB: 8bp, UMI: 8bp'
		self.autoDetectable = True
		self.description = 'R1 starts with a longer 8bp UMI  followed by a 8bp cell barcode'

class CELSeq2_c8_u8_NNLAIII(UmiBarcodeDemuxMethod):
	def __init__(self, barcodeFileParser,  **kwargs  ):
		self.barcodeFileAlias = 'celseq2_noNla'
		UmiBarcodeDemuxMethod.__init__(self,
		umiRead=0, umiStart = 0, umiLength=8,
		barcodeRead=0, barcodeStart = 8, barcodeLength=8,
	 	barcodeFileAlias = self.barcodeFileAlias ,barcodeFileParser=barcodeFileParser,  **kwargs   )
		self.shortName = 'CS2C8U8NNLA'
		self.longName = 'CELSeq 2, CB: 8bp, UMI: 8bp, NLAIII free'
		self.autoDetectable = True
		self.description = 'CEL-Seq2 without NLAIII digestable barcodes '
