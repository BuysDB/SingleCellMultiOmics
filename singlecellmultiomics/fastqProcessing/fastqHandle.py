#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gzip
from singlecellmultiomics.pyutils.handlelimiter import HandleLimiter

class FastqHandle:
	def __init__(self, path, pairedEnd=False, single_cell=False, maxHandles=500 ):
		self.pe = pairedEnd
		self.sc = single_cell
		self.path = path
		if not self.sc:
			if pairedEnd:
				self.handles = [ gzip.open(path+'R1.fastq.gz', 'wt'),   gzip.open(path+'R2.fastq.gz', 'wt') ]
			else:
				self.handles = [ gzip.open(path+'reads.fastq.gz', 'wt') ]
		else:

			self.handles = HandleLimiter(compressionLevel=1, maxHandles=maxHandles)

	def write(self, records ):
		if self.sc:
			for readIdx, record in zip(('R1','R2'), records):
				# Obtain cell from record:
				cell = f"{record.tags.get('BI','no_cell_id')}.{record.tags.get('MX','unk')}"
				self.handles.write(f'{self.path}.{cell}.{readIdx}.fastq.gz',str(record), method=1)
		else:
			for handle, record in zip(self.handles, records):
				handle.write(str(record))
	def close(self):
		if self.sc:
			self.handles.close()
		else:
			for handle in self.handles:
				handle.close()
