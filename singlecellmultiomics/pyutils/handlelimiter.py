#Handle limiter written by Buys de Barbanson, Hubrecht 2017
#This class allows for writing many files at the same time without the hassle of thinking about handle limitations.
import gzip
import time
class HandleLimiter(object):

	def __init__(self, maxHandles = 32,pruneEvery=10000, compressionLevel=1):
		self.openHandles = {}
		self.seen = set() #Which files have been opened before
		self.maxHandles = maxHandles
		self.pruneEvery = pruneEvery
		self.pruneIntervalCounter = 0
		self.compressionLevel = compressionLevel

	def write(self, path, string, method=None, forceAppend=False): #0= plain, 1:gzip

		if not path in self.openHandles:

			self.openHandles[path]={}
			failedOpening= True
			while failedOpening:
				try:
					if path in self.seen or forceAppend:
						#Append when we already wrote to the file
						if method==1:
							self.openHandles[path]['handle'] = gzip.open(path,'ab',self.compressionLevel)
						else:
							self.openHandles[path]['handle'] = open(path,'a')
					else:
						#Open as new file when it is the first write
						if method==1:
							self.openHandles[path]['handle'] = gzip.open(path,'wb',self.compressionLevel)
						else:
							self.openHandles[path]['handle'] = open(path,'w')
						self.seen.add(path) # Remember that we accessed this file
					failedOpening= False
				except Exception as e:
					failedOpening = True
					if len(self.openHandles)>1:
						self.close()
					else:
						#This is mayorly bad...
						print('Failed writing to %s, even after closing all other open file-handles. Out of options...' % path)
						print(e)
						#Raise the error to the parent method
						raise
		if method==0:
			self.openHandles[path]['handle'].write(string)
		else:
			self.openHandles[path]['handle'].write(bytes(string, 'UTF-8'))
		self.openHandles[path]['lastw'] = time.time()
		self.pruneIntervalCounter+=1
		if self.pruneIntervalCounter>=self.pruneEvery:
			self.prune()

	def prune(self):
		if len(self.openHandles)>self.maxHandles:
			toPrune = len(self.openHandles)-self.maxHandles
			pathsToPrune = sorted( self.openHandles.keys(), key=lambda path: (self.openHandles[path]['lastw']) )[:toPrune] #,reverse=True
			for path in pathsToPrune:
				if 'handle' in self.openHandles[path]:
					try:
						self.openHandles[path]['handle'].close()
					except Exception as e:
						pass

				self.openHandles.pop(path)
		self.pruneIntervalCounter=0

	def close(self):

		k = self.openHandles.keys()
		destroyed = []
		for path in k:
			if 'handle' in self.openHandles[path]:
				try:
					self.openHandles[path]['handle'].close()
				except:
					pass
			else:
				print('Closed broken file handle for %s' % path)
			destroyed.append(path)
		for delete in destroyed:
			self.openHandles.pop(delete)
