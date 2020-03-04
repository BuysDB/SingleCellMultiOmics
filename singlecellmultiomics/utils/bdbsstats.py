# Buys's statsistics doodles
import scipy.stats
import numpy as np
import collections
import sklearn.metrics
import matplotlib.pyplot as plt

def sampleFromCounter(counter, n, replace=True, outputCounter=False):
	counterSum = sum(counter.values())
	if replace:
		samples = np.random.choice(list(counter.keys()), n, p=[ float(counter[className])/float(counterSum) for className in counter ])
		if outputCounter:
			#@todo optimize
			collections.Counter( samples )
		else:
			return(  samples )
	if counterSum<n:
		raise ValueError('Cannot pick %s samples from a list of %s samples' % (n, counterSum))

	samples = []

	if outputCounter:
		samples=collections.Counter()
	for i in range(0,n):
		counterSum = sum(counter.values())
		classProbs = [ float(counter[className])/float(counterSum) for className in counter ]
		s = list( np.random.choice( list(counter.keys()), 1, p=classProbs) )[0]
		if outputCounter:
			samples[s]+=1
		else:
			samples.append( s )
		counter[s]-=1
	return(samples)


class IterativeHistogram(object):

	def __init__(self, minV, maxV, bins):
		self.binCount = bins
		self.min = minV
		self.max = maxV
		self.counts = np.zeros( self.binCount )

		self.binStarts = np.linspace(self.min, self.max, self.binCount,endpoint=False)
		self.binSize = self.binStarts[1] - self.binStarts[0]

	def add(self, vector):
		try:
			locs = np.searchsorted(self.binStarts, vector,'right')-1
			#print('input: %s' %vector)
			#print('point: %s'%locs)
			#print('starts:%s' % self.binStarts)

			np.add.at( self.counts, locs, 1)
			print(self.counts)
		except:

			exit('Failed adding vector to histogram')
		return(self)

	def plot(self,path=None):
		fig, ax = plt.subplots() #figsize=(120, 10))
		for i in range(0,self.binCount):
			plt.bar(self.binStarts[i],self.counts[i], float(self.max-self.min)/float(self.binCount))

		plt.ylabel("Density")
		plt.xlabel("PIC score")
		plt.title('PIC score distribution')
		ax.legend(loc='upper right')
		plt.yscale('log', nonposy='clip')
		if path==None:
			plt.show()
		else:
			plt.savefig(path)

	def __add__(self, histogramToAdd):
		n = IterativeHistogram(self.min, self.max, self.binCount)
		n.counts += self.counts
		n.counts += histogramToAdd.counts
		return(n)

	def p(self, value):
		return( sum( self.counts[ np.searchsorted(self.binStarts, value,'right'): ] )/sum(self.counts) )


	def getScipyRandomVariable(self):

		values = self.binStarts
		probabilities = self.counts/np.sum(self.counts)
		X = scipy.stats.rv_discrete(values=(values, probabilities))
		return(X)


	def __iter__(self):
		self.iterIndex = 0
		return(self)

	def __next__(self):
		value = self.min
		index = 0
		for i in range(0,self.binCount):
			index+=self.counts[i]

			if index> self.iterIndex:
				self.iterIndex+=1
				return( self.binStarts[i]  ) #+ 0.5*self.binSize
		raise StopIteration

def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = sklearn.metrics.mutual_info_score(None, None, contingency=c_xy)
    return mi

## Example of Chi-squared test between two distributions:
#
#
# observedValues = [16, 18, 16, 14, 12, 12]
# expectedValues = [16, 16, 16, 16, 16, 8]
# binCount = 20
#
#
# expectedDistribution = list( IterativeHistogram(0,20, binCount).add(expectedValues) )
# print(expectedDistribution)
#
# observedDistribution = list( IterativeHistogram(0,20, binCount).add(observedValues) )
#
# print( scipy.stats.chisquare( observedDistribution, f_exp=expectedDistribution).pvalue )
#
# print( scipy.stats.ks_2samp( observedDistribution, expectedDistribution ))
# expectedDistribution = IterativeHistogram(0,20, binCount).add(expectedValues).getScipyRandomVariable()
# observedDistribution = IterativeHistogram(0,20, binCount).add(observedValues).getScipyRandomVariable()
#
# print( scipy.stats.kstest( observedDistribution, expectedDistribution ))
