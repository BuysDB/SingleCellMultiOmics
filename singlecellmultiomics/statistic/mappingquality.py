
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

class MappingQualityHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        self.histogram[read.mapping_quality]+=1
    def __repr__(self):
        return f'The average mapping quality is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'
