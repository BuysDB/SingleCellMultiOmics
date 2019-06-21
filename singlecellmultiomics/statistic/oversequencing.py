#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

class OversequencingHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        if read.has_tag('RC'):
            overseq = read.get_tag('RC')
            self.histogram[overseq]+=1
            if overseq>1:
                self.histogram[overseq-1]-=1
    def __repr__(self):
        return f'The average oversequencing is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'
