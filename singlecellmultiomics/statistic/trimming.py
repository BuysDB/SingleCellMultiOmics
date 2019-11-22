#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils


class TrimmingStats(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.totalFragmentsTrimmed = 0

    def processRead(self, read):
        if read.has_tag('a1') or read.has_tag(
                'eB') or read.has_tag('A2') or read.has_tag('EB'):
            self.totalFragmentsTrimmed += 1

    def __repr__(self):
        return f'Trimmed fragments: {self.totalFragmentsTrimmed}'

    def __iter__(self):
        yield 'Trimmed fragments', self.totalFragmentsTrimmed
