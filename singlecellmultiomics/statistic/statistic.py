#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import collections
import pandas as pd
import numpy as np
import singlecellmultiomics.pyutils as pyutils


class Statistic(object):

    """
    Statistic object, initialised with arguments

    Parameters
    ----------
    args : argparse object

    """

    def __init__(self, args):
        self.args = args

    def processRead(self, R1,R2=None):
        """
        Update the statistic with information from READ

        Parameters
        ----------
        read : PySAM Aligned segment

        Returns
        ----------
        None
        """
        pass

    def __repr__(self):
        return 'dummy'


class StatisticHistogram(Statistic):
    def __init__(self, args):
        Statistic.__init__(self, args)
        self.histogram = collections.Counter()

    def __repr__(self):
        return f'Mean {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'

    def __iter__(self):
        return iter(self.histogram.most_common())

    def to_csv(self, path):
        pd.DataFrame({__class__.__name__: self.histogram}).to_csv(path)
