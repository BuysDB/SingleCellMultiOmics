#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pyBigWig
import numpy as np
import pandas as pd


if __name__=='__main__':

    argparser = argparse.ArgumentParser(
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          description='Calculate absolute difference between two bigwig files')

    argparser.add_argument(
        'file_A',
        type=str)

    argparser.add_argument(
        'file_B',
        type=str)

    argparser.add_argument(
        '-o',
        type=str,
        help="Output path (.bw)", required=True)


    args = argparser.parse_args()

    with pyBigWig.open(args.o,'w') as out, \
         pyBigWig.open(args.file_A) as a, \
         pyBigWig.open(args.file_B) as b:

        out.addHeader(list( a.chroms().items()))
        for contig in a.chroms():
            try:
                diff = pd.DataFrame( np.array( a.intervals(contig) )).set_index([0,1] ) - \
                   pd.DataFrame( np.array( b.intervals(contig) )).set_index([0,1] )
            except ValueError:
                continue

            chroms = [contig]*diff.shape[0]
            starts = list([int(i) for i in diff.index.get_level_values(0)] )
            ends = list([int(i) for i in diff.index.get_level_values(1)] )
            values = list(diff.iloc[:,0].values)
            out.addEntries(chroms,starts,ends=ends,values=values)
