#!/usr/bin/env python
# -*- coding: utf-8 -*-
import singlecellmultiomics
from singlecellmultiomics.molecule import MoleculeIterator, NlaIIIMolecule
from singlecellmultiomics.fragment import NlaIIIFragment
import pysam
import collections
import pysamiterators
import pandas as pd
import numpy as np
import re
import more_itertools


def get_random_primer_histogram(
        molecule_source,
        min_mq,
        max_size,
        size_bin_size,
        head=None):
    """
    get_random_primer_histogram
    This method counts the frequencies of the random primers which present in the supplied molecule_source.

    Args:
        molecule_source : iterable yielding Molecules
        min_mq (int) : minimum mean mapping quality of a fragment to be taken into account
        max_size(int) : maximum fragment size
        size_bin_size (int) : bin size in basepairs of the histogram
        head (int) : amount of random sequences to process. When not supplied all random primer sequences in the iterable are counted (expect for the ones with a mapping quality which is too low).

    """
    used = 0
    random_primer_obs = [collections.Counter() for fs in range(
        int(max_size / size_bin_size) + 1)]  # RS

    for i, molecule in enumerate(molecule_source):
        if molecule.get_mean_mapping_qual() < args.min_mq:
            continue

        for frag in molecule:
            if not frag.has_R2():
                continue

            fs = abs(frag.span[2] - frag.span[1])
            size = int(min(max_size - 1, fs) / size_bin_size)
            rp = frag.random_primer_sequence
            if 'N' in rp:
                continue
            random_primer_obs[size][rp] += 1
            used += 1

        if args.head is not None and used > (args.head - 1):
            break

    qf = pd.DataFrame(random_primer_obs).fillna(0).T
    qf[-1] = qf.sum(1)
    qf = qf.T.sort_values(-1, 1).drop(-1)
    qf.index = [bin_index * size_bin_size for bin_index in qf.index]
    return qf


if __name__ == "__main__":
    import argparse
    f'Please use an environment with python 3.6 or higher!'
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
    Known variant locations extraction tool
    """)
    argparser.add_argument('bamfile')
    argparser.add_argument(
        '-max_size',
        type=int,
        default=800,
        help='Maximum fragment size to clip to')
    argparser.add_argument(
        '-size_bin_size',
        help="Bin size in basepairs",
        type=int,
        default=5)
    argparser.add_argument(
        '-head',
        type=int,
        help='Amount of random sequences to count, when not specified all random primers are counted')
    argparser.add_argument('-min_mq', type=int, default=50)
    argparser.add_argument(
        '-o',
        type=str,
        default='./randomer_usage.pickle.gz',
        help='Output pickle/csv path')
    args = argparser.parse_args()

    with pysam.AlignmentFile(args.bamfile) as alignments:
        molecule_source = MoleculeIterator(
            alignments,
            molecule_class=NlaIIIMolecule,
            fragment_class=NlaIIIFragment,
        )

        qf = get_random_primer_histogram(
            molecule_source,
            args.min_mq,
            args.max_size,
            args.size_bin_size,
            head=args.head)

        print('Writing dataframe to disk')
        if args.o.endswith('csv') or args.o.endswith('csv.gz'):
            qf.to_csv(args.o)
        else:
            qf.to_pickle(args.o)
        print('All done')
