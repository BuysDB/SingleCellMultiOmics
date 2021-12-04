#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyBigWig
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import get_binned_counts
import pandas as pd
import argparse
from singlecellmultiomics.bamProcessing import get_contig_sizes
import numpy as np
from singlecellmultiomics.utils.path import get_valid_filename

def bam_to_wig(bam_paths, write_path, bin_size, method='sum', verbose=False, n_threads=None, sample_mapping_function=None):
    if verbose:
        print('Counting ...')
    counts = get_binned_counts(bam_paths, bin_size=bin_size, n_threads=n_threads)

    if verbose:
        print('Writing ...')

    if sample_mapping_function is not None:

        handles = {}
        with pysam.AlignmentFile(bam_paths[0]) as alignments:
            cs = get_contig_sizes(alignments)

            targets =  dict( zip(counts.columns, map(sample_mapping_function, counts.columns) ))
            for target in set(targets.values()):
                if target is None:
                    continue
                # Select only cells which have the current target label:
                subset = counts[ [cell for cell,t in targets.items() if t==target] ]
                # And write:
                values = subset.sum(1).sort_index()

                # Write values
                with pyBigWig.open(write_path.replace('.bw',f'_{str(target)}.bw'),'w') as out:
                    out.addHeader(list(zip(alignments.references, alignments.lengths)))

                    for contig in alignments.references:

                        if contig not in values.index.get_level_values(0):
                            continue

                        print(f'Writing data for {contig}, for {target}')
                        v = values.loc[[contig]]
                        out.addEntries(
                            list(v.index.get_level_values(0).values), #Contig
                            list(v.index.get_level_values(1).values), #Start
                            ends= list( np.clip(  (v.index.get_level_values(1)+bin_size) .values, 0, cs[contig]-1) ) ,  #end
                            values= np.array(v.values, dtype=np.float32))

    else:
        with pysam.AlignmentFile(bam_paths[0]) as alignments, pyBigWig.open(write_path,'w') as out:

            cs = get_contig_sizes(alignments)
            # Write header
            out.addHeader(list(zip(alignments.references, alignments.lengths)))
            values = counts.sum(1).sort_index()
            print(values)
            # Write values
            for contig in alignments.references:

                if contig not in values.index.get_level_values(0):
                    continue
                print(f'Writing data for {contig}')
                v = values.loc[[contig]]
                out.addEntries(
                    list(v.index.get_level_values(0).values), #Contig
                    list(v.index.get_level_values(1).values), #Start
                    ends= list( np.clip(  (v.index.get_level_values(1)+bin_size) .values, 0, cs[contig]-1) ) ,  #end
                    values= np.array(v.values, dtype=np.float32))


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Bam file to bigwig. Counts based on DS tag by default, otherwise falls back to the start coordinate of R1. Does not qcfail flagged reads and duplicate flagged reads.')

    argparser.add_argument('alignmentfiles', type=str, nargs='+')
    argparser.add_argument('-o', type=str, required=True, help='Output path (.bw)')
    argparser.add_argument('-t', type=int,  help='Amount of threads. Uses all available if not supplied')

    argparser.add_argument('-bin_size', type=int, required=True)

    pseudobulk_gr = argparser.add_argument_group('Pseudobulk settings')
    pseudobulk_gr.add_argument(
        '-pseudobulk_SM_csv',
        type=str,
        help="""Path to a CSV file which contains for every barcode index (SM tag) to what group it belongs.
         The CSV file has no header and two columns, the first column contains the sample name,
        the second the target sample name. Multiple barcode indices can share the same sample name, this will create a pseudobulk signal.
        Expects a comma as delimiter."""
        )

    args = argparser.parse_args()
    assert args.o.endswith('.bw')

    sample_mapping_function = None
    if args.pseudobulk_SM_csv is not None:
        sm_sample_map = {str(sm):get_valid_filename(str(sample))
            for sm, sample in pd.read_csv(args.pseudobulk_SM_csv,header=None,index_col=0).iloc[:,0].to_dict().items() }
        def sample_mapping_function(s):
            return sm_sample_map.get(s)


    bam_to_wig(args.alignmentfiles, args.o, args.bin_size, verbose=True, sample_mapping_function=sample_mapping_function,n_threads=args.t)
