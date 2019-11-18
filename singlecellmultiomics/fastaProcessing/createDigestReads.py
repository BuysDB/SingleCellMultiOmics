#!/usr/bin/env python3
import singlecellmultiomics.utils
import pysam
import collections
import numpy as np
import re
import more_itertools
import gzip
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""In silico digest genome""")
argparser.add_argument('fasta', metavar='fastafile', type=str, help='Fasta file to digest')

argparser.add_argument('-minlen', type=int, default=20, help='minimum read length')
argparser.add_argument('-maxlen', type=int, default=69, help='maximum read length')
args = argparser.parse_args()

r1_read_length = args.maxlen
minlen = args.minlen

ref = pysam.FastaFile(args.fasta)

def seq_to_fastq(seq, header):
    return f'@{header}\n{seq}\n+\n{"E"*len(seq)}\n'

with gzip.open(f'simulated_nlaIII_single_{r1_read_length}.fastq.gz','wt') as out:
    processed = 0


    for contig, contig_len in zip(ref.references, ref.lengths):
        print(f'{contig}\t{contig_len}')
        seq = ref.fetch(contig).upper()
        frag_locations = np.diff( [m.start() for m in re.finditer('CATG', seq)])
        # Generate fragments:
        for start, end_excluding_catg in more_itertools.windowed(  (m.start() for m in re.finditer('CATG', seq)), 2 ):
            if end_excluding_catg is None:
                continue
            end = end_excluding_catg+4


            forward_read = seq[start:end][:r1_read_length]
            if len(forward_read)>=minlen:
                out.write(
                    seq_to_fastq(forward_read, f'DR:fwd;ct:{contig};ds:{start};qn:{len(forward_read)}')
                )
            reverse_read = singlecellmultiomics.utils.reverse_complement(seq[start:end])[:r1_read_length]
            if len(reverse_read)>=minlen:
                out.write(
                    seq_to_fastq(reverse_read, f'DR:rev;ct:{contig};ds:{start};qn:{len(reverse_read)}')
                )

        processed += contig_len
