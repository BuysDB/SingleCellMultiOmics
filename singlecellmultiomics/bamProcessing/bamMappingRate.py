#!/usr/bin/env python3
import os
import sys
import glob
import pysam
import colorama

import argparse

f'Please source or use an environment with python 3.6 or higher!'
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
    Obtain mapping rate from libraries
    """)
    argparser.add_argument('bamfiles', nargs='+')
    args = argparser.parse_args()
    for bamfile in args.bamfiles:
        mapped_reads_total = 0
        unmapped_reads_total = 0
        unpaired_reads_total = 0
        total_reads = 0

        try:
            os.system(f'samtools idxstats {bamfile} > {bamfile}.stats.tsv')
            with open(f'{bamfile}.stats.tsv') as f:
                for line in f:
                    chrom, chlen, mapped_reads, unmapped_reads = line.strip().split()
                    mapped_reads = int(mapped_reads)
                    unmapped_reads = int(unmapped_reads)

                    total_reads += unmapped_reads
                    total_reads += mapped_reads

                    if chrom == '*':
                        unmapped_reads_total += unmapped_reads

                    else:
                        mapped_reads_total += mapped_reads
                        unpaired_reads_total += unmapped_reads
            print(f'{bamfile}\t{colorama.Fore.GREEN}{mapped_reads_total}\t{colorama.Fore.RED}{unmapped_reads_total+unpaired_reads_total}\t{colorama.Style.RESET_ALL}{mapped_reads_total/total_reads}')
        except Exception as e:
            continue
