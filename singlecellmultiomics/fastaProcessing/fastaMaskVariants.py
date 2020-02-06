#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
import os
import itertools
import gzip
import pandas as pd
import multiprocessing

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Convert genome FASTA file to convert alternative bases to Ns""")
    argparser.add_argument(
        'fasta',
        metavar='fastafile',
        type=str,
        help='Fasta file to mask')
    argparser.add_argument(
        'vcf',
        metavar='vcffile',
        type=str,
        help='VCF file(s) to extract variants from')
    argparser.add_argument(
        '-o',
        type=str,
        default='./masked.fasta.gz',
        help='output file')
    argparser.add_argument(
        '-t',
        type=int,
        default=10,
        help='Amount of contigs to process in parallel')
    args = argparser.parse_args()

    try:
        faIn = pysam.FastaFile(args.fasta)
    except ValueError as e:
        print('FAI index missing. Now creating...')
        os.system(f'samtools faidx {args.fasta}')
    except Exception as e:
        raise

    try:
        faIn = pysam.FastaFile(args.fasta)
    except Exception as e:
        raise

    if args.o.endswith('.gz'):
        outputHandle = gzip.open(args.o, 'wb')
    else:
        outputHandle = open(args.o, 'wb')
    referenceNames = faIn.references
    referenceLengths = dict(zip(referenceNames, faIn.lengths))
    totalMasked = 0

    def get_masked_bytearray( jargs ):

        totalMasked = 0
        (chrom, fasta_file_path, variant_file_path) = jargs
        with pysam.FastaFile(fasta_file_path) as faIn, pysam.VariantFile(variant_file_path, threads=4) as  variants:

            chrom_seq = bytearray(faIn.fetch(chrom), 'ascii')
            if chrom in variants.index:
                print(f'masking {chrom}')
                for rec in variants.fetch(chrom):
                    if rec.start >= (referenceLengths[rec.chrom]):
                        print(
                            f"WARNING: record {rec.chrom} {rec.pos} defines a variant outside the supplied fasta file!")
                        continue

                    if len(rec.alleles) == 1 and len(rec.alleles[0]) == 1:
                        chrom_seq[rec.pos-1] = 78  # ord N
                        totalMasked += 1
                    elif len(rec.alleles[0]) == 1 and len(rec.alleles[1]) == 1:

                        chrom_seq[rec.pos-1] = 78  # ord N
                        totalMasked += 1
            print(f'Masked {totalMasked} bases of {chrom}')
            return chrom_seq


    workers = multiprocessing.Pool(args.t)


    for chrom, chrom_seq in zip(
                referenceNames,
                workers.imap(
                    get_masked_bytearray,
                    ((chrom, args.fasta, args.vcf )
                    for chrom in referenceNames ))):
        # Write chromsome
        outputHandle.write(f'>{chrom}\n'.encode('ascii'))
        outputHandle.write(chrom_seq)
        outputHandle.write('\n'.encode('ascii'))

    outputHandle.close()
