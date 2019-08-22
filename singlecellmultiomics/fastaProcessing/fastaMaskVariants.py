#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
import os
import itertools
import gzip

import pandas as pd
argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""Convert genome FASTA file to convert alternative bases to Ns""")
argparser.add_argument('fasta', metavar='fastafile', type=str, help='Fasta file to mask')
argparser.add_argument('vcf', metavar='vcffile', type=str, nargs='+', help='VCF file(s) to extract variants from')
argparser.add_argument('-o', type=str, default='./masked.fasta.gz', help='output file')
args = argparser.parse_args()

try:
    faIn = pysam.FastaFile(args.fasta)
except ValueError as e:
    print('FAI index missing. Now creating...')
    os.system( f'samtools faidx {args.fasta}' )
except Exception as e:
    raise

try:
    faIn = pysam.FastaFile(args.fasta)
except Exception as e:
    raise

outputHandle = gzip.open(args.o,'wt')
referenceNames = faIn.references
referenceLengths = faIn.lengths
wholeGenome = {
    ref:faIn.fetch(ref) for ref in referenceNames
}
masked = 0
warnedMissing = set()
for rec in itertools.chain(*(pysam.VariantFile(vcfPath) for vcfPath in args.vcf)):

    if not rec.chrom in wholeGenome:
        if not rec.chrom in warnedMissing:
            warnedMissing.add(rec.chrom)
            print(f"WARNING {rec.chrom} is missing in your genome FASTA file. All variants on this chromosome will be skipped!")
        continue
    if rec.start>=(len(wholeGenome[rec.chrom])):
        print(f"WARNING record {rec.chrom} {rec.pos} defines a variant outside the supplied fasta file!")
        continue

    if wholeGenome[rec.chrom][rec.pos] not in [rec.ref, 'N']:
        print(f"WARNING record {rec.chrom} {rec.pos} [{rec.ref}] defines a reference base which doesn't match your fasta file!")
    wholeGenome[rec.chrom][rec.pos] = 'N'
    masked+=1
print(f'Finished conversion. \nMasked {masked} locations')
# WRite to disk
for referenceName in referenceNames:
    print(f'Writing {referenceName}')
    outputHandle.write(f'>{referenceName}\n')
    outputHandle.write(wholeGenome[referenceName])
    outputHandle.write('\n')
outputHandle.close()
