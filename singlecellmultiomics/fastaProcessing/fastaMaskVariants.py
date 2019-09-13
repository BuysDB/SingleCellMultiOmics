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
argparser.add_argument('vcf', metavar='vcffile', type=str, help='VCF file(s) to extract variants from')
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

outputHandle = gzip.open(args.o,'wb')
referenceNames = faIn.references
referenceLengths = dict( zip(referenceNames, faIn.lengths))

#wholeGenome = {
#    ref:faIn.fetch(ref) for ref in referenceNames
#}
warnedMissing = set()
totalMasked=0
variants = pysam.VariantFile(args.vcf)
for chrom in referenceNames:
    try:
        chrom_seq = bytearray(faIn.fetch(chrom),'ascii' )
        for rec in variants.fetch(chrom):
            if rec.start>=(referenceLengths[rec.chrom]):
                print(f"WARNING record {rec.chrom} {rec.pos} defines a variant outside the supplied fasta file!")
                continue

            if len(rec.alleles)==1 and len(rec.alleles[0])==1:
                chrom_seq[rec.pos-1] = 78 #ord N
                totalMasked+=1
            elif len(rec.alleles[0])==1 and len(rec.alleles[1])==1:

                chrom_seq[rec.pos-1] = 78 #ord N
                totalMasked+=1

    except ValueError:
        print(f"No variants for {chrom}")
        pass
    # Write chromsome
    outputHandle.write(f'>{chrom}\n'.encode('ascii'))
    outputHandle.write(chrom_seq)
    outputHandle.write('\n'.encode('ascii'))

outputHandle.close()
