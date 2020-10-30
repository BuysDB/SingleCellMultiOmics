#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:04:46 2020

@author: rvanesch
"""
#%%
import os
import pysam
import argparse
#import sys

#%%
if __name__ == '__main__':
    
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Filter bam file for containing less than an indicated ratio
        of soft clips (-sc)
    """)
    argparser.add_argument(
        'bamfile',
        metavar='bamfile',
        type=str,
        help='input bam file, make sure bam index (.bai) is located in same folder')
    argparser.add_argument(
        '-sc',
        type=float,
        default=0.3,
        help='maximum ratio of soft-clips to total reads (soft-clips + matched reads)')
    argparser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output bam path')
    argparser.add_argument(
        '-f',
        action='store_true',
        help='Force overwrite of existing files')
    argparser.add_argument(
        'l',
        action='store_true',
        help='flag to output logfile')
    #argparser.add_argument('-head', type=int)
    args = argparser.parse_args()

    if os.path.exists(args.o) and not args.f:
        raise ValueError(
            f'The output file {args.o} already exists! Use --f or remove the file')
    if (args.sc > 1) or (args.sc < 0):
        raise ValueError(
            'sc is a ratio and should be a float between 0 and 1')
    bamFile = pysam.AlignmentFile(args.bamfile, "rb")
    header = bamFile.header.copy()
    outputFile = pysam.AlignmentFile(args.o, "wb", header=header)
    if args.l:
        log = open('sc_filter_log.txt', 'w')
    
    max_sc = args.sc
    sc_pass_total = 0
    sc_fail_total = 0
    sc_error = 0
    for r in bamFile:
        sc_count = 0
        m_count= 0
        #r.cigar returns a list of tuples describing the cigar, with tuple[0] indicating the type and 
        #tuple[1] containing the length.
        #matches are indicated by tuple[0] == 0, softclips by tuple[0] == 4 
        try: 
            for tup in r.cigartuples:
                if tup[0]==0:
                    m_count += tup[1]
                elif tup[0]==4:
                    sc_count += tup[1]

            sc_ratio = (sc_count/(sc_count+m_count))
        except TypeError:
            sc_ratio = 1.1
            if args.l:
                try:
                    log.write("typeError at: "+str(r.reference_name)+"\n")
                    sc_error += 1
                except TypeError:
                    log.write("no read in read\n")
                    sc_error += 1

            
        if sc_ratio < max_sc:
            outputFile.write(r)
            sc_pass_total+=1
        else:
            sc_fail_total+=1
    
    if args.l:
        log.write("filtering done \n")
        log.write("bamfile: "+args.bamfile+"\n")
        log.write("soft_clip filter value: "+str(max_sc)+"\n")
        log.write("passed: "+str(sc_pass_total)+" reads \n")
        log.write("failed: "+str(sc_fail_total)+" reads \n")
        log.write("errors: "+str(sc_error)+" reads \n")
        log.write("output file: "+args.o)
        log.close()
    
    outputFile.close()
