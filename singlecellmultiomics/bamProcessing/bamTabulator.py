#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import singlecellmultiomics
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions
if __name__=='__main__':
    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Tabulate a bam file to a file where every line corresponds to the features of a single read')
    argparser.add_argument('-o',  type=str, help="output csv path", required=False)
    argparser.add_argument('alignmentfile',  type=str)
    argparser.add_argument('featureTags',  type=str)
    argparser.add_argument('-head',  type=int, help='Run the algorithm only on the first N reads to check if the result looks like what you expect.')
    argparser.add_argument('-e',  action='store_true', help='Skip reads where any attribute is missing')
    argparser.add_argument('--dedup', action='store_true', help='Count only the first occurence of a molecule. Requires RC tag to be set. Reads without RC tag will be ignored!')
    argparser.add_argument('--showtags',action='store_true', help='Show a list of commonly used tags, and tags present in your bam file' )
    args = argparser.parse_args()

    if args.featureTags is None:
        print('Please Supply feature tags. Select from:')
        args.o = None
        args.featureTags = None
        args.showtags=True

    if args.showtags:
        # Find which tags are available in the file:
        head = 1000
        tagObs = collections.Counter()

        with pysam.AlignmentFile(args.alignmentfile, ignore_truncation=True ) as f:
            for i,read in enumerate(f):
                tagObs += collections.Counter([ k for k,v in   read.get_tags(with_value_type=False)] )
                if i==(head-1):
                    break
        import colorama

        print(f'{colorama.Style.BRIGHT}Tags seen in the supplied bam file(s):{colorama.Style.RESET_ALL}')
        for observedTag in tagObs:
            tag = observedTag
            if observedTag in TagDefinitions:
                t = TagDefinitions[observedTag]
                humanName = t.humanName
                isPhred = t.isPhred
            else:
                t = observedTag
                isPhred = False
                humanName=f'{colorama.Style.RESET_ALL}<No information available>'

            allReadsHaveTag = ( tagObs[tag]==head )
            color = colorama.Fore.GREEN if allReadsHaveTag else colorama.Fore.YELLOW
            print(f'{color}{ colorama.Style.BRIGHT}{tag}{colorama.Style.RESET_ALL}\t{color+humanName}\t{colorama.Style.DIM}{"PHRED" if isPhred else ""}{colorama.Style.RESET_ALL}' + f'\t{"" if allReadsHaveTag else "Not all reads have this tag"}')

        print(f'{colorama.Style.BRIGHT}\nAVO lab tag list:{colorama.Style.RESET_ALL}')
        for tag,t in TagDefinitions.items():
            print(f'{colorama.Style.BRIGHT}{tag}{colorama.Style.RESET_ALL}\t{t.humanName}\t{colorama.Style.DIM}{"PHRED" if t.isPhred else ""}{colorama.Style.RESET_ALL}')
        exit()

    if args.alignmentfile is None:
        raise ValueError('Supply alignment (BAM) files')



    featureTags= args.featureTags.split(',')
    countTable = collections.defaultdict(collections.Counter) # cell->feature->count
    def tagToHumanName(tag,TagDefinitions ):
        if not tag in TagDefinitions:
            return tag
        return TagDefinitions[tag].humanName

    if args.o is not None:
        o = args.o
        tf = open(args.o,'w')
        # Write header:
        tf.write( '\t'.join([tagToHumanName(t, TagDefinitions) for t in featureTags])+'\n' )
    wrote=0
    try:
        with pysam.AlignmentFile(args.alignmentfile, ignore_truncation=True  ) as f:
            for i,read in enumerate(f):
                if args.dedup and read.is_duplicate:
                    continue
                values = [
                    str(singlecellmultiomics.modularDemultiplexer.metaFromRead(read,tag))
                    for tag in featureTags
                ]
                if args.e and any( (v is None for v in values )):
                    continue
                line = '%s\n' % '\t'.join(values)
                if args.o is None:
                    print(line,end="")
                else:
                    tf.write( line )
                wrote+=1
                if args.head and wrote>args.head:
                    break
    except (KeyboardInterrupt,BrokenPipeError) as e:
        pass

    if args.o is not None:
        tf.close()
