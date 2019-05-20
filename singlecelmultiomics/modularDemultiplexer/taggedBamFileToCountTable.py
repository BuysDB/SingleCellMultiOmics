#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import baseDemultiplexMethods
TagDefinitions = baseDemultiplexMethods.TagDefinitions

argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Convert BAM file(s) which has been annotated using featureCounts or otherwise to a count matrix')
argparser.add_argument('-o',  type=str, help="output csv path, or pandas dataframe if path ends with pickle.gz", required=False)
argparser.add_argument('-featureTags',  type=str, default=None, help='Tag(s) used for defining the COLUMNS of the matrix. Single dimension.')
argparser.add_argument('-joinedFeatureTags',  type=str, default=None, help='These define the COLUMNS of your matrix. For example if you want allele (DA) and restriction site (DS) use DA,DS. If you want a column containing the chromosome mapped to use "chrom" as feature. Use this argument if you want to use a multidimensional index')
argparser.add_argument('-sampleTags',  type=str, default='SM', help='Comma separated tags defining the names of the ROWS in the output matrix')
argparser.add_argument('alignmentfiles',  type=str, nargs='*')
argparser.add_argument('-head',  type=int, help='Run the algorithm only on the first N reads to check if the result looks like what you expect.')
argparser.add_argument('--divideMultimapping', action='store_true', help='Divide multimapping reads over all targets. Requires the XA or NH tag to be set.')

argparser.add_argument('-offset', type=int, default=0, help="Add offset to bin. If bin=1000, offset=200, f=1999 -> 1200. f=4199 -> 3200")
argparser.add_argument('-bin', type=int, help="Devide and floor to bin features. If bin=1000, f=1999 -> 1000." )
argparser.add_argument('--showBinEnd', action='store_true', help="If True, then show DS column as 120000-220000, otherwise 120000 only. This specifies the bin range in which the read was counted" )
argparser.add_argument('-binTag',default='DS' )

argparser.add_argument('--dedup', action='store_true', help='Count only the first occurence of a molecule. Requires RC tag to be set. Reads without RC tag will be ignored!')
argparser.add_argument('--showtags',action='store_true', help='Show a list of commonly used tags' )

args = argparser.parse_args()

if args.o is None and  args.alignmentfiles is not None:
    args.showtags=True

if args.showtags:
    # Find which tags are available in the file:
    head = 1000
    tagObs = collections.Counter()
    for bamFile in args.alignmentfiles:
        with pysam.AlignmentFile(bamFile) as f:
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
if args.o is None:
    raise ValueError('Supply an output file')
if args.alignmentfiles is None:
    raise ValueError('Supply alignment (BAM) files')

joinFeatures = False
if args.featureTags is not None:
    featureTags= args.featureTags.split(',')

if args.joinedFeatureTags is not None:
    featureTags= args.joinedFeatureTags.split(',')
    joinFeatures=True

sampleTags= args.sampleTags.split(',')
countTable = collections.defaultdict(collections.Counter) # cell->feature->count

def readTag(read, tag, missing='Missing', defective='Defective'):
    value=None
    if tag=='chrom':
        return str(read.reference_name)
    if not read.has_tag(tag):
        return missing
    else:
        try:
            value = read.get_tag(tag)
        except Exception as e:
            value = defective
    return value



for bamFile in args.alignmentfiles:
    assigned = 0
    with pysam.AlignmentFile(bamFile) as f:

        for i,read in enumerate(f):
            if i%1_000_000==0:
                print(f"{bamFile} Processed {i} reads, assigned {assigned}, completion:{100*(i/(0.001+f.mapped+f.unmapped+f.nocoordinate))}%")
            if read.is_unmapped or args.dedup and ( not read.has_tag('RC') or (read.has_tag('RC') and read.get_tag('RC')!=1)):
                continue
            sample =tuple( readTag(read,tag) for tag in sampleTags )

            countToAdd = (0.5 if read.is_paired else 1)
            assigned+=1
            if args.divideMultimapping:
                    if read.has_tag('XA'):
                        countToAdd = countToAdd/len( read.get_tag('XA').split(';') )
                    elif read.has_tag('NH'):
                        countToAdd = countToAdd/int(read.get_tag('NH') )
                    else:
                        countToAdd = countToAdd

            if not joinFeatures:
                for tag in featureTags:
                    if args.bin is not None and tag == args.binTag and read.has_tag(tag):
                        # feature = str( int( readTag(read,tag)/args.bin))
                        feat = str( int( readTag(read,tag) ) / (args.bin + args.offset) * (args.bin) + args.offset)
                        if args.showBinEnd:
                            binend = int(feat) + args.bin
                            feat = ''.join[feat, "-", binend]
                    else:
                        feat = str(readTag(read,tag))
                    countTable[sample][feat]+=countToAdd
            else:
                feature =[]
                for tag in featureTags:
                    if args.bin is not None and tag == args.binTag and read.has_tag(tag):
                        feat = str( int( readTag(read,tag)/ (args.bin + args.offset) ) * (args.bin) + args.offset)
                        if args.showBinEnd:
                            binend = str(int(feat) + args.bin)
                            feat = ''.join([feat, "-", binend])
                    else:
                        feat =str(readTag(read,tag))
                    feature.append( feat )
                countTable[sample][tuple(feature)]+=countToAdd
            if args.head is not None and i>args.head:
                break

def tagToHumanName(tag,TagDefinitions ):
    if not tag in TagDefinitions:
        return tag
    return TagDefinitions[tag].humanName

print("Finished counting, now exporting...")
df = pd.DataFrame.from_dict( countTable )
df.columns.set_names([tagToHumanName(t,TagDefinitions ) for t in sampleTags], inplace=True)
if joinFeatures:
    df.index.set_names([tagToHumanName(t, TagDefinitions) for t in featureTags], inplace=True)
else:
    df.index.set_names(','.join([tagToHumanName(t, TagDefinitions) for t in featureTags]), inplace=True)

if args.o.endswith('.pickle') or args.o.endswith('.pickle.gz'):
    df.to_pickle(args.o)
else:
    df.to_csv(args.o)
print("Finished export.")
