#!/usr/bin/env python
'''
DESCRIPTION

    Split bams into pseudobulks based on metadata

FOR HELP

    python split_bam_by_cluster.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-01-09
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv
import os

def reformat_header(header, add_prefix="chr"):
    new_header = header.to_dict()
    contigs = new_header['SQ']
    new_contigs = []
    for contig in contigs:
        contig['SN'] = ''.join([add_prefix, contig['SN']])
        new_contigs.append(contig)
    new_header['SQ'] = new_contigs
    return pysam.AlignmentHeader.from_dict(new_header)

def scan_bam_get_sampnames(infile, tag="SM"):
    '''
    Read bam file and get sampnames
    '''
    tags = set()
    with pysam.AlignmentFile(infile, "rb") as inbam:
        for read in inbam:
            tags.add(readTag(read, tag))
    return(tags)

def readTag(read, tag, missing='Missing', defective='Defective'):
    '''From /hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/taggedBamFileToCountTable.py'''
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

def main():
    parser = argparse.ArgumentParser(description='Split bams into pseudobulks based on metadata')
    parser.add_argument('-infile', metavar='INFILE', required=True,
                        help='Merged bam file')
    parser.add_argument('-annotfile', metavar='INFILE', required=True,
                        help='Annotations of merged clusters. Expect first row are colum names (skipped), first column are cell names (should match SM tag in bam) and second column is cluster name, used for bam file outputs')
    parser.add_argument('-outdir', metavar='OUTDIR', required=True,
                        help='Output directory')
    parser.add_argument('-tagid', metavar='TAG NAME', required=False, default="SM", 
                        help='Tag name to match first column in annot file (default is SM)')
    parser.add_argument('--annot_no_colnames', action='store_true', help='Set if annot has no column name')
    parser.add_argument('-mapq', metavar="MAPQ value", required=True, default=40, type=int)
    parser.add_argument('--add_chr_prefix', action='store_true', help="Add chr prefix to chromosome name")
    parser.add_argument('--overwrite', action='store_true', help="Does not check if bam files exists, will just write to file")
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('-logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    samp_tag_id = args.tagid
    samp_cname_indx = 0
    samp_cluster_indx = 1

    # initialize outbam objects
    bname = os.path.splitext(os.path.basename(args.infile))[0]
    writebamdic = {}  # bam write objs
    sorteddic={}  # output files for sorted bam files
    unsorteddic={}  # tmp files unsorted, will be deleted afterwards
    

    # # get sampnames_all by reading the bam file
    # sampnames_all = scan_bam_get_sampnames(args.infile, tag = samp_tag_id)  # sampnames_all are sample names
    # print(str(len(sampnames_all)) + " sampnames_all found in bam file")
    # # samp_ignored = set()

    # read annot file
    print("Reading annot file...")
    samp2clstr_dic = {}
    samps_good = []
    clstrs_lst = []
    clstrs = set()
    with open(args.annotfile, "r") as annotfile:
        jreader = csv.reader(annotfile, delimiter = "\t")
        if not args.annot_no_colnames:
            jheader = next(jreader)
            print("header: %s" % jheader)
        for row in jreader:
            samp, clstr = row[samp_cname_indx], row[samp_cluster_indx]
            if samp not in samp2clstr_dic:
                samp2clstr_dic[samp] = clstr
            else:
                raise Exception('Samp %s is duplicated' % samp)
            samps_good.append(samp)
            clstrs_lst.append(clstr)
            clstrs.add(clstr)
    print("Reading annot file... DONE")

    # # Check that every sample is in sampnames_all
    # print("Checking samps are found in bam...")
    # for samp in samps_good:
    #     assert samp in sampnames_all
    # print("Checking samps are found in bam... DONE")

    # check that no output bam files will clash with existing bams, exit if no overwrite ooption 
    if not args.overwrite:
        for clstr in clstrs:
            checktmppath = os.path.join(args.outdir, '.'.join([bname, clstr, "unsorted", "bam"]))
            checkoutpath = os.path.join(args.outdir, '.'.join([bname, clstr, "sorted", "bam"]))
            assert not os.path.exists(checktmppath)
            assert not os.path.exists(checkoutpath)

    dup_count = 0
    assign_count = 0
    unassign_count = 0

    print("Splitting bams into clusters ...")
    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        if args.add_chr_prefix:
            new_header = reformat_header(inbam.header, add_prefix = "chr")
        for clstr in clstrs:
            tmppath = os.path.join(args.outdir, '.'.join([bname, clstr, "unsorted", "bam"]))
            outpath = os.path.join(args.outdir, '.'.join([bname, clstr, "sorted", "bam"]))
            unsorteddic[clstr] = tmppath
            sorteddic[clstr] = outpath
            if not args.add_chr_prefix:
                writebamdic[clstr] = pysam.AlignmentFile(tmppath, "wb", template = inbam)
            else:
                writebamdic[clstr] = pysam.AlignmentFile(tmppath, "wb", header = new_header)
        for total_count, read in enumerate(inbam):
            if read.is_duplicate:
                dup_count += 1
                continue
            readsamp = readTag(read, samp_tag_id)
            if readsamp in samp2clstr_dic:
                # write to outpuet 
                clstr = samp2clstr_dic[readsamp]
                writebamdic[clstr].write(read)
                assign_count += 1
            else:
                unassign_count += 1
        # close bam files, sort, index and cleanup
        for clstr in clstrs:
            writebamdic[clstr].close()
            pysam.sort('-o', sorteddic[clstr], unsorteddic[clstr])
            pysam.index(sorteddic[clstr])
        # clean up
        for clstr in clstrs:
            os.remove(unsorteddic[clstr])
    print("Splitting bams into clusters ... DONE")

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code finished on: %c'))

if __name__ == '__main__':
    main()
