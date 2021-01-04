#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import argparse

def convert_scmo_to_cellranger(scmo_bam, cellranger_bam, n_threads=4):
    """
    Convert a SCMO bam file to a cellranger bam file
    Adds these tags:
        BC to CB
        RX to UB
    """

    with pysam.AlignmentFile(scmo_bam,threads=n_threads) as al:
        with pysam.AlignmentFile(cellranger_bam,  'wb', header=al.header, threads=n_threads) as ao:
            for read in al:
                if read.has_tag('BC'):
                    read.set_tag('CB',read.get_tag('BC'))
                if read.has_tag('RX'):
                    read.set_tag('UB',read.get_tag('RX'))

                ao.write(read)

    pysam.index(cellranger_bam)


def convert_scmo_to_dropseq(scmo_bam, dropseq_bam, n_threads=4):
    """
    Convert a SCMO bam file to a dropseq bam file
    Adds these tags:
        BC to XC
        RX to XM
    """
    with pysam.AlignmentFile(scmo_bam,threads=n_threads) as al:
        with pysam.AlignmentFile(dropseq_bam,  'wb', header=al.header, threads=n_threads) as ao:
            for read in al:
                if read.has_tag('BC'):
                    read.set_tag('XC',read.get_tag('BC'))
                if read.has_tag('RX'):
                    read.set_tag('XM',read.get_tag('RX'))
                ao.write(read)

    pysam.index(dropseq_bam)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Convert SCMO bam file to cellranger or dropseq bam file
    """)
    argparser.add_argument('inputbamfile', type=str)
    argparser.add_argument('convertedbamfile', type=str)

    og = argparser.add_argument_group("Output")
    #og.add_argument('-bed', type=str, help='Bed file to write methylation calls to')
    og.add_argument('-fmt', type=str, help='Format to convert to cellranger/dropseq', required=True)

    args = argparser.parse_args()

    if args.fmt == 'cellranger':
        convert_scmo_to_cellranger(args.inputbamfile, args.convertedbamfile)
    elif args.fmt == 'dropseq':
        convert_scmo_to_dropseq(args.inputbamfile, args.convertedbamfile)
