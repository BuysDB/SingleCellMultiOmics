#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Helpers to create small test bam files with pysam, shared between test modules."""
import os
import pysam


def write_misc_alt_bam(bam_path: str, include_unmapped: bool = True):
    """Write a sorted and indexed bam file containing one main chromosome, two alt
    contigs and (optionally) one unmapped read.

    Contents: main1 at chr1:1000, main2 at chr1:2000, alt1/alt2 on an alt contig,
    alt3 on the second alt contig, unmapped1 (when include_unmapped is True).
    All reads are tagged with SM and RX directly; the query name flagger leaves
    already tagged reads as-is.
    """
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'LN': 100000, 'SN': 'chr1'},                        # main chromosome
                     {'LN': 100000, 'SN': 'chrUn_KI270442v1'},            # alt contig
                     {'LN': 100000, 'SN': 'GL000256.3'}]}                 # alt contig
    sam_path = f'{bam_path}.tmp.sam'
    with pysam.AlignmentFile(sam_path, 'w', header=header) as out:
        def add(name, ref_id, pos, sm, rx, unmapped=False):
            read = pysam.AlignedSegment(out.header)
            read.query_name = name
            if unmapped:
                read.flag = 4
            else:
                read.reference_id = ref_id
                read.reference_start = pos
                read.cigarstring = '50M'
                read.query_sequence = 'A' * 50
                read.mapping_quality = 60
            read.is_read1 = True
            read.set_tag('SM', sm)
            read.set_tag('RX', rx)
            out.write(read)
        add('main1', 0, 1000, 'CELL_A', 'AAA')
        add('main2', 0, 2000, 'CELL_B', 'BBB')
        add('alt1', 1, 500, 'CELL_A', 'CCC')
        add('alt2', 1, 600, 'CELL_B', 'DDD')
        add('alt3', 2, 700, 'CELL_C', 'EEE')
        if include_unmapped:
            add('unmapped1', None, None, 'CELL_A', 'FFF', unmapped=True)
    pysam.sort('-o', bam_path, sam_path)
    pysam.index(bam_path)
    os.remove(sam_path)
    return bam_path


def count_reads_by_name(bam_path: str, names):
    """Return how often each of the given read names occurs in bam_path"""
    counts = {name: 0 for name in names}
    with pysam.AlignmentFile(bam_path) as f:
        for read in f:
            if read.query_name in counts:
                counts[read.query_name] += 1
    return counts
