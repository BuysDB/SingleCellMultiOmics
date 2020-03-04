#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pysam.libcalignmentfile import IteratorColumnRegion
from collections import Counter

def pileup_truncated(bam,contig, start, stop,**kwargs):
    """
    Obtain Pysam columns only at selected region

    contig(str) : contig to perform pileup
    start(int) : coordinate of first colummn (zero based)
    stop(int) : coordinate of last column (exclusive)
    **kwargs : arguments to pass to pysam.AlignmentFile.IteratorColumn
    """
    has_coord, rtid, rstart, rstop = bam.parse_region(contig, start, stop )
    yield from IteratorColumnRegion(bam,
                                    tid=rtid,
                                    start=rstart,
                                    stop=rstop,truncate=True,**kwargs)

def has_variant_reads(pysam_alignment_file, chrom, pos, alt, min_reads=2, stepper='nofilter'):
    """
    Check if the alignment file contains evidence for the supplied base

    Args:
        pysam_alignment_file(pysam.AlignmentFile) : file to check location

        chrom(str): name of contig

        pos(int) : position to check (zero based)

        alt(str): base to check
    """

    obs = Counter()
    for pile in pileup_truncated(pysam_alignment_file,chrom,pos,pos+1,stepper=stepper):
        if pos!=pile.reference_pos:
            continue
        for read in pile.pileups:
            if not read.is_del and not read.is_refskip:
                obs[read.alignment.query_sequence[read.query_position]]+=1
    return obs[alt]>=min_reads

def mate_pileup(alignments, contig, position,**kwargs):
    """
    Extract all fragments (R1, and R2) which overlap with the supplied position

    Example:
        >>> alignments = pysam.AlignmentFile('example.bam')
        >>> piled_reads = dict()
        >>> obs = collections.Counter()
        >>> pile_pos = 8774
        >>> pile_contig = '1'
        >>> add_pile_mates(alignments, piled_reads, pile_contig, pile_pos, obs )

        [[<pysam.libcalignedsegment.AlignedSegment at 0x7f8a7eca4ac8>,
          <pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecc83a8>],
         [<pysam.libcalignedsegment.AlignedSegment at 0x7f8a7d10ce28>,
          <pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecc8468>],
         [<pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecc8348>,
          <pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecae4c8>],
         [<pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecc8948>,
          <pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecae588>],
         [<pysam.libcalignedsegment.AlignedSegment at 0x7f8a7ecc4ee8>, ...

    Args:
        alignments(pysam.AlignmentFile) : Alignment file to extract reads from

        contig(str) : contig to perform pileup

        position(int) : coordinate to perform pileup

    Returns:
        list
    """

    piled_reads = dict()
    _mate_pileup(alignments=alignments,
                 piled_reads=piled_reads,
                 contig=contig,
                 position=position, obs=None, add_missing_only=False,**kwargs)
    return list(piled_reads.values())

def _mate_pileup(alignments, piled_reads, contig, position, obs=None, add_missing_only=False, max_depth=50000, stepper='nofilter',ignore_overlaps=False,ignore_orphans=False):
    """
    Extract all fragments (R1, and R2) which overlap with the supplied position

    Args:
        alignments(pysam.AlignmentFile) : Alignment file to extract reads from

        piled_reads(dict) : (Empty) Dictionary to which the reads are written

        contig(str) : contig to perform pileup

        position(int) : coordinate to perform pileup

        obs(collections.Counter) : Store base observation frequencies at the location of the pileup in this dictionary

        add_missing_only(bool): Only add missing mates to the existing dictionary (piled_reads)

    """
    for pile in alignments.pileup(contig,position,position+1,stepper=stepper,ignore_overlaps=ignore_overlaps,ignore_orphans=ignore_orphans,max_depth=max_depth):
        if position!=pile.reference_pos:
            continue

        for read in pile.pileups:

            if not add_missing_only and (read.is_del or read.is_refskip):
                continue

            if read.alignment.is_supplementary:
                continue

            if add_missing_only and not read.alignment.qname in piled_reads:
                continue

            if not read.alignment.qname in piled_reads:
                piled_reads[read.alignment.qname] = [None, None]

            piled_reads[read.alignment.qname][read.alignment.is_read2] = read.alignment

            if obs is not None:
                obs[read.alignment.query_sequence[read.query_position]]+=1

    if not add_missing_only:
        find_missing_mates(alignments, piled_reads)


def find_missing_mates(alignments, piled_reads):
    """
    Obtain missing mates

    Args:

        piled_reads(dict) : { query_name:[R1, R2], ... }

    """
    added=True
    piled_locations = set()
    while added:
        added=False
        for qname, reads in piled_reads.items():
            R1, R2 = reads

            if R1 is None and not R2.mate_is_unmapped:
                # Find R1 ..
                if look_for_mate(alignments, R2, piled_reads, piled_locations):
                    added=True
                    break

            if R2 is None and not R1.mate_is_unmapped:
                # Find R2 ..
                if look_for_mate(alignments, R1, piled_reads, piled_locations):
                    added=True
                    break


def look_for_mate(alignments, read, piled_reads, piled_locations):
    """

    Given a read, find the mate and possibly other missing reads in the piled_reads dictionary.
         Skips coordinates in piled_locations set.

     Args:
        alignments(pysam.AlignmentFile) : Alignment file to extract reads from

        read(pysam.AlignedSegment) : read to find mate for

        piled_reads(dict) : Query dictionary  { query_name:[R1, R2], ... }

        piled_locations(set) : Locations which have already been searched

    Returns:
        has_searched(bool) : Boolean which indicates if the location of the mate was already present in piled_locations and the location was not checked again.
    """
    location_descriptor = (read.next_reference_name, read.next_reference_start)
    if location_descriptor in piled_locations:
        return False

    piled_locations.add( location_descriptor )
    _mate_pileup(
                    alignments,
                    piled_reads,
                    contig = location_descriptor[0],
                    position = location_descriptor[1],
                    add_missing_only=True
                )

    if piled_reads[read.qname ][ read.is_read1] is None:
        raise(ValueError( f'Failed to find {read.qname} at {location_descriptor}' ))
    return True
