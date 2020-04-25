#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Bio import bgzf
import os


class BlockZip():


    def verify(self):

        prev_contig = None
        prev_pos = None
        for line in self.bgzf_handle:
            if len(line) == 0:
                continue
            line_contig, line_pos, line_strand, rest = self.read_file_line(line)

            if prev_pos is not None and line_pos<prev_pos and line_contig==prev_contig:
                raise ValueError('Corrupted mapfile')

            prev_contig = line_contig
            prev_pos = line_pos

    def __init__(self, path, mode='r', read_all=False):
        """
        Store tabular information tied to genomic locations in a bgzipped file
        Args:
            path (str) : path to file
            mode (str) : mode, r: read, w: write

            read_all(bool) : when enabled all data is read from the file and the handles are closed
        """
        self.path = path
        self.index_path = f'{path}.idx'
        self.prev_contig = None
        self.mode = mode
        self.index = {}
        self.cache = {}

        if self.mode == 'w':
            self.bgzf_handle = bgzf.BgzfWriter(self.path, 'w')
            self.index_handle = open(self.index_path, 'wt')
        elif self.mode == 'r':
            if not os.path.exists(self.path):
                raise ValueError(f'BGZIP index file missing at {self.path}')
            self.bgzf_handle = bgzf.BgzfReader(self.path, 'rt')
            if not os.path.exists(self.index_path):
                raise ValueError(
                    f'BGZIP index file missing at {self.index_path}')
            self.index_handle = open(self.index_path, 'rt')

            for line in self.index_handle:
                contig, start = line.strip().split()
                self.index[contig] = int(start)

            if read_all:

                for line in self.bgzf_handle:
                    if len(line) == 0:
                        continue
                    line_contig, line_pos, line_strand, rest = self.read_file_line(line)
                    #print((line_pos, line_strand,rest))
                    if not line_contig in self.cache:
                        self.cache[line_contig] = {}
                    self.cache[line_contig][(line_pos, line_strand)] = rest
                    cpos = line_pos
                self.bgzf_handle.close()
                self.bgzf_handle=None
                self.index_handle.close()
                self.index_handle=None

        else:
            raise ValueError('Mode can be r or w')




    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.bgzf_handle.close()
        self.index_handle.close()

    def write(self, contig, position, strand, data):
        """
        Write information for location contig/postion/strand
        !! Write the contig data per contig, random mixing of contigs will
        result in a corrupted file
        Args:
            contig(str)
            postion(int)
            strand(bool)
            data(str)
        """
        assert(self.mode == 'w')
        if self.prev_contig is None or self.prev_contig != contig:
            self.index_handle.write(
                f'{contig}\t{int(self.bgzf_handle.tell())}\n')

        self.bgzf_handle.write(
            f'{contig}\t{position}\t{"+-"[strand]}\t{data}\n')
        self.prev_contig = contig

    def __iter__(self):
        "Get iterator going over all lines in the file"
        assert(self.mode == 'r')
        yield from iter(self.bgzf_handle)

    def read_file_line(self, line):
        line_contig, line_pos, line_strand, rest = line.strip().split(None, 3)
        return line_contig, int(line_pos), line_strand == '-', rest

    def __getitem__(self, contig_position_strand):
        """Obtain data at the supplied contig position and strand
        Args:
            contig_position_strand : tuple of (
                contig(str)
                postion(int)
                strand(bool))
        Returns:
            result (str) : data stored for the genomic location, returns None
            when no data is available
        """
        contig, position, strand = contig_position_strand
        if contig not in self.cache and contig in self.index:
            self.read_contig_to_cache(contig)
        if contig in self.cache:
            return self.cache[contig].get((position, strand), None)



    def read_contig_to_cache(self, contig, region_start=None, region_end=None):
        if contig not in self.index:
            return

        self.cache[contig] = {}
        # Seek to the start:
        self.bgzf_handle.seek(self.index[contig])

        while True:
            try:
                line = self.bgzf_handle.readline()
                if len(line) == 0:
                    break
                line_contig, line_pos, line_strand, rest = self.read_file_line(
                    line)
                if line_contig != contig:
                    break
                #print((line_pos, line_strand,rest))
                if region_start is not None and line_pos<region_start:
                    continue
                if region_end is not None and line_pos>region_end:
                    break
                self.cache[contig][(line_pos, line_strand)] = rest

            except Exception as e:
                raise

            if line_contig != contig:
                break
