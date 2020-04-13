#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class VariantWrapper:
    def __init__(self, variant, pos=None,contig=None,ref=None,alts=None,qual=0):
        if pos is not None:
            self.pos = pos
            self.contig = contig
            self.chrom= contig
            self.ref = ref
            self.alts = alts
            self.qual = qual

        else:
            self.pos = variant.pos
            self.contig = variant.contig
            self.chrom= variant.contig
            self.ref = variant.ref
            self.alts = variant.alts
            self.qual = variant.qual
    def __repr__(self):
        return f'VARIANT: {self.contig}:{self.pos} {self.ref} {self.alts} {self.qual}'
