#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysamiterators.iterators as pysamiterators
import itertools



"""
Convert molecule to  HTML representation

"""
visualise_molecule( molecule, reference_handle=None, show_reference=False, margin_bases=None, start_span=None, end_span = None  ):

    html_str = f'<div style="font-family:monospace;line-height: 10px; font-size:12px; white-space: nowrap;">'

    reference_sequence = None

    molecule_flattened = itertools.chain.from_iterable( molecule )
    chrom, _start_span, _end_span = pysamiterators.getListSpanningCoordinates(moleculeFlat)

    # Calculate start and end of visualisation if not supplied by user:
    if end_span is None:
        end_span = _end_span
        end_span+=margin_bases

    if start_span is None:
        start_span = _start_span
        start_span-=margin_bases

    start_span = max(0,start_span)
    #todo : end span

    # Obtain reference sequence
    if reference_handle is not None and show_reference:
        reference_sequence = reference.fetch(chrom,start_span, end_span)
