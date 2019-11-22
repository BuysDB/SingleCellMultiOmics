#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysamiterators.iterators as pysamiterators
import itertools
import html


def cstr(s, color='black', weight=300):
    return f'<text style="color:{color}; font-weight:{weight}" >{s}</text>'


"""
Convert molecule to  HTML representation

"""


def visualise_molecule(molecule,
                       reference=None,  # pysam handle to reference
                       show_reference=False,  # Show reference sequence
                       margin_bases=0,
                       start_span=None,
                       end_span=None,
                       show_quals=True,
                       R1PrimerLength=4,
                       R2PrimerLength=6,
                       highlight=set()
                       ):

    html_str = f'<div style="font-family:monospace;line-height: 10px; font-size:12px; white-space: nowrap;">'

    reference_sequence = None

    molecule_flattened = list(itertools.chain.from_iterable(molecule))
    chrom, _start_span, _end_span = pysamiterators.getListSpanningCoordinates(
        molecule_flattened)

    # Calculate start and end of visualisation if not supplied by user:
    if end_span is None:
        end_span = _end_span
        end_span += margin_bases

    if start_span is None:
        start_span = _start_span
        start_span -= margin_bases

    start_span = max(0, start_span)
    # todo : end span
    end = end_span
    start = start_span

    # Obtain reference sequence
    if reference is not None and show_reference:
        reference_sequence = reference.fetch(chrom, start_span, end_span)
        html_str += f'ref : {reference_sequence}<br />'

    for R1, R2 in molecule:
        try:
            s, e = pysamiterators.getPairGenomicLocations(
                R1, R2,
                R1PrimerLength=4, R2PrimerLength=6)
            keep = range(s, e + 1)
        except Exception as e:
            keep = None
        for read in [R1, R2]:

            visBases = ['.'] * (end - start)
            visQuals = ['.'] * (end - start)
            emittedRef = ['.'] * (end - start)
            it = pysamiterators.ReadCycleIterator(
                read, with_seq=True, matches_only=True, reference=reference)

            for i, (cycle, qpos, refpos, refbase) in enumerate(it):
                if i == 0:
                    firstCycle = cycle

                visPos = refpos - start
                if not (visPos >= 0 and visPos < (end - start)):
                    continue

                refbase = refbase.upper()
                emittedRef[visPos] = refbase

                #print(qpos,refpos,refbase ,visPos,start, end)
                qbase = read.query_sequence[qpos]
                qqual = read.qual[qpos]

                if read.is_read2 and cycle < R2PrimerLength:
                    weight = 100
                    visBases[visPos] = cstr(qbase, color='blue', weight=weight)
                    visQuals[visPos] = cstr(
                        html.escape(
                            str(qqual)),
                        color='black',
                        weight=weight)

                elif keep is not None and refpos not in keep:

                    weight = 200
                    c = '#CCC'
                    visBases[visPos] = cstr(qbase, color=c, weight=weight)
                    visQuals[visPos] = cstr(
                        html.escape(
                            str(qqual)),
                        color=c,
                        weight=weight)

                elif refpos in highlight:

                    if qbase != refbase:
                        weight = 700
                        visBases[visPos] = cstr(
                            qbase, color='red', weight=weight)
                        weight = 400
                        visQuals[visPos] = cstr(
                            html.escape(
                                str(qqual)),
                            color='black',
                            weight=weight)
                    else:
                        weight = 700
                        visBases[visPos] = cstr(
                            qbase, color='green', weight=weight)
                        weight = 400
                        visQuals[visPos] = cstr(
                            html.escape(
                                str(qqual)),
                            color='black',
                            weight=weight)

                else:
                    if qbase != refbase:
                        weight = 700
                        # +'/'+cstr(refbase,color='red',weight=weight)
                        visBases[visPos] = cstr(
                            f'{qbase}', color='orange', weight=weight)
                        weight = 400
                        visQuals[visPos] = cstr(
                            html.escape(
                                str(qqual)),
                            color='black',
                            weight=weight)
                    else:
                        weight = 100
                        visBases[visPos] = cstr(
                            qbase, color='grey', weight=weight)
                        visQuals[visPos] = cstr(
                            html.escape(
                                str(qqual)),
                            color='#CCC',
                            weight=weight)

            html_str += f'{"RD1 " if read.is_read1 else "RD2 "}: ' + \
                ''.join(visBases) + f' cycle:{firstCycle}:{cycle}  {it.start} {it.len}'

            html_str += ' ' + cstr(
                f'SEQ:{read.get_tag("Is")}:{read.get_tag("Fc")}:{read.get_tag("La")}:{read.get_tag("Ti")}',
                'brown') + cstr(
                f' MD:{read.get_tag("MD")}',
                'black')

            html_str += '<br />'
            if show_quals:
                html_str += 'rdQ : ' + ''.join(visQuals) + '<br />'
            # if showEmittedRef:
        #        html_str+='eref: ' +''.join(emittedRef)+'<br />'

            #html_str+=f'<br />pred: {"".join(visPredict)}<br />conf: {"".join(visConf)}<br /><br /><br />'
    return html_str
