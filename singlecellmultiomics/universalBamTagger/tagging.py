#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from datetime import datetime
from os import remove
from pysam import AlignmentFile
from singlecellmultiomics.bamProcessing import sorted_bam_file
from uuid import uuid4
from copy import copy


def prefetch(contig, start, end, fetch_start,fetch_end,molecule_iterator_args):
    """ Prefetch selected region
    Prefetches
        AlleleResolver
    """
    if 'molecule_class_args' in molecule_iterator_args:

        if 'allele_resolver' in molecule_iterator_args['molecule_class_args']:
            molecule_iterator_args['molecule_class_args']['allele_resolver'] = molecule_iterator_args['molecule_class_args']['allele_resolver']
        if 'features' in molecule_iterator_args['molecule_class_args']:
            molecule_iterator_args['molecule_class_args']['features'].prefetch(contig,start,end)



def run_tagging_task(alignments, output,
                    contig=None, start =None, end=None, fetch_start=None, fetch_end=None,
                    molecule_iterator_class=None,  molecule_iterator_args={},
                    read_groups=None, timeout_time=None ):
    """ Run tagging task for the supplied region

    Args:
        alignments (pysam.AlignmentFile) : handle to fetch reads from
        output (pysam.AlignmentFile or similar) : handle to write tagged reads to
        contig : (str)  : Contig to run the tagging for

        start (int) : only return molecules with a site identifier (DS) falling between start and end
        end (int) : see start
        fetch_start (int) : Start fetching reads from this coordinate onwards
        fetch_end (int) : Stop fetching reads at this coordinate (exclusive)

        molecule_iterator_class (class) : Class of the molecule iterator (not initialised, will be constructed using **molecule_iterator_args )
        molecule_iterator_args  (dict) : Arguments for the molecule iterator

    Returns:
        statistics : {'molecules_written':molecules_written}

    Raises:
        TimeoutError when the molecule_iterator_class decides tagging is taking too long

    """
    # Check some of the input arguments:
    assert alignments is not None
    assert output is not None
    assert molecule_iterator_class is not None
    # There is two options: either supply start and end coordinates, or not supplying any coordinates:
    fetching = not all((x is None for x in ( start, end, fetch_start, fetch_end)))
    if fetching:
        if start is not None and end is not None and fetch_start is None and fetch_end is None:
            fetch_start = start
            fetch_end = end
        else:
            assert not any((x is not None for x in (contig, start, end, fetch_start, fetch_end))), 'supply all these: contig, start, end, fetch_start, fetch_end'
        prefetch(contig, start, end, fetch_start,fetch_end,molecule_iterator_args)

    time_start = datetime.now()

    def timeout_check_function(iteration, mol_iter, reads ):
        nonlocal time_start
        nonlocal timeout_time
        if timeout_time is None:
            return
        if (datetime.now()-time_start).total_seconds() > timeout_time:
            raise TimeoutError()


    total_molecules_written = 0
    for i,molecule in enumerate(
            molecule_iterator_class(alignments,  # Input alignments
                            contig=contig, start=fetch_start, end=fetch_end, # Region
                            # Set a callback function used to check if a problematic region is reached
                            progress_callback_function = timeout_check_function,
                            **molecule_iterator_args
                            )
        ):

        if fetching:
            cut_site_contig, cut_site_pos = molecule[0].get_site_location() # @todo: refactor to molecule function

            # Stopping criteria:
            if cut_site_pos>=fetch_end:
                break

            # Skip molecules outside the desired region:
            if cut_site_contig!=contig or cut_site_pos<start or cut_site_pos>=end: # End is exclusive
                continue

        molecule.write_tags()

        # Update read groups @todo: reduce LOC here?
        if read_groups is not None:
            for fragment in molecule:
                rgid = fragment.get_read_group()
                if not rgid in read_groups:
                    read_groups[rgid] = fragment.get_read_group(True)[1]

        molecule.write_pysam(output)
        total_molecules_written+=1

    return {'total_molecules_written':total_molecules_written,
            'time_start':time_start}



def run_tagging_tasks(args):
    """ Run tagging for one or more tasks

    """

    (alignments_path, temp_dir, timeout_time), arglist = args

    target_file = f"{temp_dir}/{uuid4()}.bam"

    timeout_tasks = []
    total_molecules = 0
    read_groups = dict()
    with AlignmentFile(alignments_path) as alignments:
        with sorted_bam_file(target_file, origin_bam=alignments, mode='wbu',
            fast_compression=True, read_groups=read_groups) as output:
            for task  in arglist:
                try:
                    statistics = run_tagging_task(alignments, output, read_groups=read_groups, **task)
                    total_molecules += statistics.get('total_molecules_written',0)
                except TimeoutError:
                    timeout_tasks.append( task )


    meta = {
        'timeout_tasks' : timeout_tasks,
        'total_molecules' : total_molecules,
    }

    if total_molecules>0:
        return target_file, meta
    else:
        # Clean up ?
        try:
            remove(target_file)
            remove(f'{target_file}.bai')
        except Exception as e:
                pass

        return None, meta
