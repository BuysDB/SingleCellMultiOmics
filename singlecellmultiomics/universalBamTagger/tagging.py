#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from datetime import datetime
from os import remove
from pysam import AlignmentFile
from singlecellmultiomics.bamProcessing import sorted_bam_file
from uuid import uuid4
from copy import copy
from typing import Generator
import os
import gzip

def prefetch(contig, start, end, fetch_start,fetch_end,molecule_iterator_args):
    """ Prefetch selected region
    Prefetches
        AlleleResolver
    """

    new_kwarg_dict = {}
    for iterator_arg, iterator_value in molecule_iterator_args.items():
        if iterator_arg in ('molecule_class_args','fragment_class_args'):
            new_args = {}
            for key, value in molecule_iterator_args[iterator_arg].items():
                if key == 'features':
                    value = value.prefetch(contig,start,end)
                if key == 'mappability_reader':
                    value = value.prefetch(contig,start,end)
                new_args[key] = value
            new_kwarg_dict[iterator_arg] = new_args
        else:
            new_kwarg_dict[iterator_arg] = iterator_value
    return new_kwarg_dict


def write_job_gen_to_bed(job_gen: list, bed_path: str):
    with (gzip.open(bed_path,'wt') if bed_path.endswith('.gz') else open(bed_path,'w')) as o:
        for job_id, tasks in enumerate(job_gen):
            for i,(contig, start, end, fetch_start,fetch_end) in enumerate(tasks):
                o.write(f'{contig}\t{fetch_start}\t{fetch_end}\t{job_id}:{i}\t1\t+\t{start}\t{end}\n')

def run_tagging_task(alignments, output,
                    contig=None, start=None, end=None, fetch_start=None, fetch_end=None,
                    molecule_iterator_class=None,  molecule_iterator_args={},
                    read_groups=None, timeout_time=None, enable_prefetch=True, consensus_mode=None, no_source_reads=False):
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

    # There is three options: either supply only contig, or  start and end coordinates, or not supplying any coordinates:
    fetching = not all((x is None for x in (contig, start, end, fetch_start, fetch_end)))
    if fetching:
        if start is None and end is None:
            assert fetch_start is None and fetch_end is None, 'start and end need to be supplied'
            assert contig is not None

        elif start is not None and end is not None and fetch_start is None and fetch_end is None:
            fetch_start = start
            fetch_end = end
        elif all((x is None for x in (start, end, fetch_start, fetch_end))):
            assert contig is not None, 'A contig is required'

        else:
            if any((x is not None for x in (contig, start, end, fetch_start, fetch_end))):
                for variable in  ('contig', 'start', 'end', 'fetch_start', 'fetch_end'):
                    if not variable in locals() or locals()[variable] is None:
                        raise ValueError(f'{variable} has to be supplied')
            #assert not any((x is not None for x in (contig, start, end, fetch_start, fetch_end))), 'supply all these: contig, start, end, fetch_start, fetch_end'
        if enable_prefetch:
            molecule_iterator_args = prefetch(contig, start, end, fetch_start, fetch_end, molecule_iterator_args)
    time_start = datetime.now()


    def timeout_check_function(iteration, mol_iter, reads ):
        nonlocal time_start
        nonlocal timeout_time
        if timeout_time is None:
            return

        if (datetime.now()-time_start).total_seconds() > timeout_time:
            raise TimeoutError()


    total_molecules_written = 0
    for i, molecule in enumerate(
            molecule_iterator_class(alignments,  # Input alignments
                            contig=contig, start=fetch_start, end=fetch_end, # Region
                            # Set a callback function used to check if a problematic region is reached
                            progress_callback_function = timeout_check_function,
                            **molecule_iterator_args
                            )
        ):

        # Do not process molecules out of the defined region
        if fetch_start is not None:
            if fetching:
                cut_site_contig, cut_site_pos = None, None
                for fragment in molecule:
                    #@todo this needs better handling
                    r = fragment.get_site_location() # @todo: refactor to molecule function
                    if r is not None:
                        cut_site_contig, cut_site_pos = r
                        break
                if cut_site_contig is None:
                    # We cannot process this molecule using multiprocessing as there is a chance of a collision.
                    continue

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

        if consensus_mode is None:
            molecule.write_pysam(output)
        elif consensus_mode=='majority':
            molecule.write_pysam(output, consensus=True, no_source_reads=no_source_reads)
        else:
            raise ValueError(f'Unknown consensus method {consensus_mode}')

        total_molecules_written+=1

    return {'total_molecules_written': total_molecules_written,
            'time_start': time_start}



def run_tagging_tasks(args: tuple):
    """ Run tagging for one or more tasks

    Args:
        args (tuple): (alignments_path, temp_dir, timeout_time), arglist

    """

    (alignments_path, temp_dir, timeout_time), arglist = args

    target_file = f"{temp_dir}/{uuid4()}.bam"
    while os.path.exists(target_file):
        print(f'Collision at {target_file}')
        target_file = f"{temp_dir}/{uuid4()}.bam"


    timeout_tasks = []
    total_molecules = 0
    read_groups = dict()

    with AlignmentFile(alignments_path) as alignments:
        with sorted_bam_file(target_file, origin_bam=alignments, mode='wb', fast_compression=False,
                             read_groups=read_groups) as output:
            for task in arglist:
                try:
                    statistics = run_tagging_task(alignments, output, read_groups=read_groups, timeout_time=timeout_time, **task)
                    total_molecules += statistics.get('total_molecules_written', 0)
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
            print(f'Cleaning up failed for {target_file}')
            print(e)
            pass

    return None, meta


def generate_tasks(input_bam_path: str, temp_folder: str, job_gen: Generator, iteration_args: dict,
                   additional_args: dict,
                   max_time_per_segment: int = None) -> Generator:
    # Task generator:
    return (
        (
            (input_bam_path, temp_folder, max_time_per_segment),

            [{
                'contig': contig,
                'start': start,
                'end': end,
                'fetch_start': fetch_start,
                'fetch_end': fetch_end,
                **iteration_args, **additional_args

            } for contig, start, end, fetch_start, fetch_end in job]) for job in job_gen)
