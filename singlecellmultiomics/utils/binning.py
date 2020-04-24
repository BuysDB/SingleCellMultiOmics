import numpy as np


def coordinate_to_sliding_bin_locations(dp, bin_size, sliding_increment):
    """
    Convert a single value to a list of overlapping bins

    Parameters
    -----   -----
    point : int
        coordinate to look up

    bin_size : int
        bin size

    sliding_increment : int
        sliding window offset, this is the increment between bins

    Returns
    -------
    start : int
        the start coordinate of the first overlapping bin
    end :int
        the end of the last overlapping bin

    start_id : int
        the index of the first overlapping bin
    end_id : int
        the index of the last overlapping bin

    """
    start_id = int(np.ceil(((dp - bin_size) / sliding_increment)))
    start = start_id * sliding_increment
    end_id = int(np.floor(dp / sliding_increment))
    end = end_id * sliding_increment + bin_size
    return start, end, start_id, end_id


def coordinate_to_bins(point, bin_size, sliding_increment):
    """
    Convert a single value to a list of overlapping bins

    Parameters
    ----------
    point : int
        coordinate to look up

    bin_size : int
        bin size

    sliding_increment : int
        sliding window offset, this is the increment between bins

    Returns
    -------
    list: [(bin_start,bin_end), .. ]

    """
    start, end, start_id, end_id = coordinate_to_sliding_bin_locations(
        point, bin_size, sliding_increment)
    return [(i * sliding_increment, i * sliding_increment + bin_size)
            for i in range(start_id, end_id + 1)]


def bp_chunked(job_generator, bp_per_job):
    """ Chunk an iterator containing coordinate sorted tasks in chunks of a total size of roughly bp_per_job

    Args:
        job_generator : iterable of commands, format (contig, start, end, *task)
        bp_per_job (int) : Amount of bp per chunk of jobs/tasks

    Yields:
        chunk(list) :  [(contig, start, end, *task),(contig, start, end, *task),..]

    @todo: contig is not used, this function expects that only bins on a single contig are supplied
    """
    bp_current = 0
    current_tasks = []
    for job in job_generator:
        start,end = job[1],job[2]
        bp_current += abs(end-start)
        current_tasks.append(job)

        if bp_current>=bp_per_job:
            yield current_tasks
            bp_current=0
            current_tasks=[]
    yield current_tasks
