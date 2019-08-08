import numpy as np

def coordinate_to_sliding_bin_locations(dp, bin_size, sliding_increment):
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
    start : int
        the start coordinate of the first overlapping bin
    end :int
        the end of the last overlapping bin

    start_id : int
        the index of the first overlapping bin
    end_id : int
        the index of the last overlapping bin

    """
    start_id = int( np.ceil(( (dp-bin_size)/sliding_increment ))   )
    start = start_id * sliding_increment
    end_id = int(np.floor(dp/sliding_increment))
    end = end_id * sliding_increment  + bin_size
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
    start,end,start_id,end_id = coordinate_to_sliding_bin_locations(point, bin_size,sliding_increment)
    return [ (i*sliding_increment,i*sliding_increment+bin_size) for i in range(start_id,end_id+1)]
