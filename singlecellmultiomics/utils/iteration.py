from more_itertools import consecutive_groups

# https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0], group[0] # modification to always return tuples
        else:
            yield group[0], group[-1]
