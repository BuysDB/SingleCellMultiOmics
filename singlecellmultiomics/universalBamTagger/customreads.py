#!/usr/bin/env python3
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.utils import split_nth


class VaninsbergheQueryNameFlagger(DigestFlagger):
    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs)

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue
            origin, mi_tag, cell_barcode, umi, cell_index = read.query_name.rsplit(
                ':', 4)
            read.set_tag('MI', mi_tag)
            read.set_tag('RX', umi)
            read.set_tag('bi', int(cell_index))
            read.set_tag('SM', cell_barcode)

class BulkFlagger(DigestFlagger):
    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs)

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue

            read.set_tag('MI', "A")
            read.set_tag('RX', "A")
            read.set_tag('bi', 0)
            read.set_tag('SM', "BULK")



class CustomAssingmentQueryNameFlagger(DigestFlagger):
    """This query name flagger converts values between colons ":"  to tags"""

    def __init__(self, block_assignments, **kwargs):
        """Initialise CustomAssingmentQueryNameFlagger

        Args:
            block_assignments(list) : list of two letter codes to assign blocks to

        """
        self.block_assignments = block_assignments
        self.origin_colons = 7  # amount of ':' in original read name
        # Verify if all of the assignments are 2 letters:
        if not all((len(b) == 2 for b in block_assignments)):
            for b in block_assignments:
                if len(b) != 2:
                    raise ValueError(f'Tag {b} is not two letters long')

        DigestFlagger.__init__(self, **kwargs)

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue

            # Split original read name from added data
            origin, rest = split_nth(read.query_name, ':', self.origin_colons)
            # Reset the read name
            read.query_name = origin
            # Write the tags
            for tag, value in zip(self.block_assignments, rest.split(':')):
                read.set_tag(tag, value)
