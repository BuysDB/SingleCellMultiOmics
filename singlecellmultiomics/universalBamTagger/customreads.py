#!/usr/bin/env python3
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger


class VaninsbergheQueryNameFlagger(DigestFlagger):
    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def digest(self, reads):
        for read in reads:
            origin, mi_tag, cell_barcode, umi, cell_index = read.query_name.rsplit(':',4)
            read.set_tag('MI', mi_tag)
            read.set_tag('RX', mi_tag)
            read.set_tag('BI', int(cell_index))
            read.set_tag('SM', cell_barcode)
