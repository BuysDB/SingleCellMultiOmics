# Fastq iterator class, Buys de Barbanson
import collections
import gzip
import os

FastqRecord = collections.namedtuple(
    'FastqRecord', 'header sequence plus qual')


class FastqIterator():
    """FastqIterator, iterates over one or more fastq files."""

    def __init__(self, *args):
        """Initialise  FastqIterator.

        Argument(s):
        path to fastq file, path to fastq file 2 , ...
        example: for rec1, rec2 in FastqIterator('./R1.fastq', './R2.fastq'):
        """

        self.handles = tuple(
            gzip.open(path, 'rt')
            # Load as GZIP when the extension is .gz
            if os.path.splitext(path)[1] == '.gz' else open(path, 'r')
            for path in args
        )
        self.readIndex = 0

    def _readFastqRecord(self, handle):
        # Read four lines and load them into a FastqRecord
        return(
            #FastqRecord(*tuple(handle.readline().rstrip() for i in range(4)))
            FastqRecord(
                handle.readline().rstrip(),
                handle.readline().rstrip(),
                handle.readline().rstrip(),
                handle.readline().rstrip()
            )
        )

    def __iter__(self):
        """Exectuted upon generator initiation."""
        return(self)

    def __next__(self):
        """Obtain the next fastq record for all opened files."""
        self.readIndex += 1  # Increment the current read counter
        records = tuple(self._readFastqRecord(handle)
                        for handle in self.handles)
        # Stop when empty records are being returned; the file end is reached
        if any((len(rec.header) == 0 for rec in records)):
            raise StopIteration
        return(records)
