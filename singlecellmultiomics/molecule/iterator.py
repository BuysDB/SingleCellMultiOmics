#!/usr/bin/env python
# -*- coding: utf-8 -*-
from singlecellmultiomics.molecule import Molecule
from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.utils.prefetch import initialise_dict, initialise
from singlecellmultiomics.universalBamTagger import QueryNameFlagger
import pysamiterators.iterators
import collections
import pysam


class ReadIterator(pysamiterators.iterators.MatePairIterator):
    def __next__(self):

        try:
            rec = next(self.iterator)
            return tuple((rec, None))
        except StopIteration:
            raise

class MoleculeIterator():
    """Iterate over molecules in pysam.AlignmentFile or reads from a generator or list

    Example:
        >>> !wget https://github.com/BuysDB/SingleCellMultiOmics/blob/master/data/mini_nla_test.bam?raw=true -O mini_nla_test.bam
        >>> !wget https://github.com/BuysDB/SingleCellMultiOmics/blob/master/data/mini_nla_test.bam.bai?raw=true -O mini_nla_test.bam.bai
        >>> import pysam
        >>> from singlecellmultiomics.molecule import NlaIIIMolecule, MoleculeIterator
        >>> from singlecellmultiomics.fragment import NlaIIIFragment
        >>> import pysamiterators
        >>> alignments = pysam.AlignmentFile('mini_nla_test.bam')
        >>> for molecule in MoleculeIterator(
        >>>             alignments=alignments,
        >>>             molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
        >>>             fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
        >>>         ):
        >>>     break
        >>> molecule
        NlaIIIMolecule
        with 1 assinged fragments
        Allele :No allele assigned
            Fragment:
            sample:APKS1P25-NLAP2L2_57
            umi:CCG
            span:chr1 164834728-164834868
            strand:+
            has R1: yes
            has R2: no
            randomer trimmed: no
            DS:164834865
        RS:0
        RZ:CAT
        Restriction site:('chr1', 164834865)


    It is also possible to supply and iterator instead of a SAM/BAM file handle
    Example:

        >>> from singlecellmultiomics.molecule import MoleculeIterator
        >>> from singlecellmultiomics.fragment import Fragment
        >>> import pysam
        >>> # Create SAM file to write some example reads to:
        >>> test_sam = pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000])
        >>> read_A = pysam.AlignedSegment(test_sam.header)
        >>> read_A.set_tag('SM','CELL_1')
        >>> read_A.set_tag('RX','CAT')
        >>> read_A.reference_name = 'chr1'
        >>> read_A.reference_start = 100
        >>> read_A.query_sequence = 'ATCGGG'
        >>> read_A.cigarstring = '6M'
        >>> read_A.mapping_quality = 60
        >>> # Create a second read which is a duplicate of the previous
        >>> read_B = pysam.AlignedSegment(test_sam.header)
        >>> read_B.set_tag('SM','CELL_1')
        >>> read_B.set_tag('RX','CAT')
        >>> read_B.reference_name = 'chr1'
        >>> read_B.reference_start = 100
        >>> read_B.query_sequence = 'ATCGG'
        >>> read_B.cigarstring = '5M'
        >>> read_B.mapping_quality = 60
        >>> # Create a thids read which is belonging to another cell
        >>> read_C = pysam.AlignedSegment(test_sam.header)
        >>> read_C.set_tag('SM','CELL_2')
        >>> read_C.set_tag('RX','CAT')
        >>> read_C.reference_name = 'chr1'
        >>> read_C.reference_start = 100
        >>> read_C.query_sequence = 'ATCGG'
        >>> read_C.cigarstring = '5M'
        >>> read_C.mapping_quality = 60
        >>> # Set up an iterable containing the reads:
        >>> reads = [  read_A,read_B,read_C ]
        >>> molecules = []
        >>> for molecule in MoleculeIterator( reads ):
        >>>     print(molecule)
        Molecule
            with 2 assinged fragments
            Allele :No allele assigned
                Fragment:
                sample:CELL_1
                umi:CAT
                span:chr1 100-106
                strand:+
                has R1: yes
                has R2: no
                randomer trimmed: no
            Fragment:
                sample:CELL_1
                umi:CAT
                span:chr1 100-105
                strand:+
                has R1: yes
                has R2: no
                randomer trimmed: no
        Molecule
                with 1 assinged fragments
                Allele :No allele assigned
                    Fragment:
                    sample:CELL_2
                    umi:CAT
                    span:chr1 100-105
                    strand:+
                    has R1: yes
                    has R2: no
                    randomer trimmed: no


    In the next example the molecules overlapping with a single location on chromosome `'1'` position `420000` are extracted
    Don't forget to supply `check_eject_every = None`, this allows non-sorted data to be passed to the MoleculeIterator.

    Example:

        >>> from singlecellmultiomics.bamProcessing import mate_pileup
        >>> from singlecellmultiomics.molecule import MoleculeIterator
        >>> with pysam.AlignmentFile('example.bam') as alignments:
        >>>     for molecule in MoleculeIterator(
        >>>         mate_pileup(alignments, contig='1', position=420000, check_eject_every=None)
        >>>     ):
        >>>         pass


    Warning:
        Make sure the reads being supplied to the MoleculeIterator sorted by genomic coordinate! If the reads are not sorted set `check_eject_every=None`
    """

    def __init__(self, alignments, molecule_class=Molecule,
                 fragment_class=Fragment,
                 check_eject_every=10_000,  # bigger sizes are very speed benificial
                 molecule_class_args={},  # because the relative amount of molecules
                 # which can be ejected will be much higher
                 fragment_class_args={},
                 perform_qflag=True,
                 pooling_method=1,
                 yield_invalid=False,
                 yield_overflow=True,
                 query_name_flagger=None,
                 ignore_collisions=True, # Ignore read-dupe collisions
                 every_fragment_as_molecule=False,
                 yield_secondary =  False,
                 yield_supplementary= False,
                 max_buffer_size=None,  #Limit the amount of stored reads, when this value is exceeded, a MemoryError is thrown
                 iterator_class = pysamiterators.iterators.MatePairIterator,
                 skip_contigs=None,
                 progress_callback_function=None,
                 min_mapping_qual = None,
                 perform_allele_clustering = False,

                 **pysamArgs):
        """Iterate over molecules in pysam.AlignmentFile

        Args:
            alignments (pysam.AlignmentFile) or iterable yielding tuples: Alignments to extract molecules from

            molecule_class (pysam.FastaFile): Class to use for molecules.

            fragment_class (pysam.FastaFile): Class to use for fragments.

            check_eject_every (int): Check for yielding every N reads. When None is supplied, all reads are kept into memory making coordinate sorted data not required.

            molecule_class_args (dict): arguments to pass to molecule_class.

            fragment_class_args (dict): arguments to pass to fragment_class.

            perform_qflag (bool):  Make sure the sample/umi etc tags are copied
                from the read name into bam tags

            pooling_method(int) : 0: no  pooling, 1: only compare molecules with the same sample id and hash

            ignore_collisions   (bool) : parameter passed to pysamIterators MatePairIterator, this setting will throw a fatal error when a duplicated read is found

            yield_invalid (bool) : When true all fragments which are invalid will be yielded as a molecule

            yield_overflow(bool) : When true overflow fragments are yielded as separate molecules

            query_name_flagger(class) : class which contains the method digest(self, reads) which accepts pysam.AlignedSegments and adds at least the SM and RX tags

            every_fragment_as_molecule(bool): When set to true all valid fragments are emitted as molecule with one associated fragment, this is a way to disable deduplication.

            yield_secondary (bool):  When true all secondary alignments will be yielded as a molecule

            iterator_class : Class name of class which generates mate-pairs out of a pysam.AlignmentFile either (pysamIterators.MatePairIterator or pysamIterators.MatePairIteratorIncludingNonProper)

            skip_contigs (set) : Contigs to skip

            min_mapping_qual(int) : Dont process reads with a mapping quality lower than this value. These reads are not yielded as molecules!

            **kwargs: arguments to pass to the pysam.AlignmentFile.fetch function

        Yields:
            molecule (Molecule): Molecule
        """
        if query_name_flagger is None:
            query_name_flagger = QueryNameFlagger()
        self.query_name_flagger = query_name_flagger
        self.skip_contigs = skip_contigs if skip_contigs is not None else set()
        self.alignments = alignments
        self.molecule_class = molecule_class
        self.fragment_class = fragment_class
        self.check_eject_every = check_eject_every
        self.molecule_class_args = initialise_dict(molecule_class_args)
        self.fragment_class_args = initialise_dict(fragment_class_args)
        self.perform_qflag = perform_qflag
        self.pysamArgs = pysamArgs
        self.matePairIterator = None
        self.ignore_collisions = ignore_collisions
        self.pooling_method = pooling_method
        self.yield_invalid = yield_invalid
        self.yield_overflow = yield_overflow
        self.every_fragment_as_molecule = every_fragment_as_molecule
        self.progress_callback_function = progress_callback_function
        self.iterator_class = iterator_class
        self.max_buffer_size=max_buffer_size
        self.min_mapping_qual = min_mapping_qual
        self.perform_allele_clustering = perform_allele_clustering

        self._clear_cache()

    def _clear_cache(self):
        """Clear cache containing non yielded molecules"""
        self.waiting_fragments = 0
        self.yielded_fragments = 0
        self.deleted_fragments = 0
        self.check_ejection_iter = 0
        if self.pooling_method == 0:
            self.molecules = []
        elif self.pooling_method == 1:
            self.molecules_per_cell = collections.defaultdict(
                list)  # {hash:[], :}
        else:
            raise NotImplementedError()

    def __repr__(self):
        return f"""Molecule Iterator, generates fragments from {self.fragment_class} into molecules based on {self.molecule_class}.
        Yielded {self.yielded_fragments} fragments, {self.waiting_fragments} fragments are waiting to be ejected. {self.deleted_fragments} fragments rejected.
        {self.get_molecule_cache_size()} molecules cached.
        Mate pair iterator: {str(self.matePairIterator)}"""

    def get_molecule_cache_size(self):
        if self.pooling_method == 0:
            return len(self.molecules)
        elif self.pooling_method == 1:
            return sum(len(cell_molecules) for cell,
                       cell_molecules in self.molecules_per_cell.items())

        else:
            raise NotImplementedError()


    def yield_func(self, molecule_to_be_emitted):
        if self.perform_allele_clustering:
            if molecule_to_be_emitted.can_be_split_into_allele_molecules:
                new_molecules = molecule_to_be_emitted.split_into_allele_molecules()
                if len(new_molecules)>1:
                    yield from new_molecules
                else:
                    yield molecule_to_be_emitted
            else:
                yield molecule_to_be_emitted
        else:
            yield molecule_to_be_emitted


    def __iter__(self):
        if self.perform_qflag:
            qf = self.query_name_flagger

        self._clear_cache()
        self.waiting_fragments = 0
        # prepare the source iterator which generates the read pairs:
        if isinstance(self.alignments, pysam.libcalignmentfile.AlignmentFile):

            # Don't pass the ignore_collisions to other classes than the matepair iterator
            if self.iterator_class == pysamiterators.iterators.MatePairIterator:
                self.matePairIterator = self.iterator_class(
                    self.alignments,
                    performProperPairCheck=False,
                    ignore_collisions=self.ignore_collisions,
                    **self.pysamArgs)
            else:
                self.matePairIterator = self.iterator_class(
                    self.alignments,
                    performProperPairCheck=False,
                    **self.pysamArgs)

        else:
            # If an iterable is provided use this as read source:
            self.matePairIterator = self.alignments

        for iteration,reads in enumerate(self.matePairIterator):

            if self.progress_callback_function is not None and iteration%500==0:
                self.progress_callback_function(iteration, self, reads)

            if isinstance(reads, pysam.AlignedSegment):
                R1 = reads
                R2 = None
            elif len(reads) == 2:
                R1, R2 = reads
            elif (isinstance(reads, list) or isinstance(reads, tuple)) and len(reads) == 1:
                R1 = reads[0]
                R2 = None
            else:
                raise ValueError(
                    'Iterable not understood, supply either  pysam.AlignedSegment or lists of  pysam.AlignedSegment')

            # skip_contigs
            if len(self.skip_contigs)>0:
                keep = False
                for read in reads:
                    if read is not None and read.reference_name not in self.skip_contigs:
                        keep = True
                if not keep:
                    continue

            if self.min_mapping_qual is not None:
                keep = True
                for read in reads:
                    if read is not None and read.mapping_quality<self.min_mapping_qual:
                        self.deleted_fragments+=1
                        keep=False
                if not keep:
                    continue

            # Make sure the sample/umi etc tags are placed:
            if self.perform_qflag:
                qf.digest([R1, R2])

            fragment = self.fragment_class([R1, R2], **self.fragment_class_args)

            if not fragment.is_valid():
                if self.yield_invalid:
                    m = self.molecule_class(
                        fragment, **self.molecule_class_args)
                    m.__finalise__()
                    yield m
                else:
                    self.deleted_fragments+=1
                continue

            if self.every_fragment_as_molecule:
                m = self.molecule_class(fragment, **self.molecule_class_args)
                m.__finalise__()
                yield m
                continue

            added = False
            try:
                if self.pooling_method == 0:
                    for molecule in self.molecules:
                        if molecule.add_fragment(fragment, use_hash=False):
                            added = True
                            break
                elif self.pooling_method == 1:
                    for molecule in self.molecules_per_cell[fragment.match_hash]:
                        if molecule.add_fragment(fragment, use_hash=True):
                            added = True
                            break
            except OverflowError:
                # This means the fragment does belong to a molecule, but the molecule does not accept any more fragments.
                if self.yield_overflow:
                    m = self.molecule_class(fragment, **self.molecule_class_args)
                    m.set_rejection_reason('overflow')
                    m.__finalise__()
                    yield from self.yield_func(m)
                else:
                    self.deleted_fragments+=1
                continue

            if not added:
                if self.pooling_method == 0:
                    self.molecules.append(self.molecule_class(
                        fragment, **self.molecule_class_args))
                else:
                    self.molecules_per_cell[fragment.match_hash].append(
                        self.molecule_class(fragment, **self.molecule_class_args)
                    )

            self.waiting_fragments += 1
            self.check_ejection_iter += 1

            if self.max_buffer_size is not None and self.waiting_fragments>self.max_buffer_size:
                raise MemoryError(f'max_buffer_size exceeded with {self.waiting_fragments} waiting fragments')

            if self.check_eject_every is not None and self.check_ejection_iter > self.check_eject_every:
                current_chrom, _, current_position = fragment.get_span()
                if current_chrom is None:
                    continue

                self.check_ejection_iter = 0
                if self.pooling_method == 0:
                    to_pop = []
                    for i, m in enumerate(self.molecules):
                        if m.can_be_yielded(current_chrom, current_position):
                            to_pop.append(i)
                            self.waiting_fragments -= len(m)
                            self.yielded_fragments += len(m)

                    for i, j in enumerate(to_pop):
                        m = self.molecules.pop(i - j)
                        m.__finalise__()
                        yield from self.yield_func(m)
                else:
                    for hash_group, molecules in self.molecules_per_cell.items():

                        to_pop = []
                        for i, m in enumerate(molecules):
                            if m.can_be_yielded(
                                    current_chrom, current_position):
                                to_pop.append(i)
                                self.waiting_fragments -= len(m)
                                self.yielded_fragments += len(m)

                        for i, j in enumerate(to_pop):
                            m = self.molecules_per_cell[hash_group].pop(i - j)
                            m.__finalise__()
                            yield from self.yield_func(m)

        # Yield remains
        if self.pooling_method == 0:
            for m in self.molecules:
                m.__finalise__()
                yield from self.yield_func(m)
            #yield from iter(self.molecules)
        else:

            for hash_group, molecules in self.molecules_per_cell.items():
                for i, m in enumerate(molecules):
                    m.__finalise__()
                    yield from self.yield_func(m)
        self._clear_cache()
