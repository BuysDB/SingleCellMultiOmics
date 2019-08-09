from singlecellmultiomics.molecule import Molecule
from singlecellmultiomics.fragment import Fragment
import singlecellmultiomics.universalBamTagger
import pysamiterators.iterators
import collections
import pysam

class MoleculeIterator():

    def __init__(self, alignments, moleculeClass=Molecule,
        fragmentClass=Fragment,
        check_eject_every=10_000, #bigger sizes are very speed benificial
        # because the relative amount of molecules which can be ejected will be much higher
        molecule_class_args={},
        fragment_class_args={},
        perform_qflag=True,
        pooling_method=1,
        yield_invalid=False,
        queryNameFlagger=None,
        **pysamArgs):
        """Iterate over molecules in pysam.AlignmentFile

        Args:
            alignments (pysam.AlignmentFile) or iterable yielding tuples: Alignments to extract molecules from

            moleculeClass (pysam.FastaFile): Class to use for molecules.

            fragmentClass (pysam.FastaFile): Class to use for fragments.

            check_eject_every (int): Check for yielding every N reads.

            molecule_class_args (dict): arguments to pass to moleculeClass.

            fragment_class_args (dict): arguments to pass to fragmentClass.

            perform_qflag (bool):  Make sure the sample/umi etc tags are copied
                from the read name into bam tags

            pooling_method(int) : 0: no  pooling, 1: only compare molecules with the same sample id.

            yield_invalid (bool) : When true all fragments which are invalid will be yielded as a molecule

            queryNameFlagger(class) : class which contains the method digest(self, reads) which accepts pysam.AlignedSegments and adds at least the SM and RX tags

            **kwargs: arguments to pass to the pysam.AlignmentFile.fetch function

        Yields:
            molecule (Molecule): Molecule
        """
        if queryNameFlagger is None:
            queryNameFlagger = singlecellmultiomics.universalBamTagger.QueryNameFlagger()
        self.queryNameFlagger = queryNameFlagger

        self.alignments = alignments
        self.moleculeClass = moleculeClass
        self.fragmentClass = fragmentClass
        self.check_eject_every = check_eject_every
        self.molecule_class_args = molecule_class_args
        self.fragment_class_args = fragment_class_args
        self.perform_qflag = perform_qflag
        self.pysamArgs = pysamArgs
        self.matePairIterator=None
        self.pooling_method=pooling_method
        self.yield_invalid = yield_invalid
        self._clear_cache()

    def _clear_cache(self):
        """Clear cache containing non yielded molecules"""
        self.waiting_fragments = 0
        self.yielded_fragments = 0
        self.yielded_molecules = 0
        self.check_ejection_iter = 0
        if self.pooling_method==0:
            self.molecules=[]
        elif self.pooling_method==1:
            self.molecules_per_cell =collections.defaultdict(list)
        else:
            raise NotImplementedError()

    def __repr__(self):
        return f"""Molecule Iterator, generates fragments from {self.fragmentClass} into molecules based on {self.moleculeClass}.
        Yielded {self.yielded_fragments} fragments, {self.waiting_fragments} fragments are waiting to be ejected.
        {self.get_molecule_cache_size()} molecules cached.
        Mate pair iterator: {str(self.matePairIterator)}"""

    def get_molecule_cache_size(self):
        if self.pooling_method==0:
            return len(self.molecules)
        elif self.pooling_method==1:
            return sum( len(cell_molecules) for cell, cell_molecules in self.molecules_per_cell.items() )

        else:
            raise NotImplementedError()

    def __iter__(self):
        if self.perform_qflag:
            qf = self.queryNameFlagger

        self._clear_cache()

        self.waiting_fragments = 0

        # prepare the source iterator which generates the read pairs:
        if type(self.alignments)==pysam.libcalignmentfile.AlignmentFile:
            self.matePairIterator = pysamiterators.iterators.MatePairIterator(
                self.alignments,
                performProperPairCheck=False,
                **self.pysamArgs)
        else:
            # If an iterable is provided use this as read source:
            self.matePairIterator = self.alignments

        for R1,R2 in self.matePairIterator:
            # Make sure the sample/umi etc tags are placed:
            if self.perform_qflag:
                qf.digest([R1,R2])

            fragment = self.fragmentClass([R1,R2], **self.fragment_class_args)

            if not fragment.is_valid() :
                if self.yield_invalid:
                    yield self.moleculeClass(fragment, **self.molecule_class_args )
                continue

            added = False
            if self.pooling_method==0:
                for molecule in self.molecules:
                    if molecule.add_fragment(fragment,use_hash=False):
                        added = True
                        break
            elif self.pooling_method==1:
                for molecule in self.molecules_per_cell[fragment.sample]:
                    if molecule.add_fragment(fragment,use_hash=True):
                        added = True
                        break

            if not added:
                if self.pooling_method==0:
                    self.molecules.append(self.moleculeClass(fragment, **self.molecule_class_args ))
                else:
                    self.molecules_per_cell[fragment.sample].append(
                        self.moleculeClass(fragment, **self.molecule_class_args )
                    )

            self.waiting_fragments+=1
            self.check_ejection_iter += 1
            if self.check_ejection_iter>self.check_eject_every:
                current_chrom, _, current_position = fragment.get_span()
                if current_chrom is None:
                    continue

                self.check_ejection_iter=0
                if self.pooling_method==0:
                    to_pop = []
                    for i,m in enumerate(self.molecules):
                        if m.can_be_yielded(current_chrom,current_position):
                            to_pop.append(i)
                            self.waiting_fragments-=len(m)
                            self.yielded_fragments+=len(m)

                    for i,j in enumerate(to_pop):
                        yield self.molecules.pop(i-j)
                else:
                    for cell, cell_molecules in self.molecules_per_cell.items():
                        to_pop=[]
                        for i,m in enumerate(cell_molecules):
                            if m.can_be_yielded(current_chrom,current_position):
                                to_pop.append(i)
                                self.waiting_fragments-=len(m)
                                self.yielded_fragments+=len(m)

                        for i,j in enumerate(to_pop):
                            yield self.molecules_per_cell[cell].pop(i-j)


        # Yield remains
        if self.pooling_method==0:
            yield from iter(self.molecules)
        else:
            for cell, cell_molecules in self.molecules_per_cell.items():
                for i,m in enumerate(cell_molecules):
                    yield m
        self._clear_cache()
