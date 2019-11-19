from singlecellmultiomics.molecule import Molecule
from singlecellmultiomics.fragment import Fragment
import singlecellmultiomics.universalBamTagger
import pysamiterators.iterators
import collections
import pysam

class MoleculeIterator():
    """Iterate over molecules in pysam.AlignmentFile or reads from a generator or list

    Example:
        >>> !wget https://github.com/BuysDB/SingleCellMultiOmics/blob/master/data/mini_nla_test.bam?raw=true -O mini_nla_test.bam
        >>> !wget https://github.com/BuysDB/SingleCellMultiOmics/blob/master/data/mini_nla_test.bam.bai?raw=true -O mini_nla_test.bam.bai
        >>> import pysam
        >>> from singlecellmultiomics.molecule import NlaIIIMolecule, MoleculeIterator
        >>> from singlecellmultiomics.fragment import NLAIIIFragment
        >>> import pysamiterators
        >>> alignments = pysam.AlignmentFile('mini_nla_test.bam')
        >>> for molecule in MoleculeIterator(
        >>>             alignments=alignments,
        >>>             moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
        >>>             fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
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

    Example:
        Using an iterator instead of a SAM/BAM file
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
            """


    """
    def __init__(self, alignments, moleculeClass=Molecule,
        fragmentClass=Fragment,
        check_eject_every=10_000, #bigger sizes are very speed benificial

        molecule_class_args={},# because the relative amount of molecules
                               # which can be ejected will be much higher
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

            pooling_method(int) : 0: no  pooling, 1: only compare molecules with the same sample id and hash

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
            self.molecules_per_cell =collections.defaultdict(list)#{hash:[], :}
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

        for reads in self.matePairIterator:

            if type(reads) is pysam.AlignedSegment:
                R1 = reads
                R2 = None
            elif len(reads)==2:
                R1,R2 = reads
            elif (type(reads) is list or type(reads) is tuple) and len(reads)==1:
                R1 = reads[0]
                R2 = None
            else:
                raise ValueError('Iterable not understood, supply either  pysam.AlignedSegment or lists of  pysam.AlignedSegment')
            # Make sure the sample/umi etc tags are placed:
            if self.perform_qflag:
                qf.digest([R1,R2])

            fragment = self.fragmentClass([R1,R2], **self.fragment_class_args)

            if not fragment.is_valid() :
                if self.yield_invalid:
                    m = self.moleculeClass(fragment, **self.molecule_class_args )
                    m.__finalise__()
                    yield m
                continue

            added = False
            if self.pooling_method==0:
                for molecule in self.molecules:
                    if molecule.add_fragment(fragment,use_hash=False):
                        added = True
                        break
            elif self.pooling_method==1:
                for molecule in self.molecules_per_cell[fragment.match_hash]:
                    if molecule.add_fragment(fragment,use_hash=True):
                        added = True
                        break

            if not added:
                if self.pooling_method==0:
                    self.molecules.append(self.moleculeClass(fragment, **self.molecule_class_args ))
                else:
                    self.molecules_per_cell[fragment.match_hash].append(
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
                        m = self.molecules.pop(i-j)
                        m.__finalise__()
                        yield m
                else:
                    for hash_group, molecules in self.molecules_per_cell.items():

                        to_pop=[]
                        for i,m in enumerate(molecules):
                            if m.can_be_yielded(current_chrom,current_position):
                                to_pop.append(i)
                                self.waiting_fragments-=len(m)
                                self.yielded_fragments+=len(m)

                        for i,j in enumerate(to_pop):
                            m = self.molecules_per_cell[hash_group].pop(i-j)
                            m.__finalise__()
                            yield m

        # Yield remains
        if self.pooling_method==0:
            for m in self.molecules:
                m.__finalise__()
            yield from iter(self.molecules)
        else:

            for hash_group,molecules in self.molecules_per_cell.items():
                for i,m in enumerate(molecules):
                    m.__finalise__()
                    yield m
        self._clear_cache()
