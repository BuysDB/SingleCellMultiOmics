from singlecellmultiomics.molecule.molecule import Molecule
from singlecellmultiomics.molecule.featureannotatedmolecule import FeatureAnnotatedMolecule

class CHICMolecule(Molecule):
    """CHIC Molecule class

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule

        **kwargs: extra args

    """
    def __init__(self,fragment,
                **kwargs):
        Molecule.__init__(self,fragment,**kwargs)

    def write_tags(self):
        Molecule.write_tags(self)

    def is_valid(self,set_rejection_reasons=False, reference=None):
        # @todo
        return True

    def get_fragment_span_sequence(self,reference=None):
            """Obtain the sequence between the start and end of the molecule
            Args:
                reference(pysam.FastaFile) : reference  to use.
                    If not specified self.reference is used
            Returns:
                sequence (str)
            """
            if reference is None:
                if self.reference is None:
                    raise ValueError('Please supply a reference (PySAM.FastaFile)')
            reference = self.reference
            return reference.fetch(self.chromosome, self.spanStart, self.spanEnd).upper()
