from singlecellmultiomics.molecule.molecule import Molecule
from singlecellmultiomics.molecule.featureannotatedmolecule import FeatureAnnotatedMolecule


class VASA(FeatureAnnotatedMolecule):
    """VASA seq molecule

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
        features (singlecellmultiomics.features.FeatureContainer) : container to use to obtain features from
        **kwargs: extra args

    """

    def __init__(self, fragment, features, **kwargs):
        FeatureAnnotatedMolecule.__init__(self, fragment, features, **kwargs)
