from singlecellmultiomics.molecule.molecule import Molecule
import collections

class FeatureAnnotatedMolecule(Molecule):
    """Molecule which is annotated with features (genes/exons/introns, .. )
    """

    def __init__(self, fragment, features, **kwargs):
        """
            Args:
                fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
                features (singlecellmultiomics.features.FeatureContainer) : container to use to obtain features from
                **kwargs: extra args

        """
        Molecule.__init__(self,fragment,**kwargs)
        self.features = features
        self.hits = collections.defaultdict(set) #feature -> hit_bases

    def annotate(self):
        for read in self.iter_reads():
            for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                for hit in self.features.findFeaturesAt(chromosome=read.reference_name,lookupCoordinate=ref_pos,strand=None):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    self.hits[hit_ids].add((read.reference_name,ref_pos))
