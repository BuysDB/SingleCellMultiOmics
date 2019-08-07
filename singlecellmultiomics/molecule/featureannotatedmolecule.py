from singlecellmultiomics.molecule.molecule import Molecule
import collections
from more_itertools import consecutive_groups

# https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]

class FeatureAnnotatedMolecule(Molecule):
    """Molecule which is annotated with features (genes/exons/introns, .. )
    """

    def __init__(self, fragment, features, stranded=None, **kwargs):
        """
            Args:
                fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
                features (singlecellmultiomics.features.FeatureContainer) : container to use to obtain features from
                stranded : None; not stranded, False: same strand as R1, True: other strand
                **kwargs: extra args

        """
        Molecule.__init__(self,fragment,**kwargs)
        self.features = features
        self.hits = collections.defaultdict(set) #feature -> hit_bases
        self.stranded = stranded

    def annotate(self, method=1):
        """
            Args:
                method (int) : 0, obtain blocks and then obtain features. 1, try to obtain features for every aligned base

        """
        # When self.stranded is None, set to None strand. If self.stranded is True reverse the strand, otherwise copy the strand
        strand = None if self.stranded is None else '+-'[( not self.strand if self.stranded else self.strand )]

        if method==0:

            # Obtain all blocks:
            for start,end in find_ranges(
                sorted(list(set(
                    (ref_pos
                    for read in self.iter_reads()
                    for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False) ))))
                ):
                for hit in self.features.findFeaturesBetween(chromosome =self.chromosome, sampleStart=start, sampleEnd=end, strand=strand):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    self.hits[hit_ids].add((self.chromosome,(hit_start,hit_end)))


        else:

            for read in self.iter_reads():
                for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                    for hit in self.features.findFeaturesAt(chromosome=read.reference_name,lookupCoordinate=ref_pos,strand=strand):
                        hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                        self.hits[hit_ids].add((read.reference_name,ref_pos))
