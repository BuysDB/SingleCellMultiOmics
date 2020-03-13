from singlecellmultiomics.molecule.molecule import Molecule
from singlecellmultiomics.molecule.featureannotatedmolecule import FeatureAnnotatedMolecule
import collections

class CHICMolecule(Molecule):
    """CHIC Molecule class

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
        **kwargs: extra args

    """

    def __init__(self, fragment,
                 **kwargs):
        Molecule.__init__(self, fragment, **kwargs)
        self.ligation_motif = None

    def write_tags(self):
        Molecule.write_tags(self)

    def is_valid(self, set_rejection_reasons=False):

        try:
            chrom, start, strand = self.get_cut_site()
        except Exception as e:
            if set_rejection_reasons:
                self.set_rejection_reason(f'no_cut_site_found')
            return False

        return True

    def get_fragment_span_sequence(self, reference=None):
        """Obtain the sequence between the start and end of the molecule

        Args:
            reference(pysam.FastaFile) : reference  to use.
                If not specified `self.reference` is used

        Returns:
            sequence (str)
        """
        if reference is None:
            if self.reference is None:
                raise ValueError('Please supply a reference (PySAM.FastaFile)')
        reference = self.reference
        return reference.fetch(
            self.chromosome,
            self.spanStart,
            self.spanEnd).upper()

    def __finalise__(self):
        # Obtain ligation motif(s)
        self.update_ligation_motif()
        Molecule.__finalise__(self)

    def update_ligation_motif(self):
        """
        Extract lh tag from associated reads and set most common one to the RZ tag
        """
        ligation_motif_counts = collections.Counter()
        for fragment in self:
            if fragment.ligation_motif is not None:
                ligation_motif_counts[fragment.ligation_motif]+=1
        if len(ligation_motif_counts)>0:
            self.ligation_motif = ligation_motif_counts.most_common(1)[0][0]
            self.set_meta('RZ', self.ligation_motif)


class CHICNLAMolecule(Molecule):
    """CHIC NLA Molecule class

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
        **kwargs: extra args

    """
    def get_upstream_site(self):
        read = None
        for fragment in self:
            if fragment[0] is not None:
                read = fragment[0]
                break

        if read is None:
            raise ValueError() # Read 1 is not mapped

        if self.strand==1:
            start = read.reference_start
            return self.reference.fetch( self.chromosome, start-4, start)
        else:
            start = read.reference_end
            return self.reference.fetch( self.chromosome, start, start+4)


    def __init__(self, fragment, reference, **kwargs):
        self.reference = reference
        CHICMolecule.__init__(self, fragment, reference=reference, **kwargs)

    def write_tags(self):
        try:
            site = self.get_upstream_site()
            if 'N' in site:
                self.set_meta('dt','UNK')
            elif site=='CATG':
                self.set_meta('dt','NLA')
            else:
                self.set_meta('dt','CHIC')
        except ValueError:
            self.set_meta('dt','UNK')
            self.set_meta('RR','No_R1')




class AnnotatedCHICMolecule(FeatureAnnotatedMolecule, CHICMolecule):
    """Chic Molecule which is annotated with features (genes/exons/introns, .. )

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
        features (singlecellmultiomics.features.FeatureContainer) : container to use to obtain features from
        **kwargs: extra args

    """

    def write_tags(self):
        CHICMolecule.write_tags(self)
        FeatureAnnotatedMolecule.write_tags(self)

    def __init__(self, fragment, features, **kwargs):
        FeatureAnnotatedMolecule.__init__(self, fragment, features, **kwargs)
        CHICMolecule.__init__(self, fragment, **kwargs)
