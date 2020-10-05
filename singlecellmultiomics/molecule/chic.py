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

        self.site_location = None
        Molecule.__init__(self, fragment, **kwargs)
        self.ligation_motif = None

        # Extracted from fragments:
        self.assignment_radius=None



    def _add_fragment(self,fragment):
        # Update the cut coordinate tho the (left most extreme value)
        self.assignment_radius =fragment.assignment_radius

        if fragment.site_location is not None:

            if self.site_location is None:
                self.site_location = [fragment.site_location[0],fragment.site_location[1]]

            elif fragment.strand: # Reverse:
                self.site_location[1] = max(fragment.site_location[1], self.site_location[1]) # this is the coordinate
            else:
                self.site_location[1] = min(fragment.site_location[1], self.site_location[1]) # this is the coordinate

        # else : writing a fragment which has no cut location associated

        Molecule._add_fragment(self, fragment)


    def get_cut_site(self):
        """For restriction based protocol data, obtain genomic location of cut site

        Returns:
            None if site is not available

            chromosome (str)
            position (int)
            strand (bool)
        """
        return (*self.site_location, self.strand)

    def write_tags(self):

        # Write DS tag when it could have been updated
        if self.assignment_radius is not None and self.assignment_radius>0 and len(self)>1:
            for frag in self:
                frag.set_meta('DS', self.site_location[1])

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
    def get_upstream_site(self, scan_extra_bp=3):
        read = None
        for fragment in self:
            if fragment[0] is not None:
                read = fragment[0]
                break

        if read is None:
            raise ValueError() # Read 1 is not mapped

        if self.strand==0:
            start = read.reference_start
            return self.reference.fetch( self.chromosome, start-4-scan_extra_bp, start)
        else:
            start = read.reference_end
            return self.reference.fetch( self.chromosome, start, start+4+scan_extra_bp)


    def __init__(self, fragment, reference, **kwargs):
        self.reference = reference
        CHICMolecule.__init__(self, fragment, reference=reference, **kwargs)

    def write_tags(self):
        try:
            site = self.get_upstream_site()
            self.set_meta('RZ',site)
            if 'N' in site:
                self.set_meta('dt','UNK')
            elif 'CATG' in site:
                self.set_meta('dt','NLA')
            else:
                self.set_meta('dt','CHIC')
        except ValueError:
            self.set_meta('dt','UNK')
            self.set_meta('RR','No_R1')




class AnnotatedCHICMolecule(CHICMolecule, FeatureAnnotatedMolecule):
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
        CHICMolecule.__init__(self, fragment, **kwargs)
        FeatureAnnotatedMolecule.__init__(self, fragment, features, **kwargs)
