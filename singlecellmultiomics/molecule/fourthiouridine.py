from singlecellmultiomics.molecule import Molecule


class FourThiouridine(Molecule):
    def __init__(self, fragments=None, classifier=None, **kwargs):
        """ FourThiouridine Molecule class

        Args:
            fragments(list) :  list of fragments to associate with the Molecule
            classifier : fitted sklearn classifier, when supplied this classifier is used to obtain a consensus from which the methylation calls are generated.

        """
        Molecule.__init__(self, fragments=fragments, **kwargs)
        self.classifier = classifier
        self.gene = None

    def __finalise__(self):
        super().__finalise__()
        self.obtain_conversions(self.classifier)

        for frag in self:
            if frag.gene is not None:
                self.gene = frag.gene

    def is_valid(self, set_rejection_reasons=False):
        if not super().is_valid(set_rejection_reasons=set_rejection_reasons):
            return False

        try:
            consensus = self.get_consensus()
        except ValueError:
            if set_rejection_reasons:
                self.set_rejection_reason('no_consensus')
            return False
        except TypeError:
            if set_rejection_reasons:
                self.set_rejection_reason('getPairGenomicLocations_failed')
            return False

        return True

    def obtain_conversions(self, classifier=None):
        """ This methods obtains the amount of converted bases and stores them to self.converted_bases and the 4U tag
            Args:
                classifier : classifier used for consensus determination
            returns:
                None
        """

        # Find all aligned positions and corresponding reference bases:
        aligned_reference_positions = {}  # (chrom,pos)->base
        for read in self.iter_reads():
            for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                    with_seq=True, matches_only=True):
                aligned_reference_positions[(
                    read.reference_name, ref_pos)] = ref_base.upper()

        # Obtain consensus:
        try:
            consensus = self.get_consensus(classifier=classifier)
        except ValueError:
            raise ValueError(
                'Cannot obtain a safe consensus for this molecule')

        # look for T > C conversions
        self.converted_bases = 0
        conversions = {}

        for location, reference_base in aligned_reference_positions.items():
            if location not in consensus:
                continue
            if (not self.strand and reference_base == 'T' and consensus[location] == 'C') or \
                    self.strand and reference_base == 'A' and consensus[location] in 'G':
                conversions[location] = {
                    'ref': reference_base, 'obs': consensus[location]}
                self.converted_bases += 1
        self.set_meta('4U', self.converted_bases)
