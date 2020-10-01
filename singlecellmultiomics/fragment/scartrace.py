from singlecellmultiomics.utils.sequtils import reverse_complement
from singlecellmultiomics.fragment import Fragment

class ScarTraceFragment(Fragment):
    """
    Fragment definition for ScarTrace
    """

    def __init__(self, reads,scartrace_r1_primers=None, **kwargs):
        Fragment.__init__(self, reads,  **kwargs)
        self.scartrace_r1_primers = scartrace_r1_primers
        assert scartrace_r1_primers is not None, 'please supply primer sequences'

    # remove the set_umi function
    def set_umi(self, **kwargs):
        pass

    def is_valid(self):
        """ Check if R1 starts with the defined primers and if R2 and R1 are mapped

        Returns:
            bool
        """
        if not self.has_R1() or not self.has_R2():
            if not self.has_R1():
                self.set_meta('fr','no_R1', as_set=True)
            if not self.has_R2():
                self.set_meta('fr','no_R2', as_set=True)
            return False

        if self[0].is_unmapped or self[1].is_unmapped:
            self.set_meta('fr','unmapped_mate', as_set=True)
            return False

        r1_seq = self.get_R1().seq

        if self.get_R1().is_reverse:
            r1_seq = reverse_complement(r1_seq)

        if any( r1_seq.startswith(primer) for primer in self.scartrace_r1_primers ):
            # Good
            return True

        self.set_meta('fr','primer_not_matching', as_set=True)
        # Bad
        return False




    # Replace the equality function
    def __eq__(self, other):  # other can also be a Molecule!
        # Make sure fragments map to the same strand, cheap comparisons
        if self.sample != other.sample:
            return False

        if self.strand != other.strand:
            return False

        if min(abs(self.span[1] -
                   other.span[1]), abs(self.span[2] -
                                       other.span[2])) > self.assignment_radius:
            return False

        # Sample matches and starting position is within the defined span
        # radius
        return True
