from singlecellmultiomics.fragment import Fragment

class NLAIIIFragment(Fragment):
    def __init__(self,reads, random_primer_length=6):
        Fragment.__init__(self, reads, assignment_radius=1_000, umi_hamming_distance=1 )

    def is_valid(self):
        # R1 needs to be mapped:
        if self.get_R1() is None:
            return False

        if self.get_R1().is_unmapped:
            return False
        
        return True
