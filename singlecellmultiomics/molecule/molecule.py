from singlecellmultiomics.utils.sequtils import hamming_distance

class Molecule():
    def __init__(self, fragments=None, assignment_radius=3, umi_hamming_distance=1 ):
        self.assignment_radius = assignment_radius
        self.umi_hamming_distance = umi_hamming_distance
        self.fragments  = []

        self.molecule_id = None

        if fragments is not None:
            if type(fragments) is list:
                for frag in fragments:
                    self.add_fragment(frag)
            else:
                self.add_fragment(fragments)


    def __repr__(self):
        frag_repr = '\n\t'.join([', '.join([str(read) for read in fragment]) for fragment in self.fragments])
        return f"""Molecule
        with {len(self.fragments)} assinged fragments
        id: {str(self.molecule_id)}\n
        """ + frag_repr

    def assignment_function(self, fragment):
        return  fragment[0].get_tag('SM'),fragment[0].get_tag('RX'),fragment[0].get_tag('RS'),int(fragment[0].get_tag('DS'))

    # Returns True if two molecule ids  are identical or close enough
    def eq_function(self, assignment_a, assignment_b):
        if assignment_a is None or assignment_b is None:
            return False

        sample_A, umi_A, strand_A, location_A = assignment_a
        sample_B, umi_B, strand_B, location_B = assignment_b

        # Make sure fragments map to the same strand, cheap comparisons
        if sample_A!=sample_B or strand_A!=strand_B:
            return False

        # Make sure fragments map close enough to eachother, cheap comparison
        if abs(location_A-location_B)>self.assignment_radius:
            return False

        # Make sure UMI's are similar enough, more expensive hamming distance calculation
        if self.umi_hamming_distance==0:
            return umi_A==umi_B
        else:
            return hamming_distance(umi_A,umi_B)<=self.umi_hamming_distance

    def fragment_eq(self, fragment, fragment_molecule_id=None): # can a fragment be added to the molecule?
        if self.molecule_id is None:
            return True
        if fragment_molecule_id is None:
            fragment_molecule_id = self.assignment_function(fragment)
        return self.eq_function(fragment_molecule_id, self.molecule_id)


    def add_fragment(self, fragment):
        fragment_molecule_id = self.assignment_function(fragment)
        if self.fragment_eq(fragment,fragment_molecule_id):
            self.molecule_id = fragment_molecule_id
            self.fragments.append(fragment)
            return True
        else:
            return False

    def _localisation_function(self, fragment):
        if not fragment[0].has_tag('DS'):
            return None
        return fragment[0].get_tag('DS')

    def sample_assignment_function(self, fragment):
        for read in fragment:
            if read is not None:
                if read.has_tag(self.sample_tag):
                    return read.get_tag(self.sample_tag)
        return None
