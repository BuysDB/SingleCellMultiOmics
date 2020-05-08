from singlecellmultiomics.molecule.molecule import Molecule
from singlecellmultiomics.utils.sequtils import phred_to_prob
import numpy as np

class ScarTraceMolecule(Molecule):
    """ScarTrace Molecule class

    Args:
        fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule

        **kwargs: extra args

    """

    def __init__(self, fragment,
                 **kwargs):
        Molecule.__init__(self, fragment, **kwargs)

        # Run code to obtain scar description
        assert len(self) == 1, 'only implemented for molecules with a single associated fragment'

        scarDescription = set()
        qualities = []
        for read in self.iter_reads():
            if read is None:
                continue
            if read.is_unmapped:
                continue
            firstCigarOperation, firstCigarLen = read.cigartuples[0]
            insertPos = 0
            lastNonNoneRefPos = read.reference_start if firstCigarOperation != 4 else read.reference_start - firstCigarLen

            expandedCigar = []
            for cigarOperation, cigarLen in read.cigartuples:
                expandedCigar += [cigarOperation] * cigarLen

            for queryPos, referencePos in read.get_aligned_pairs(
                    matches_only=False):
                if queryPos is not None:
                    qualities.append(phred_to_prob(read.qual[queryPos]))

                if queryPos is None and referencePos is None:
                    continue

                if referencePos is not None:
                    lastNonNoneRefPos = referencePos
                    insertPos = 0
                    if queryPos is not None:  # If both the reference and query match, there is not scar information
                        continue

                if queryPos is None:
                    scarDescription.add(f'{referencePos}.D')
                elif referencePos is None:  # insert or clip:
                    operation = expandedCigar[queryPos]
                    if operation == 1:
                        if lastNonNoneRefPos is None:
                            raise ValueError('Unsolvable :(')
                        queryBase = read.seq[queryPos]
                        scarDescription.add(
                            f'{queryBase}.{insertPos+lastNonNoneRefPos}.I')
                        insertPos += 1

        scarDescription = ','.join(sorted(list(scarDescription)))

        # Add average base calling quality excluding primers:
        meanQ = np.mean(qualities) if len(qualities)>0 else 0
        for read in self.iter_reads():
            read.set_tag('SQ', 1 - meanQ)
            if len(scarDescription) == 0:
                scarDescription = 'WT'
            #@todo when collapsing into molecule:
            read.set_tag('SD', scarDescription)
            #@todo when collapsing into molecule:
        if self[0].get_R1() is not None:
            self.set_meta('DS', self[0].get_R1().reference_start)
            self.set_meta('RZ', 'SCAR')
            self.set_meta('DT', 'SCAR')


    def write_tags(self):
        Molecule.write_tags(self)

    def is_valid(self, set_rejection_reasons=False, reference=None):
        if not any( fragment.is_valid() for fragment in self ):
            self.set_rejection_reason('all_fragments_rejected')
            return False

        return True
