import itertools
from singlecellmultiomics.molecule import Molecule
from singlecellmultiomics.molecule.nlaIII import NlaIIIMolecule
from singlecellmultiomics.molecule.nlaIII import AnnotatedNLAIIIMolecule
from singlecellmultiomics.molecule.chic import CHICMolecule
from singlecellmultiomics.molecule.chic import AnnotatedCHICMolecule

complement = str.maketrans('ATGC', 'TACG')

class TAPS():
    # Methylated Cs get converted into T readout
    def __init__(self, reference, reference_variants=None,**kwargs):
        """
        Intialise the TAPS class

        Args:
            reference (pysam.FastaFile) : reference fasta file

            reference_variants(pysam.VariantFile) : variants in comparison to supplied reference. @todo

        """
        if reference_variants is not None:
            raise NotImplementedError()

        self.overlap_tag = 'XM'
        self.reference = reference
        if self.reference is None:
            raise ValueError('A reference fasta file is required!')

        """
        z unmethylated C in CpG context (CG)
        Z methylated C in CpG context (CG)
        x unmethylated C in CHG context ( C[ACT]G )
        X methylated C in CHG context   ( C[ACT]G )
        h unmethylated C in CHH context ( C[ACT][ACT] )
        H methylated C in CHH context ( C[ACT][ACT] )
        """
        self.context_mapping={}
        self.context_mapping[False] = {
            'CGA':'z',
            'CGT':'z',
            'CGC':'z',
            'CGG':'z',
            'CAG':'x',
            'CCG':'x',
            'CTG':'x',
        }

        self.context_mapping[True] = { context:letter.upper()
                                      for context,letter in self.context_mapping[False].items() }


        for x in list( itertools.product('ACT',repeat=2) ):
            self.context_mapping[True][ ''.join( ['C']+list(x)) ] =  'H'
            self.context_mapping[False][ ''.join( ['C']+list(x)) ] =  'h'

    def position_to_context(self,chromosome, position, observed_base='N', strand=0):
        """Extract bismark call letter from a chromosomal location given the observed base

        Args:

            chromosome (str) : chromosome / contig of location to test

            position (int) : genomic location of location to test

            observed_base(str) : base observed in the read

            strand(str) : mapping strand

        Returns:
            context(str) : 3 basepair context
            bismark_letter(str) : bismark call
        """


        qbase = observed_base.upper()

        ref_base= self.reference.fetch(chromosome, position,position+1).upper()

        # if a vcf file is supplied  we can extract the possible reference bases
        # @todo

        context=None
        methylated = False
        rpos = position
        if ref_base=='C' and strand==0:
            context = self.reference.fetch(chromosome, rpos, rpos+3).upper()

            if qbase=='T':
                methylated=True
            methylationStateString = self.context_mapping[methylated].get(context,'uU'[methylated])

        elif ref_base=='G' and strand==1:
            origin = self.reference.fetch(chromosome, rpos-2, rpos+1).upper()
            context = origin.translate(complement)[::-1]
            if qbase=='A':
                methylated=True

        return context,self.context_mapping[methylated].get(context,'.')

    def molecule_to_context_call_dict(self, molecule):
        """Extract bismark call_string dictionary from a molecule

        Args:
            molecule(singlecellmultiomics.molecule.Molecule) : molecule to extract calls from

        Returns:
            context_call_dict(dict) : {(chrom,pos):bismark call letter}
        """
        call_dict = {} # (chrom,pos)-> "bismark" call_string
        for (chrom,pos),base in molecule.get_consensus().items():
            context, letter = self.position_to_context(chromosome=molecule.chromosome,
                                     position=pos,
                                     observed_base=base,
                                     strand=molecule.get_strand())
            if letter is not None:
                call_dict[(chrom,pos)] = letter
        return call_dict


class TAPSMolecule(Molecule):
    def __init__(self, fragments=None, taps=None, classifier=None,  **kwargs):
        """ TAPSMolecule

        Args:
            fragments(list) :  list of fragments to associate with the Molecule
            taps(singlecellmultiomics.molecule.TAPS) : TAPS class which performs the methylation calling
            classifier : fitted sklearn classifier, when supplied this classifier is used to obtain a consensus from which the methylation calls are generated.

        """
        Molecule.__init__(self, fragments=fragments, **kwargs)
        if taps is None:
            raise ValueError("""Supply initialised TAPS class
                taps = singlecellmultiomics.molecule.TAPS( reference )
            """)
        self.taps = taps #initialised TAPS class
        self.methylation_call_dict = None
        self.classifier = classifier
        try:
            self.obtain_methylation_calls(classifier=self.classifier)
        except ValueError:
            pass

    def is_valid(self,set_rejection_reasons=False):
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

        if self.methylation_call_dict is None:
            if set_rejection_reasons:
                self.set_rejection_reason('methylation_calls_failed')
            return False

        return True



    def obtain_methylation_calls(self, classifier=None):
        """ This methods returns a methylation call dictionary
            Args:
                classifier : classifier used for consensus determination
            returns:
                mcalls(dict) : (chrom,pos) : {'consensus': consensusBase, 'reference': referenceBase, 'call': call}
        """

        # Find all aligned positions and corresponding reference bases:
        aligned_reference_positions = {} #(chrom,pos)->base
        for read in self.iter_reads():
            for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True, matches_only=True):
                aligned_reference_positions[(read.reference_name,ref_pos)] = ref_base.upper()

        # Obtain consensus:
        try:
            consensus = self.get_consensus(classifier=classifier)
        except ValueError:
            raise ValueError('Cannot obtain a safe consensus for this molecule')


        # find all locations where a C/G was converted into A/T, now strand specific
        converted_bases = 0
        conversions = {}
        for location, reference_base in aligned_reference_positions.items():
            if not location in consensus:
                continue
            if (not self.strand and reference_base=='C' and consensus[location] in 'CT') or \
                self.strand and reference_base=='G' and consensus[location] in 'AG':
                conversions[location] = {'ref':reference_base, 'obs':consensus[location]}
                if consensus[location] in 'TA':
                    converted_bases+=1

        # obtain the context of the conversions:
        conversion_contexts  ={
                location:
                {'consensus': consensus[location],
                 'reference_base':conversions[location]['ref'],
                 'context':self.taps.position_to_context(
                        *location,
                        observed_base = observations['obs'],
                        strand = self.strand)[1]}
            for location, observations in conversions.items()}

        # Write bismark tags:
        self.set_methylation_call_tags(conversion_contexts)

class TAPSNlaIIIMolecule(NlaIIIMolecule,TAPSMolecule):
    """Molecule class for combined TAPS and NLAIII """
    def __init__(self, fragments=None, taps=None, **kwargs):
        NlaIIIMolecule.__init__(self, fragments, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps=taps, **kwargs)

    def write_tags(self):
        NlaIIIMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self,set_rejection_reasons=False):
        return NlaIIIMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons) and \
               TAPSMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons)


class AnnotatedTAPSNlaIIIMolecule(AnnotatedNLAIIIMolecule,TAPSMolecule):
    """Molecule class for combined TAPS, NLAIII and transcriptome """
    def __init__(self, fragments=None, features=None, taps=None, **kwargs):
        assert features is not None, "Supply features!"
        AnnotatedNLAIIIMolecule.__init__(self, fragments, features=features, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps=taps, **kwargs)

    def write_tags(self):
        AnnotatedNLAIIIMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self,set_rejection_reasons=False):
        return AnnotatedNLAIIIMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons) and \
               TAPSMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons)

class TAPSCHICMolecule(CHICMolecule,TAPSMolecule):
    """Molecule class for combined TAPS and CHIC """
    def __init__(self, fragments=None, taps=None, **kwargs):
        CHICMolecule.__init__(self, fragments, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps=taps, **kwargs)

    def write_tags(self):
        CHICMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self,set_rejection_reasons=False):
        return CHICMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons) and \
               TAPSMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons)


class AnnotatedTAPSCHICMolecule(AnnotatedCHICMolecule,TAPSMolecule):
    """Molecule class for combined TAPS, CHIC and transcriptome """
    def __init__(self, fragments=None, features=None, taps=None, **kwargs):
        assert features is not None, "Supply features!"
        AnnotatedCHICMolecule.__init__(self, fragments, features=features, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps=taps, **kwargs)

    def write_tags(self):
        AnnotatedCHICMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self,set_rejection_reasons=False):
        return AnnotatedCHICMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons) and \
               TAPSMolecule.is_valid(self,set_rejection_reasons=set_rejection_reasons)
