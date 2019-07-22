import itertools
complement = str.maketrans('ATGC', 'TACG')

class TAPS():
    # Methylated Cs get converted into T readout
    def __init__(self, reference, **kwargs):
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
