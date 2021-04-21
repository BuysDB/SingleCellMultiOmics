import itertools
from singlecellmultiomics.molecule import Molecule
from singlecellmultiomics.molecule.nlaIII import NlaIIIMolecule
from singlecellmultiomics.molecule.nlaIII import AnnotatedNLAIIIMolecule
from singlecellmultiomics.molecule.chic import CHICMolecule
from singlecellmultiomics.molecule.chic import AnnotatedCHICMolecule
from collections import Counter
from singlecellmultiomics.utils.sequtils import complement
from itertools import product
from matplotlib.pyplot import get_cmap
from copy import copy
import numpy as np
complement_trans = str.maketrans('ATGC', 'TACG')


class TAPS:
    # Methylated Cs get converted into T readout
    def __init__(self, reference=None, reference_variants=None, taps_strand='F', **kwargs):
        """
        Intialise the TAPS class

        Args:
            reference (pysam.FastaFile) : reference fasta file

            reference_variants(pysam.VariantFile) : variants in comparison to supplied reference. @todo

        """
        assert reference is None, 'reference is now supplied by the molecule'
        if reference_variants is not None:
            raise NotImplementedError()

        self.overlap_tag = 'XM'

        """
        z unmethylated C in CpG context (CG)
        Z methylated C in CpG context (CG)
        x unmethylated C in CHG context ( C[ACT]G )
        X methylated C in CHG context   ( C[ACT]G )
        h unmethylated C in CHH context ( C[ACT][ACT] )
        H methylated C in CHH context ( C[ACT][ACT] )
        """
        self.context_mapping = dict()
        self.context_mapping[False] = {
            'CGA': 'z',
            'CGT': 'z',
            'CGC': 'z',
            'CGG': 'z',
            'CAG': 'x',
            'CCG': 'x',
            'CTG': 'x',
        }

        self.context_mapping[True] = {
            context: letter.upper() for context,
            letter in self.context_mapping[False].items()}

        self.taps_strand = taps_strand

        for x in list(itertools.product('ACT', repeat=2)):
            self.context_mapping[True][''.join(['C'] + list(x))] = 'H'
            self.context_mapping[False][''.join(['C'] + list(x))] = 'h'

        self.colormap = copy(get_cmap('RdYlBu_r')) # Make a copy
        self.colormap.set_bad((0,0,0)) # For reads without C's

    def position_to_context(
            self,
            chromosome,
            position,
            ref_base,
            observed_base='N',
            strand=None,
            reference=None
        ):
        """Extract bismark call letter from a chromosomal location given the observed base

        Args:

            chromosome (str) : chromosome / contig of location to test

            position (int) : genomic location of location to test

            observed_base(str) : base observed in the read

            strand(int) : mapping strand, False: forward, True: Reverse

        Returns:
            context(str) : 3 basepair context
            bismark_letter(str) : bismark call
        """

        assert reference is not None
        assert strand is not None

        qbase = observed_base.upper()

        #ref_base = reference.fetch( chromosome, position, position + 1).upper()


        # if a vcf file is supplied  we can extract the possible reference bases
        # @todo

        context = None
        methylated = None

        try:
            if ref_base == 'C':
                context = reference.fetch(chromosome, position, position + 3).upper()
                if qbase == 'T':
                    methylated = True
                elif qbase=='C':
                    methylated = False

            elif ref_base == 'G':
                origin = reference.fetch( chromosome, position - 2, position + 1).upper()
                context = origin.translate(complement_trans)[::-1]
                if qbase == 'A':
                    methylated = True
                elif qbase == 'G':
                    methylated = False

            else:
                raise ValueError('Only supply reference C or G')
        except ValueError: # Happens when the coordinates are outstide the reference:
            context = None
            methylated = None
        except FileNotFoundError:
            raise ValueError('Got a file not found error. This is probably caused by the reference handle being shared between processes. Stop doing that (quick solution: disable multithreading/multiprocess).')
        #print(ref_base, qbase,  strand, position, chromosome,context, '?' if methylated is None  else ('methylated' if methylated else  'not methylated'))

        if methylated is None:
            symbol='.'
        else:
            symbol = self.context_mapping[methylated].get(context, '.')
        return context, symbol

    def molecule_to_context_call_dict(self, molecule):
        """Extract bismark call_string dictionary from a molecule

        Args:
            molecule(singlecellmultiomics.molecule.Molecule) : molecule to extract calls from

        Returns:
            context_call_dict(dict) : {(chrom,pos):bismark call letter}
        """
        call_dict = {}  # (chrom,pos)-> "bismark" call_string
        for (chrom, pos), base in molecule.get_consensus().items():
            context, letter = self.position_to_context(chromosome=molecule.chromosome,
                                                       position=pos,
                                                       reference=molecule.reference,
                                                       observed_base=base,
                                                       strand=molecule.strand)
            if letter is not None:
                call_dict[(chrom, pos)] = letter
        return call_dict


class TAPSMolecule(Molecule):
    def __init__(self, fragments=None, taps=None, classifier=None, taps_strand='F', allow_unsafe_base_calls=False, methylation_consensus_kwargs=None, **kwargs):
        """ TAPSMolecule

        Args:
            fragments(list) :  list of fragments to associate with the Molecule
            taps(singlecellmultiomics.molecule.TAPS) : TAPS class which performs the methylation calling
            classifier : fitted sklearn classifier, when supplied this classifier is used to obtain a consensus from which the methylation calls are generated.

        """
        Molecule.__init__(self, fragments=fragments, **kwargs)
        if taps is None:
            raise ValueError("""Supply initialised TAPS class
                taps = singlecellmultiomics.molecule.TAPS( )
            """)
        self.taps = taps  # initialised TAPS class**self.
        self.methylation_call_dict = None
        self.classifier = classifier
        self.taps_strand = taps_strand
        self.allow_unsafe_base_calls = allow_unsafe_base_calls
        self.get_consensus_dictionaries_kwargs = methylation_consensus_kwargs if methylation_consensus_kwargs is not None else dict()

    def add_cpg_color_tag_to_read(self, read):
        try:
            CpG_methylation_rate = read.get_tag('sZ') / (read.get_tag('sZ') + read.get_tag('sz'))
            cfloat = self.taps.colormap(CpG_methylation_rate)[:3]
        except Exception as e:
            CpG_methylation_rate = None
            cfloat = self.taps.colormap._rgba_bad[:3]
        read.set_tag('YC', '%s,%s,%s' % tuple((int(x * 255) for x in cfloat)))

    def __finalise__(self):
        super().__finalise__()

        self.obtain_methylation_calls(**self.get_consensus_dictionaries_kwargs)

        for read in self.iter_reads():
            try:
                CpG_methylation_rate = read.get_tag('sZ')/(read.get_tag('sZ')+read.get_tag('sz'))
                cfloat = self.taps.colormap(CpG_methylation_rate)[:3]
            except Exception as e:
                CpG_methylation_rate = None
                cfloat = self.taps.colormap._rgba_bad[:3]
            read.set_tag('YC', '%s,%s,%s' % tuple( (int(x*255) for x in cfloat)))

    def write_tags_to_psuedoreads(self, reads, call_super=True):
        if call_super:
            Molecule.write_tags_to_psuedoreads(self,reads)
        for read in reads:
            self.add_cpg_color_tag_to_read(read)


    def is_valid(self, set_rejection_reasons=False):
        if not super().is_valid(set_rejection_reasons=set_rejection_reasons):
            return False
        """
        try:
            consensus = self.get_consensus(allow_unsafe=self.allow_unsafe_base_calls)
        except ValueError:
            if set_rejection_reasons:
                self.set_rejection_reason('no_consensus')
            return False

        except TypeError:
            if set_rejection_reasons:
                self.set_rejection_reason('getPairGenomicLocations_failed')
            return False
        """
        if self.methylation_call_dict is None:
            if set_rejection_reasons:
                self.set_rejection_reason('methylation_calls_failed')
            return False

        return True


    #def obtain_methylation_calls_experimental(self):




    def obtain_methylation_calls(self, **get_consensus_dictionaries_kwargs):
        """ This methods returns a methylation call dictionary

            returns:
                mcalls(dict) : (chrom,pos) : {'consensus': consensusBase, 'reference': referenceBase, 'call': call}
        """
        expected_base_to_be_converted = ('G' if self.strand else 'C') \
            if self.taps_strand=='F' else ('C' if self.strand else 'G')

        try:
            if True:
                c_pos_consensus, phreds, coverage = self.get_consensus(dove_safe = not self.allow_unsafe_base_calls,
                                                     only_include_refbase=expected_base_to_be_converted,
                                                     with_probs_and_obs=True, # Include phred scores and coverage
                                                     **get_consensus_dictionaries_kwargs)
            else:
                c_pos_consensus = self.get_consensus(dove_safe = not self.allow_unsafe_base_calls,
                                                 only_include_refbase=expected_base_to_be_converted,
                                                 with_probs_and_obs=False, # Include phred scores and coverage
                                                 **get_consensus_dictionaries_kwargs)

        except ValueError as e:
            if 'MD tag not present' in str(e):
                self.set_rejection_reason("MD_TAG_MISSING")
                return None
            raise

        # obtain the context of the conversions:
        conversion_contexts = {
            (contig, position):
            {'consensus': base_call,
             'reference_base': expected_base_to_be_converted,
             'cov': len(phreds[(contig, position)][base_call]),
             'qual': np.mean( phreds[(contig, position)][base_call] ),
             'context': self.taps.position_to_context(
                chromosome=contig,
                position=position,
                reference=self.reference,
                observed_base=base_call,
                ref_base=expected_base_to_be_converted,
                strand=self.strand)[1]}
            for (contig, position), base_call in c_pos_consensus.items()}

        # Write bismark tags:
        self.set_methylation_call_tags(conversion_contexts)
        return conversion_contexts


class TAPSNlaIIIMolecule(NlaIIIMolecule, TAPSMolecule):
    """Molecule class for combined TAPS and NLAIII """

    def __init__(self, fragments=None, taps=None,  taps_strand='R', **kwargs):
        NlaIIIMolecule.__init__(self, fragments, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps_strand=taps_strand, taps=taps, **kwargs)

    def write_tags(self):
        NlaIIIMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def __finalise__(self):
        NlaIIIMolecule.__finalise__(self)
        TAPSMolecule.__finalise__(self)

    def is_valid(self, set_rejection_reasons=False):
        return NlaIIIMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons) and TAPSMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons)

    def write_tags_to_psuedoreads(self, reads):
        NlaIIIMolecule.write_tags_to_psuedoreads(self, reads)
        TAPSMolecule.write_tags_to_psuedoreads(self, reads, call_super=False)


class AnnotatedTAPSNlaIIIMolecule(AnnotatedNLAIIIMolecule, TAPSMolecule):
    """Molecule class for combined TAPS, NLAIII and transcriptome """

    def __init__(self, fragments=None, features=None, taps=None, **kwargs):
        assert features is not None, "Supply features!"
        AnnotatedNLAIIIMolecule.__init__(
            self, fragments, features=features, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps=taps, **kwargs)

    def write_tags(self):
        AnnotatedNLAIIIMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self, set_rejection_reasons=False):
        return AnnotatedNLAIIIMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons) and TAPSMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons)

    def write_tags_to_psuedoreads(self, reads):
        AnnotatedNLAIIIMolecule.write_tags_to_psuedoreads(self, reads)
        TAPSMolecule.write_tags_to_psuedoreads(self, reads, call_super=False)

class TAPSCHICMolecule(CHICMolecule, TAPSMolecule):
    """Molecule class for combined TAPS and CHIC """

    def __init__(self, fragments=None, taps=None,  taps_strand='R', **kwargs):
        CHICMolecule.__init__(self, fragments, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps_strand=taps_strand,taps=taps, **kwargs)

    def write_tags(self):
        CHICMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def is_valid(self, set_rejection_reasons=False):
        return CHICMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons) and TAPSMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons)

    def __finalise__(self):
        CHICMolecule.__finalise__(self)
        TAPSMolecule.__finalise__(self)

    def write_tags_to_psuedoreads(self, reads):
        CHICMolecule.write_tags_to_psuedoreads(self, reads)
        TAPSMolecule.write_tags_to_psuedoreads(self, reads, call_super=False)


class AnnotatedTAPSCHICMolecule(AnnotatedCHICMolecule, TAPSMolecule):
    """Molecule class for combined TAPS, CHIC and transcriptome """

    def __init__(self, fragments=None, features=None, taps_strand='R', taps=None, **kwargs):
        assert features is not None, "Supply features!"
        AnnotatedCHICMolecule.__init__(
            self, fragments, features=features, **kwargs)
        TAPSMolecule.__init__(self, fragments=fragments, taps_strand=taps_strand,taps=taps, **kwargs)

    def write_tags(self):
        AnnotatedCHICMolecule.write_tags(self)
        TAPSMolecule.write_tags(self)

    def __finalise__(self):
        AnnotatedCHICMolecule.__finalise__(self)
        TAPSMolecule.__finalise__(self)

    def is_valid(self, set_rejection_reasons=False):
        return AnnotatedCHICMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons) and TAPSMolecule.is_valid(
            self, set_rejection_reasons=set_rejection_reasons)

    def write_tags_to_psuedoreads(self, reads):
        AnnotatedCHICMolecule.write_tags_to_psuedoreads(self, reads)
        TAPSMolecule.write_tags_to_psuedoreads(self, reads, call_super=False)


def strip_array_tags(molecule):
    for read in molecule.iter_reads():
        read.set_tag('jM',None)
        read.set_tag('jI',None)
    return molecule



idmap = [''.join(p) for p in product('0123','0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')]
context_ids = {''.join(ctx):idmap[i] for i,ctx in enumerate(product('ACTG',repeat=3))}

class TAPSPTaggedMolecule(AnnotatedTAPSNlaIIIMolecule):


    def __finalise__(self):
        AnnotatedTAPSNlaIIIMolecule.__finalise__(self)

        # Additional finalisation steps:
        conversions_count = 0
        forward_counts = 0
        rev_counts = 0
        conversion_locations = []
        mismatches = 0
        ambig_count = 0
        molecule_variants=set()

        # Obtain consensus call for the molecule
        consensus = self.get_consensus(allow_unsafe=True)

        context_obs = Counter()

        for read in self.iter_reads():
            for qpos, refpos, reference_base in read.get_aligned_pairs(with_seq=True, matches_only=True):

                location = (read.reference_name, refpos)

                #if (location[0], location[1]) in variant_locations:
                #    skipped+=1
                #    continue

                query = read.seq[qpos]
                context = self.reference.fetch(self.chromosome, refpos-1, refpos+2).upper()

                # Check if call is the same as the consensus call:
                if query!=consensus.get((read.reference_name,refpos), 'X'):
                    continue

                reference_base = reference_base.upper()
                if read.is_reverse: # forward, we sequence the other side
                    reference_base = complement(reference_base)
                    query = complement(query)
                    context=complement(context)

                if query!=reference_base:
                    mismatches += 1
                    molecule_variants.add(refpos)

                if (reference_base, query) in ( ('A','G'), ('G','A'), ('T','C'), ('C','T')):
                    conversion_locations.append(f'{location[0]}:{location[1]+1} {reference_base}>{query}' )
                    conversions_count+=1
                    context_obs[context_ids.get(context,'99')]+=1
                    if (reference_base, query) in ( ('A','G'), ('G','A') ):
                        forward_counts+=1
                    elif (reference_base, query) in (('T','C'), ('C','T')):
                        rev_counts+=1

        self.set_meta('mM',mismatches)
        self.set_meta('pP',conversions_count)
        self.set_meta('pR',rev_counts)
        self.set_meta('pF',forward_counts)
        self.set_meta('pE',','.join(conversion_locations))
        for context,obs in context_obs.most_common():
            self.set_meta(context,obs)

        self.mismatches = mismatches
        self.conversions_count = conversions_count
        self.context_obs = context_obs
        self.rev_counts = rev_counts
        self.forward_counts = forward_counts
        self.conversion_locations=conversion_locations

    def write_tags_to_psuedoreads(self, reads):
        AnnotatedTAPSNlaIIIMolecule.write_tags_to_psuedoreads(self,reads)
