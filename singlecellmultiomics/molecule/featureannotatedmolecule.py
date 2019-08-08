from singlecellmultiomics.molecule.molecule import Molecule
import collections
from more_itertools import consecutive_groups

# https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]

class FeatureAnnotatedMolecule(Molecule):
    """Molecule which is annotated with features (genes/exons/introns, .. )
    """

    def __init__(self, fragment, features, stranded=None, **kwargs):
        """
            Args:
                fragments (singlecellmultiomics.fragment.Fragment): Fragments to associate to the molecule
                features (singlecellmultiomics.features.FeatureContainer) : container to use to obtain features from
                stranded : None; not stranded, False: same strand as R1, True: other strand
                **kwargs: extra args

        """
        Molecule.__init__(self,fragment,**kwargs)
        self.features = features
        self.hits = collections.defaultdict(set) #feature -> hit_bases
        self.stranded = stranded
        self.is_annotated = False

        self.junctions = set()
        self.genes = set()
        self.introns = set()
        self.exons = set()

    def set_intron_exon_features(self):
        if not self.is_annotated:
            self.annotate()

        # Collect all hits:
        hits = self.hits.keys()
        f_hits = collections.defaultdict(collections.Counter)
        for hit in hits:
            if hit.startswith('type:exon'):
                gene = hit.split(',')[-1].replace('gene_id:','')
                self.genes.add(gene)
                self.exons.add(gene)
                #if allele is not None:
                #    gene = f'{allele}_{gene}_{self.chromosome}'
                f_hits[gene]['exon']+=1
            elif hit.startswith('type:intron'):
                gene = hit.split(',')[-1].replace('gene_id:','')
                self.genes.add(gene)
                self.introns.add(gene)
                #if allele is not None:
                #       gene = f'{allele}_{gene}_{self.chromosome}'
                f_hits[gene]['intron']+=1

        # Find junctions and add all annotations to annotation sets
        for gene, intron_exon_hits in f_hits.items():

            spliced=True
            if 'intron' in intron_exon_hits:
                spliced=False

            # If two exons are detected from the same gene we detected a junction:
            if intron_exon_hits['exon']>=2:
                self.junctions.add(gene)


    def write_tags(self):
        Molecule.write_tags(self)

        if len(self.exons)>0:
            self.set_meta('EX',  ','.join(sorted([ str(x) for x in self.exons] )))
        if len(self.introns)>0:
            self.set_meta('IN',  ','.join(sorted([ str(x) for x in self.introns] )))
        if len(self.genes)>0:
            self.set_meta('GN',  ','.join(sorted([ str(x) for x in self.genes] )))
        if  len(self.junctions)>0:
            self.set_meta('JN',  ','.join(sorted([ str(x) for x in self.junctions] )))
            # Is transcriptome
            self.set_meta('IT',1)
        elif len(self.genes)>0:
            # Maps to gene but not junction
            self.set_meta('IT',0.5)
        else:
            # Doesn't map to gene
            self.set_meta('IT',0)


    def annotate(self, method=0):
        """
            Args:
                method (int) : 0, obtain blocks and then obtain features. 1, try to obtain features for every aligned base

        """
        # When self.stranded is None, set to None strand. If self.stranded is True reverse the strand, otherwise copy the strand
        strand = None if self.stranded is None else '+-'[( not self.strand if self.stranded else self.strand )]
        self.is_annotated = True
        if method==0:

            # Obtain all blocks:
            try:
                for start,end in find_ranges(
                    sorted(list(set(
                        (ref_pos
                        for read in self.iter_reads()
                        for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False) ))))
                    ):
                    for hit in self.features.findFeaturesBetween(chromosome =self.chromosome, sampleStart=start, sampleEnd=end, strand=strand):
                        hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                        self.hits[hit_ids].add((self.chromosome,(hit_start,hit_end)))

            except TypeError:
                # This happens when no reads map
                pass
        else:

            for read in self.iter_reads():
                for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                    for hit in self.features.findFeaturesAt(chromosome=read.reference_name,lookupCoordinate=ref_pos,strand=strand):
                        hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                        self.hits[hit_ids].add((read.reference_name,ref_pos))
