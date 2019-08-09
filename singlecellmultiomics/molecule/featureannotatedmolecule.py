from singlecellmultiomics.molecule.molecule import Molecule
import collections



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
        self.is_spliced=None

    def set_spliced(self, is_spliced):
        """ Set wether the transcript is spliced, False has priority over True """
        if self.is_spliced == True and not is_spliced:
            # has already been set
            self.is_spliced = False
        else:
            self.is_spliced = is_spliced

    def set_intron_exon_features(self):
        if not self.is_annotated:
            self.annotate()

        # Collect all hits:
        hits = self.hits.keys()

        # (gene, transcript) -> set( exon_id  .. )
        exon_hits = collections.defaultdict(set)
        intron_hits = collections.defaultdict(set)

        for hit, locations in self.hits.items():
            if type(hit) is not tuple:
                continue
            meta = dict(list(hit))
            if not 'gene_id' in meta:
                continue
            if meta.get('type')=='exon':
                if not 'transcript_id' in meta:
                    continue
                self.genes.add(meta['gene_id'])
                self.exons.add(meta['exon_id'])
                exon_hits[(meta['gene_id'],meta['transcript_id'])].add(meta['exon_id'])
            elif meta.get('type')=='intron':
                self.genes.add(meta['gene_id'])
                self.introns.add(meta['gene_id'])
                exon_hits[(meta['gene_id'],meta['transcript_id'])].add(meta['exon_id'])

        # Find junctions and add all annotations to annotation sets
        debug = []

        for (gene, transcript), exons_overlapping in exon_hits.items():
            # If two exons are detected from the same gene we detected a junction:
            if len(exons_overlapping)>1:
                self.junctions.add(gene)

                # We found two exons and an intron:
                if gene in self.introns:
                    self.set_spliced(False)
                else:
                    self.set_spliced(True)

            debug.append( f'{gene}_{transcript}:' + ','.join(list(exons_overlapping)) )

        # Write exon dictionary:
        self.set_meta('DB', ';'.join(debug) )

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
        if self.is_spliced is True:
            self.set_meta('SP',True)
        elif self.is_spliced is False:
            self.set_meta('SP',False)

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
                for start,end in self.get_aligned_blocks():
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
