import singlecellmultiomics.features
import collections


class RNA_Flagger():

    def __init__(
            self,
            reference=None,
            alleleResolver=None,
            moleculeRadius=0,
            verbose=False,
            exon_gtf=None,
            intron_gtf=None,
            **kwargs):

        self.annotations = {}
        self.annotations['EX'] = singlecellmultiomics.features.FeatureContainer()
        self.annotations['EX'].loadGTF(exon_gtf, select_feature_type=['exon'])
        self.annotations['IN'] = singlecellmultiomics.features.FeatureContainer()
        self.annotations['IN'].loadGTF(
            intron_gtf, select_feature_type=['intron'])

        self.exon_hit_tag = 'eH'
        self.intron_hit_tag = 'iH'
        self.assinged_exons_tag = 'AE'
        self.assinged_introns_tag = 'AI'
        self.is_spliced_tag = 'SP'  # unused

        self.overlap_tag = 'FO'

    def digest(self, reads):

        feature_overlap = collections.Counter()  # feature->overlap
        # Go over reads and calculate overlap with features

        exon_count = 0
        intron_count = 0

        for read in reads:
            if read is None:
                continue

            states = ['.'] * read.query_length
            for q_pos, ref_pos in read.get_aligned_pairs(
                    matches_only=True, with_seq=False):

                overlaps_with_intron = False
                overlaps_with_exon = False
                exon_hits = set()

                for hit in self.annotations['EX'].findFeaturesAt(
                        chromosome=read.reference_name, lookupCoordinate=ref_pos, strand=None):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    overlaps_with_exon = True
                    exon_hits.add(hit_id)

                intron_hits = set()
                for hit in self.annotations['IN'].findFeaturesAt(
                        chromosome=read.reference_name, lookupCoordinate=ref_pos, strand=None):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    overlaps_with_intron = True
                    intron_hits.add(hit_id)

                if overlaps_with_exon and not overlaps_with_intron:
                    exon_count += 1
                    if self.overlap_tag is not None:
                        states[q_pos] = 'E'
                elif not overlaps_with_exon and overlaps_with_intron:
                    intron_count += 1
                    if self.overlap_tag is not None:
                        states[q_pos] = 'I'

            if self.overlap_tag is not None:
                read.set_tag(self.overlap_tag, ''.join(states))

        for read in reads:
            if read is not None:
                read.set_tag(self.exon_hit_tag, exon_count)
                read.set_tag(self.intron_hit_tag, intron_count)

                read.set_tag(
                    self.assinged_exons_tag, ','.join(
                        (str(e) for e in exon_hits)))
                read.set_tag(
                    self.assinged_introns_tag, ','.join(
                        (str(e) for e in intron_hits)))
