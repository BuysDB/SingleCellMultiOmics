#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import gzip
import itertools
import re
import functools
import pysam
from singlecellmultiomics.utils import Prefetcher
from copy import copy
import collections
import pandas as pd

def get_gene_id_to_gene_name_conversion_table(annotation_path_exons,
                                              featureTypes=['gene_name']):
    """Create a dictionary converting a gene id to other gene features,
        such as gene_name/gene_biotype etc.

    Arguments:
        annotation_path_exons(str) : path to GTF file (can be gzipped)
        featureTypes(list) : list of features to convert to, for example ['gene_name','gene_biotype']

    Returns:
        conversion_dict(dict) : { gene_id : 'firstFeature_secondFeature'}
    """
    conversion_table = {}
    with (gzip.open(annotation_path_exons, 'rt') if annotation_path_exons.endswith('.gz') else open(annotation_path_exons, 'r')) as t:
        for i, line in enumerate(t):
            parts = line.rstrip().split(None, 8)
            keyValues = {}
            for part in parts[-1].split(';'):
                kv = part.strip().split()
                if len(kv) == 2:
                    key = kv[0]
                    value = kv[1].replace('"', '')
                    keyValues[key] = value
            # determine the conversion name:
            if 'gene_id' in keyValues and any(
                    [feat in keyValues for feat in featureTypes]):
                conversion_table[keyValues['gene_id']] = '_'.join([
                    keyValues.get(feature, 'None')
                    for feature in featureTypes])

    return conversion_table


class FeatureContainer(Prefetcher):

    def __init__(self, verbose=False):
        self.args = locals().copy()
        self.startCoordinates = {}  # dict of np.array()
        self.features = {}
        self.endCoordinates = {}
        self.sorted = True  # Flag containing if the features are all coordinate sorted
        # When set to true, the class will be (very) verbose
        self.debug = False
        self.remapKeys = {}  # {'chrMT':'chrM','MT':'chrM'}  Use this to convert chromosome names between the GTF and requested locations
        self.verbose = verbose
        self.preload_list = []

    def debugMsg(self, msg):
        if self.verbose:
            print(msg)

    def __repr__(self):
        s= 'FeatureContainer,'
        if len(self.preload_list):
            s+= ' Preloaded files:\n'
            for f in self.preload_list:
                s+=str(f)+'\n'
        if len(self.features):
            s+=f'Features in memory for {len(self.features)} contigs\n'
        else:
            s+='No features in memory\n'
        return s

    def __len__(self):
        return sum(len(f) for f in self.features.values())

    def preload_GTF(self, **kwargs):
        self.preload_list.append( {'gtf':kwargs} )

    def get_gene_to_location_dict(self, meta_key='gene_name', with_strand=False):
        """
        generate dictionary, {gene_name: contig,start,end}

        Args:
            meta_key(str): key of the meta information used to use as primary key for the returned gene_locations

        Returns:
            gene_locations(dict)
        """
        location_map = {}

        for contig, start, end, name, strand, meta in self:
            meta =dict(meta)
            try:
                if with_strand:
                    location_map[ meta[meta_key]] = (contig, start,end,strand)
                else:
                    location_map[ meta[meta_key]] = (contig, start,end)
            except Exception as e:
                pass

        return location_map



    def __iter__(self):
        for contig, contig_features in self.features.items():
            for feature in contig_features:
                yield (contig,)+ feature


    def instance(self, arg_update):
        if 'self' in self.args:
            del self.args['self']
        clone = FeatureContainer(**self.args)
        for cmd in self.preload_list:
            for preload_type, kwargs in cmd.items():
                kwargs_copy = copy(kwargs)
                kwargs_copy.update(arg_update)
                if preload_type=='gtf':
                    clone.loadGTF(**kwargs_copy)
                else:
                    raise ValueError()
        return clone


    def prefetch(self, contig, start, end):
        return self.instance( {'contig':contig, 'region_start':start,  'region_end':end})


    def loadGTF(self, path, thirdOnly=None, identifierFields=['gene_id'],
                ignChr=False, select_feature_type=None, exon_select=None,
                head=None, store_all=False, contig=None, offset=-1,
                region_start=None, region_end=None):
        """Load annotations from a GTF file.
        ignChr: ignore the chr part of the Annotation chromosome
        """
        if region_end is not None or region_start is not None:
            assert contig is not None and region_end is not None and region_start is not None

        self.loadedGtfFeatures = thirdOnly
        #pattern = '^(.*) "(.*).*"'
        #prog = re.compile(pattern)
        if self.verbose:
            if contig is None:
                print(f"Loading {path} completely")
            elif region_start is None:
                print(f"Loading {path}, for contig {contig}")
            else:
                print(f"Loading {path}, for contig {contig}:{region_start}-{region_end}")
        added = 0
        with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')) as f:
            for line_id, line in enumerate(f):
                if head is not None and added > head:
                    break
                if line[0] != '#':
                    parts = line.rstrip().split(None, 8)

                    chrom = parts[0]
                    if contig is not None and chrom != contig:
                        continue

                    if thirdOnly is not None:
                        if parts[2] not in thirdOnly:
                            continue

                    # Example line: (gene)
                    # 1  havana  gene    11869   14409   .   +   .   gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; havana_gene "OTTHUMG00000000961"; havana_gene_version "2";
                    # Example line (exon)
                    # 1  havana  exon    11869   12227   .   +   .   gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; havana_gene "OTTHUMG00000000961"; havana_gene_version "2"; transcript_name "DDX11L1-002"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTHUMT00000362751"; havana_transcript_version "1"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
                    #Part, info
                    # 0  Chromosome
                    # 1  Source
                    # 2  Feature type (matches thirdOnly)
                    # 3 Feature Start
                    # 4 Feature End
                    # 5 .
                    # 6 Strand
                    # 7 .
                    # 8 KEY "VALUE";
                    #keyValues = { part.strip().split()[0]:part.strip().split()[1].replace('"','') for part in parts[-1].split(';') }

                    #keyValues = {i.group(1) : i.group(2) for i in (prog.match(j) for j in parts[-1].split('; '))}
                    if select_feature_type is not None and not parts[2] in select_feature_type:
                        continue

                    exon = parts[7]
                    if exon_select is not None and exon not in exon_select:
                        continue

                    keyValues = {}
                    for part in parts[-1].split(';'):
                        kv = part.strip().split()
                        if len(kv) == 2:
                            key = kv[0]
                            value = kv[1].replace('"', '')
                            keyValues[key] = value

                    #self.addFeature( self.remapKeys.get(parts[0], parts[0]), int(parts[3]), int(parts[4]), parts[9].replace('"','').replace(';',''))

                    chrom = self.remapKeys.get(chrom, chrom)
                    chromosome = chrom if ignChr == False else chrom.replace(
                        'chr', '')

                    if identifierFields is None:
                        if parts[2] == 'exon':
                            featureName = keyValues['exon_id']
                            #featureName = ','.join([keyValues['exon_id'],keyValues['transcript_id']])
                        elif parts[2] == 'gene':
                            featureName = keyValues['gene_id']
                        elif parts[2] == 'transcript':
                            featureName = keyValues['transcript_id']
                        else:
                            featureName = ','.join(
                                [parts[2], parts[3], parts[4], keyValues['transcript_id']])
                    else:
                        featureName = ','.join(
                            [keyValues.get(i, 'none') for i in identifierFields if i in keyValues])


                    start = int( parts[3] ) + offset
                    end = int( parts[4] ) + offset

                    if region_end is not None and region_start is not None and ( end<region_start or start>region_end):
                        continue

                    if store_all:
                        keyValues['type'] = parts[2]
                        self.addFeature(
                            self.remapKeys.get(
                                chromosome, chromosome),start,end, strand=parts[6], name=featureName, data=tuple(
                                keyValues.items()))

                    else:
                        self.addFeature(
                            self.remapKeys.get(
                                chromosome, chromosome), start,end, strand=parts[6], name=featureName, data=','.join(
                                (':'.join(
                                    ('type', parts[2])), ':'.join(
                                    ('gene_id', keyValues['gene_id'])))))
                    added += 1

            if self.verbose:
                print("Loaded %s features, now sorting" %
                      sum([len(self.features[c]) for c in self.features]))
            self.sort()
            if self.verbose:
                print("done sorting")
        if self.verbose:
            print("The following chromosomes are available:")
            print(', '.join(sorted(list(self.startCoordinates.keys()))))

    def annotateUTRs(self, utrs=['three_prime_utr', 'five_prime_utr']):
        """flag the exons that contain a utr"""

        chromosomes = self.features.keys()
        typeRegex = re.compile('.*type:([^,]*).*')

        for c in chromosomes:
            print('chromosome:%s' % (c))
            types = np.array([typeRegex.match(f[-1]).group(1)
                              for f in self.features[c]])
            featStartEndStrand = np.array([[f[0] for f in self.features[c]], [
                                          f[1] for f in self.features[c]], [f[3] for f in self.features[c]]])

            fe = np.where(types == 'exon')[0]
            names = np.array([f[2] for f in self.features[c]])[fe]
            starts = featStartEndStrand[0, :][fe]

            for utr in utrs:
                isUtrF = 'is_' + utr + ':False'
                isUtrT = 'is_' + utr + ':True'
                utrRegex = re.compile(isUtrF)

                for position in fe:
                    (hitStart, hitEnd, name, hitStrand,
                     data) = self.features[c][position]
                    data = ','.join((data, isUtrF))
                    self.features[c][position] = (
                        hitStart, hitEnd, name, hitStrand, data)
                lu = np.where(types == utr)[0]
                if lu.size:
                    foo, idx = np.unique(
                        featStartEndStrand[:, lu], axis=1, return_index=True)
                    lu = lu[idx]

                    for i in lu:
                        hits = self.findFeaturesBetween(
                            c, self.features[c][i][0], self.features[c][i][0], self.features[c][i][3])
                        hitNames = np.array([h[2] for h in hits])
                        hitStarts = np.array([h[0] for h in hits])
                        for index, n in enumerate(hitNames):
                            fl = starts == hitStarts[index]
                            fn = names[fl] == n
                            positions = fe[fl][fn]
                            for position in positions:
                                (hitStart, hitEnd, name, hitStrand,
                                 data) = self.features[c][position]
                                data = utrRegex.sub(isUtrT, data)
                                self.features[c][position] = (
                                    hitStart, hitEnd, name, hitStrand, data)

                        # hitTypes = np.array([typeRegex.match(h[-1]).group(1) for h in hits])#np.array([h[-1]['type'] for h in hits])
                        #hitNames = np.array([h[2] for h in hits])
                        #fh = hitTypes=='exon'

                        # if np.any(fh):
                        #   hitNames=np.unique(hitNames[fh])
                        #   for n in hitNames:
                        #       fn = names==n
                        #       positions = fe[fn]
                        #       for position in positions:
                        #           (hitStart, hitEnd, name, hitStrand, data) = self.features[c][position]
                        #           data = utrRegex.sub(isUtrT,data)
                        #           self.features[c][position] = (hitStart, hitEnd, name, hitStrand, data)

    def loadBED(self, path, ignChr=False, parseBlocks=True):
        """Load UCSC based table. """

        """
        parseBlocks ; parse the defined blocks as separate features
        """
        with (gzip.open(path, 'rt') if '.gz' in path else open(path, 'r')) as f:
            # chr1   14662   187832  calJac3:ENSCJAT00000061222.1-1.1    943 -
            # 187832  187832  0   9   5,129,69,110,42,23,133,202,78,
            # 0,38,307,1133,1243,1944,1970,172713,173092,
            for line in f:
                if line.split()[0] == "track":
                    continue
                blockAvail = False

                parts = line.strip().split()
                strand = None
                name = None
                score = None
                if len(parts)==12:
                    chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRGB, blockCount, blockSizes, blockStarts = parts
                    blockAvail = True
                elif len(parts)==10:
                    chrom, chromStart, chromEnd, name, score, strand, itemRGB, blockCount, blockSizes, blockStarts = parts
                    blockAvail = True
                elif len(parts)==6:
                    chrom, chromStart, chromEnd, name, score, strand = parts
                elif len(parts)>=4:
                    chrom, chromStart, chromEnd, value = parts[:4]
                    name = value
                else:
                    raise ValueError('Could not read the supplied bed file, it has too little columns, expecting at least 4 columns: contig,start,end,value')

                chrom = self.remapKeys.get(chrom, chrom)
                chrom = chrom if ignChr == False else chrom.replace('chr', '')

                if blockAvail and parseBlocks:
                    blockStarts = blockStarts.split(',')
                    blocks = [
                        (int(
                            blockStarts[index]) +
                            int(chromStart),
                            int(blockSize)) for index,
                        blockSize in enumerate(
                            blockSizes.split(',')) if len(blockSize) > 0]
                    if len(blocks) != int(blockCount):
                        raise ValueError(
                            'BlockCount at line %s do not match actual amount of blocks present' %
                            line)
                else:
                    blocks = [
                        (int(chromStart), int(chromEnd) - int(chromStart),)]

                # Report very big annotations?
                if self.debug:
                    if max((e - s for s, e in blocks)) > 10000:
                        print('')
                        print(line, chromStart, chromEnd, blocks)
                for start, length in blocks:
                    self.addFeature(
                        chrom,
                        start,
                        start + length,
                        name,
                        strand=strand,
                        data=None)
            print("Loaded %s features, now sorting" %
                  sum([len(self.features[c]) for c in self.features]))
            self.sort()
            print("done sorting")

    def getReferenceList(self):
        return(list(self.startCoordinates.keys()))

    def addFeature(self, chromosome, start, end, name, strand=None, data=None):
        if strand is not None and strand not in ['+', '-']:
            raise ValueError('Invalid strand specified: %s' % strand)
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('Start and end coordinates should be integers')

        if self.debug:
            self.debugMsg(
                "Adding feature: chromosome:%s, start:%s, end:%s, name:%s, strand:%s, data:%s" %
                (chromosome, start, end, name, strand, data))

        if chromosome not in self.features:
            self.features[chromosome] = list()

        self.features[chromosome].append((start, end, name, strand, data))
        self.sorted = False

    def getCentroids(self):
        centroids = {}
        for chromosome in self.startCoordinates:
            centroids[chromosome] = {}

            for start, end, name, strand, data in self.features[chromosome]:
                centroids[chromosome][name] = (end - start) / 2 + start
        return(centroids)

    def findFeaturesBetween(
            self,
            chromosome,
            sampleStart,
            sampleEnd,
            strand=None):

        if chromosome not in self.startCoordinates:
            if self.debug:
                self.debugMsg(
                    "Chromosome %s is not present in the annotations" %
                    chromosome)
            return([])
        startIndex = max(
            0,
            np.searchsorted(
                self.startCoordinates[chromosome],
                sampleStart,
                'left') - 1)
        startIndex = min(
            startIndex, max(
                0, np.searchsorted(
                    self.endCoordinates[chromosome], sampleEnd, 'left')))
        hits = set()
        x = True
        while(x and startIndex < len(self.features[chromosome])):
            d = self.features[chromosome][startIndex]
            (hitStart, hitEnd, name, hitStrand, data) = d
            if hitStart > sampleEnd:
                x = False
            else:
                # Does the feature overlap the sampling region
                if(max(sampleStart, hitStart) <= min(sampleEnd, hitEnd)):
                    #print(sampleStart,sampleEnd, hitStart,hitEnd, name)
                    # Check the strand
                    if strand is None or (strand == hitStrand):
                        hits.add(d)
            startIndex += 1

        hits.update(set(self.findFeaturesAt(chromosome, sampleStart, strand)))
        hits.update(set(self.findFeaturesAt(chromosome, sampleEnd, strand)))

        return(list(hits))

    def sort(self):
        """ Build coordinate sorted datastructure to perform fast lookups."""
        self.endIndexes = {}
        self.endIndexLookup = {}
        self.fastIndex = {}
        self.maxFeatureSizes = {}
        self.sorted = True
        for chromosome in self.features.keys():
            # Sort in place and return new indices
            self.features[chromosome].sort()
            self.startCoordinates[chromosome] = np.fromiter(
                (tup[0] for tup in self.features[chromosome]), dtype=np.int64)
            self.endCoordinates[chromosome] = np.fromiter(
                (tup[1] for tup in self.features[chromosome]), dtype=np.uint64)
            self.endIndexes[chromosome] = np.argsort(
                self.endCoordinates[chromosome])
            self.endCoordinates[chromosome] = self.endCoordinates[chromosome][self.endIndexes[chromosome]]
            self.endIndexLookup[chromosome] = {
                inR: orig for orig, inR in enumerate(
                    self.endIndexes[chromosome])}

            ####################### perform magic indexing ####################
            maxLengthFeature = np.max([tup[1] - tup[0]
                                       for tup in self.features[chromosome]])
            self.maxFeatureSizes[chromosome] = maxLengthFeature

            lowestStarts = np.fromiter(
                (min(
                    (f[0] for f in self.findFeaturesAt(
                        chromosome,
                        feature[0],
                        optim='nb'))) for feature in self.features[chromosome]),
                dtype=np.int64)

            self.fastIndex[chromosome] = np.searchsorted(
                self.startCoordinates[chromosome], lowestStarts, 'left')
            ########################

        # find the longest feature

    """Return a feature left of the lookupCoordinate"""

    def findNearestLeftFeature(
            self,
            chromosome,
            lookupCoordinate,
            strand=None):
        if chromosome not in self.features:
            return([])
        if not self.sorted:
            self.sort()
        """ Find closest feature left of the supplied coordinate """
        index = np.clip(
            np.searchsorted(
                self.endCoordinates[chromosome], lookupCoordinate, side='left'), 0, len(
                self.endCoordinates))
        self.debugMsg(
            "Looking up %s, Initial index is %s (start: %s, end %s)" %
            (lookupCoordinate,
             index,
             self.features[chromosome][index][0] if index < len(
                 self.features[chromosome]) else 'No feature hit',
                self.features[chromosome][index][1] if index < len(
                 self.features[chromosome]) else 'No feature hit'))
        #index = np.clip(index, 1, len(self.endCoordinates)-1)
        # Check if the index is zero, and this feature is actually more right:
        if index > (len(self.features[chromosome]) - 1):
            self.debugMsg("overflow INDEX condition, decreasing index")
            index -= 1

        if self.features[chromosome][index][0] > lookupCoordinate:
            self.debugMsg("Zero INDEX condition. Rejected feature.")
            return([])

        hitStrand = self.features[chromosome][index][3]
        # Find first feature which is on the same strand...
        while strand is not None and (hitStrand != strand) and index > 0:
            index -= 1
            if index > 0:
                hitStrand = self.features[chromosome][index][3]
        if index < 0:
            return([])
        return([self.features[chromosome][index]])

    def findNearestRightFeature(
            self,
            chromosome,
            lookupCoordinate,
            strand=None):
        if chromosome not in self.features:
            return([])
        if not self.sorted:
            self.sort()
        """ Find closest feature left of the supplied coordinate """
        index = np.searchsorted(
            self.startCoordinates[chromosome],
            lookupCoordinate + 1,
            side='left')
        self.debugMsg(
            "Looking up %s, Initial index is %s (start: %s, end %s), strand %s" %
            (lookupCoordinate,
             index,
             self.features[chromosome][index][0] if index < len(
                 self.features[chromosome]) else 'No feature hit',
                self.features[chromosome][index][1] if index < len(
                 self.features[chromosome]) else 'No feature hit',
                strand))

        # Check if the index is maxlen, and this feature is actually more right
        while not index < len(self.features[chromosome]):
            self.debugMsg("No feature on the right")
            return([])
            index -= 1

        self.debugMsg("Index is now %s" % index)
        if self.features[chromosome][index][1] < lookupCoordinate:
            self.debugMsg("Zero INDEX condition. Rejected feature.")
            return([])

        # print(self.features[chromosome][index])
        hitStrand = self.features[chromosome][index][3]
        # Find first feature which is on the same strand...
        while strand is not None and (hitStrand != strand):
            index += 1
            if not index < len(self.features[chromosome]):
                self.debugMsg("No feature on the right")
                return([])
            hitStrand = self.features[chromosome][index][3]
        if index < 0:
            return([])
        return([self.features[chromosome][index]])

    @functools.lru_cache(maxsize=512)
    def findNearestFeature(self, chromosome, lookupCoordinate, strand=None):

        s = self.findFeaturesAt(chromosome, lookupCoordinate, strand=None)
        if len(s):
            self.debugMsg(
                'Issued nearest feature search, but the coordinate lies within %s feature(s), returning those' %
                len(s))
            return(s)

        fr = self.findNearestRightFeature(
            chromosome, lookupCoordinate=lookupCoordinate, strand=strand)
        fl = self.findNearestLeftFeature(
            chromosome, lookupCoordinate=lookupCoordinate, strand=strand)
        self.debugMsg(
            'Feature R presence %s, feature L presence: %s' %
            (len(fr), len(fl)))
        if len(fr) == 0 and len(fl) == 0:
            self.debugMsg("Returning no hits")
            return([])
        elif len(fr) == 0:
            self.debugMsg("Returning Left as right is empty")
            return(fl)
        elif len(fl) == 0:
            self.debugMsg("Returning Right as left is empty")
            return(fr)
        else:
            distanceR, distanceL = fr[0][0] - \
                lookupCoordinate, lookupCoordinate - fl[0][1]
            self.debugMsg("Distances: %s and %s" % (distanceR, distanceL))
            if distanceR < distanceL:
                return([fr[0]])
            elif distanceL < distanceR:
                return([fl[0]])
            else:
                return([fl[0], fr[0]])

    @functools.lru_cache(maxsize=512)
    def findFeaturesAt(
            self,
            chromosome,
            lookupCoordinate,
            strand=None,
            optim='bdbnb'):
        return self._findFeaturesAt(
                chromosome,
                lookupCoordinate,
                strand=strand,
                optim=optim)

    def _findFeaturesAt(
            self,
            chromosome,
            lookupCoordinate,
            strand=None,
            optim='bdbnb'):
        if not self.sorted:
            self.sort()
        """Obtain the features at a give coordinate and optionally strand."""

        if chromosome not in self.startCoordinates:
            if self.debug:
                self.debugMsg(
                    "Chromosome %s is not present in the annotations" %
                    chromosome)
            return([])

        s = np.searchsorted(
            self.startCoordinates[chromosome],
            lookupCoordinate + 1,
            side='left')
        if optim == 'bdbnb':
            startRange = self.fastIndex[chromosome][s - 1]
            # We stop looking at the rightmost h
            endRange = min(s, len(self.features[chromosome]))
            if self.debug:
                self.debugMsg("index: %s" % self.fastIndex[chromosome])
                self.debugMsg("Fast index result: %s" %
                              self.fastIndex[chromosome][s - 1])
                ml = lookupCoordinate - self.maxFeatureSizes[chromosome]
                self.debugMsg(
                    "Vanilla index result: %s" %
                    np.searchsorted(
                        self.startCoordinates[chromosome],
                        ml,
                        side='left'))
                self.debugMsg("Start looking from %s" % startRange)
                self.debugMsg("To %s" % endRange)

            return([self.features[chromosome][i] for i in range(startRange, endRange) if (self.features[chromosome][i][1] >= lookupCoordinate and (strand is None or self.features[chromosome][i][3] == strand))])
        elif optim == 'nb':
            # Be smarter: take only segments where the end coordinate is bigger
            # than the lookupCoordinate
            # We start looking from the most left index possible given our
            # knowledge of the longest feature:
            ml = lookupCoordinate - self.maxFeatureSizes[chromosome]
            startRange = np.searchsorted(
                self.startCoordinates[chromosome],
                ml,
                side='left')  # Find where the left most index lies
            # For more optimalisation we want to skip the searchsorted and
            # prefetch the lookup array?

            # We stop looking at the leftmost h
            endRange = min(s, len(self.features[chromosome]))
            if self.debug:
                self.debugMsg("Start looking from %s" % startRange)
                self.debugMsg("To %s" % endRange)
                #self.debugMsg("start is %s"%self.startIndexLookup[chromosome][k])
            candidates = set(self.features[chromosome][i] for i in range(
                startRange, endRange) if self.features[chromosome][i][1] >= lookupCoordinate)

        elif optim == 'optim':
            s = np.searchsorted(
                self.startCoordinates[chromosome],
                lookupCoordinate + 1,
                side='left')
            # Be smarter: take only segments where the end coordinate is bigger
            # than the lookupCoordinate
            candidates = set(
                self.features[chromosome][i] for i in range(
                    0, min(
                        s, len(
                            self.features[chromosome]))) if self.features[chromosome][i][1] >= lookupCoordinate)
        else:
            e = np.searchsorted(
                self.endCoordinates[chromosome],
                lookupCoordinate - 1,
                side='right')
            leftSet = set(
                self.features[chromosome][i] for i in range(
                    0, min(
                        s, len(
                            self.features[chromosome]))))
            rightSet = set(self.features[chromosome][self.endIndexLookup[chromosome][i]] for i in range(
                min(e, len(self.features[chromosome])), len(self.features[chromosome])))
            candidates = leftSet.intersection(rightSet)
        return([candidate for candidate in candidates if (strand is None or candidate[3] == strand)])
        """
        if self.debug:
            self.debugMsg("Looking up %s %s %s Hit %s %s" % (chromosome, lookupCoordinate, strand, s, e))
            self.debugMsg(self.startCoordinates[chromosome])
        hits = []
        for i in range(s, self.endIndexLookup[chromosome][e] ): #min(s+1, len(self.features[chromosome]))
            if self.features[chromosome][i][0]<=lookupCoordinate and self.features[chromosome][i][1]>=lookupCoordinate:
                if self.debug:
                    self.debugMsg("  Strandless Hit feature %s %s %s [%s %s]" % (self.features[chromosome][i][0], self.features[chromosome][i][1], self.features[chromosome][i][2], self.features[chromosome][i][3],  self.features[chromosome][i][4]))
                if strand is None or strand==self.features[chromosome][i][3]:
                    hits.append(self.features[chromosome][i])
                elif self.debug:
                    self.debugMsg("   Strand miss %s for %s %s" % (strand, i, self.features[chromosome][i][3] ))
            elif self.debug:
                self.debugMsg("  Missed feature %s %s %s [%s %s]" % (self.features[chromosome][i][0], self.features[chromosome][i][1], self.features[chromosome][i][2], self.features[chromosome][i][3],  self.features[chromosome][i][4]))
                self.debugMsg("    reason:"  "%s <= coord" % self.features[chromosome][i][0] if  self.features[chromosome][i][0]<=lookupCoordinate else "%s > coord" % self.features[chromosome][i][1] )
        return(hits)
        """

    def findFeaturesAtPysamAlign(self, pysamRead, strand=None, method=1):
        """Obtain all features mapping the pysam aligned segment.
        method 0: Query EVERY base
        method 1: Query every subsequent block of reads (pysam aligned segment .get_blocks)
        """
        if pysamRead.reference_name not in self.startCoordinates:
            if self.debug:
                self.debugMsg(
                    "Chromosome %s is not present in the annotations" %
                    pysamRead.reference_name)
            return([])
        if method == 0:
            hits = set()
            for queryPos, referencePos in pysamRead.get_aligned_pairs(
                    matches_only=True, with_seq=False):
                hits.update(set(self.findFeaturesAt(
                    pysamRead.reference_name, referencePos, strand=strand)))
            return(hits)
        else:
            return(set(itertools.chain.from_iterable([
                self.findFeaturesBetween(
                    pysamRead.reference_name,
                    lookupCoordinateStart,
                    lookupCoordinateEnd,
                    strand=strand)
                for lookupCoordinateStart, lookupCoordinateEnd
                in pysamRead.get_blocks()])))

    def findFeaturesBetweenBRK(
            self,
            chromosome,
            lookupCoordinateStart,
            lookupCoordinateEnd,
            strand=None):
        """Obtain all features between Start and end coordinate."""
        if chromosome not in self.startCoordinates:
            if self.debug:
                self.debugMsg(
                    "Chromosome %s is not present in the annotations" %
                    chromosome)
            return([])

        """Find all features which are present at both the supplied start and end coordinate"""
        startHits = self.findFeaturesAt(
            chromosome, lookupCoordinateStart, strand)
        endHits = self.findFeaturesAt(chromosome, lookupCoordinateEnd, strand)
        if self.debug:
            self.debugMsg("Hits at %s %s %s" %
                          (chromosome, lookupCoordinateEnd, strand))
            self.debugMsg("  %s" % startHits)
            self.debugMsg("  %s" % endHits)
        return(list(set(startHits).intersection(set(endHits))))

    def loadSNPSFromVcf(self, vcfFilePath, locations=None):
        for rec in pysam.VariantFile(vcfFilePath):
            if locations is None or (rec.chrom, rec.pos) in locations:
                for sample in rec.samples:
                    # get genotypes:
                    for allele in rec.samples[sample].alleles:
                        try:
                            self.addVariant(
                                chromosome=rec.chrom,
                                start=rec.pos,
                                value=allele,
                                name='SNP',
                                variantType='SNP',
                                end=None)
                        except BaseException:
                            pass

    def addVariant(
            self,
            chromosome,
            start,
            value=None,
            name='SNP',
            variantType='SNP',
            end=None):
        end = end if end is not None else start + 1
        if value not in ['A', 'T', 'C', 'G']:
            raise ValueError('%s is not a base' % value)
        self.addFeature(chromosome, start, end, name=name, data=('SNP', value))


def massIdConvert(
        baseIds,
        pathToIdMapping='/media/sf_data/references/human/HUMAN_9606_idmapping_selected.tab.gz',
        targetCol=1):
    """Convert GENE identifiers into another format.
    Get a conversion table from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
    """
    converted = {}
    baseIds = set(baseIds)
    h = gzip.open(pathToIdMapping)
    for l in h:
        line = l.decode('utf8')
        for identifier in baseIds:
            if identifier in line:
                convertedTo = line.split()[targetCol]
                if len(convertedTo):
                    if identifier not in converted:
                        converted[identifier] = []
                    converted[identifier].append(convertedTo)

    return(converted)


class FeatureAnnotatedObject():

    def __init__(self, features, stranded, capture_locations, auto_set_intron_exon_features ):

        self.features = features
        self.hits = collections.defaultdict(set)  # feature -> hit_bases
        self.stranded = stranded
        self.is_annotated = False
        self.capture_locations = capture_locations
        if capture_locations:
            self.feature_locations = {} #feature->locations (chrom,start,end, strand)

        self.junctions = set()
        self.genes = set()
        self.introns = set()
        self.exons = set()
        self.exon_hit_gene_names = set()  # readable names
        self.is_spliced = None

        if auto_set_intron_exon_features:
            self.set_intron_exon_features()


    def set_spliced(self, is_spliced):
        """ Set wether the transcript is spliced, False has priority over True """
        if self.is_spliced and not is_spliced:
            # has already been set
            self.is_spliced = False
        else:
            self.is_spliced = is_spliced

    def get_hit_df(self):
        """Obtain dataframe with hits
        Returns:
            pd.DataFrame
        """
        if not self.is_annotated:
            self.set_intron_exon_features()

        d = {}
        tabulated_hits = []
        for hit, locations in self.hits.items():
            if not isinstance(hit, tuple):
                continue
            meta = dict(list(hit))
            for location in locations:
                location_dict = {
                    'chromosome': location[0],
                    'start': location[1][0],
                    'end': location[1][1]}
                location_dict.update(meta)
                tabulated_hits.append(location_dict)

        return pd.DataFrame(tabulated_hits)

    def write_tags(self):

        if len(self.exons) > 0:
            self.set_meta('EX', ','.join(sorted([str(x) for x in self.exons])))
        else:
            self.set_meta('EX',None)

        if len(self.introns) > 0:
            self.set_meta('IN', ','.join(
                sorted([str(x) for x in self.introns])))
        else:
            self.set_meta('IN',None)

        if len(self.genes) > 0:
            self.set_meta('GN', ','.join(sorted([str(x) for x in self.genes])))
        else:
            self.set_meta('GN',None)

        if len(self.junctions) > 0:
            self.set_meta('JN', ','.join(
                sorted([str(x) for x in self.junctions])))
            # Is transcriptome
            self.set_meta('IT', 1)
        elif len(self.genes) > 0:
            # Maps to gene but not junction
            self.set_meta('IT', 0.5)
            self.set_meta('JN',None)
        else:
            # Doesn't map to gene
            self.set_meta('IT', 0)
            self.set_meta('JN', None)

        if self.is_spliced is True:
            self.set_meta('SP', True)
        elif self.is_spliced is False:
            self.set_meta('SP', False)
        if len(self.exon_hit_gene_names):
            self.set_meta('gn', ';'.join(list(self.exon_hit_gene_names)))
        else:
            self.set_meta('gn',None)

    def set_intron_exon_features(self):
        if not self.is_annotated:
            self.annotate()

        # Collect all hits:
        hits = self.hits.keys()

        # (gene, transcript) -> set( exon_id  .. )
        exon_hits = collections.defaultdict(set)
        intron_hits = collections.defaultdict(set)

        for hit, locations in self.hits.items():
            if not isinstance(hit, tuple):
                continue

            meta = dict(list(hit))
            if 'gene_id' not in meta:
                continue
            if meta.get('type') == 'exon':
                if 'transcript_id' not in meta:
                    continue
                self.genes.add(meta['gene_id'])
                self.exons.add(meta['exon_id'])
                if 'transcript_id' not in meta:
                    raise ValueError(
                        "Please use an Intron GTF file generated using -id 'transcript_id'")
                exon_hits[(meta['gene_id'], meta['transcript_id'])].add(
                    meta['exon_id'])
                if 'gene_name' in meta:
                    self.exon_hit_gene_names.add(meta['gene_name'])
            elif meta.get('type') == 'intron':
                self.genes.add(meta['gene_id'])
                self.introns.add(meta['gene_id'])

        # Find junctions and add all annotations to annotation sets
        debug = []

        for (gene, transcript), exons_overlapping in exon_hits.items():
            # If two exons are detected from the same gene we detected a
            # junction:
            if len(exons_overlapping) > 1:
                self.junctions.add(gene)

                # We found two exons and an intron:
                if gene in self.introns:
                    self.set_spliced(False)
                else:
                    self.set_spliced(True)

            debug.append(
                f'{gene}_{transcript}:' +
                ','.join(
                    list(exons_overlapping)))

        # Write exon dictionary:
        self.set_meta('DB', ';'.join(debug))


if __name__ == "__main__":
    """The following are all test functions for the annotation class"""

    from colorama import Fore  # ,Back, Style
    from colorama import Back
    from colorama import Style
    from colorama import init
    init(autoreset=True)

    def formatColor(string):
        return(string.replace("[GREEN]", Fore.GREEN).replace("[RED]", Fore.RED).replace("[DIM]", Style.DIM).replace("[RESET]", Style.RESET_ALL).replace("[BRIGHT]", Style.BRIGHT).replace("[NORMAL]", Style.NORMAL))

    def printFormatted(string):
        print(formatColor(str(string)))

    def printFormattedDim(string):
        print(formatColor("   [DIM]%s" % str(string)))
    """Self tests:"""

    # addFeature( chromosome, start, end, name, strand=None, data=None):
    # Build the reference:
    print(
        """
    .......(1)---------------------->
    ..............<----------------(2)
    chr1   100    110             200


    .......(3)---------------------->
    .......(4)----->
    .......(5)<-------------
    chr2   100    110     150     200
    .
    """
    )

    def expect(result, desired, presenceTestOnly=False):
        if presenceTestOnly:
            printFormatted(
                "Expecting %s, %s" %
                (desired,
                 ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" %
                  result) if any(
                     r[2] in desired for r in result) else '[RED]FAIL %s\n' %
                    result))
            return
        if desired is None:
            printFormatted(
                "Expecting %s, %s" %
                (desired,
                 ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" %
                  result) if len(result) == 0 else '[RED]FAIL %s' %
                    result))
        elif isinstance(desired, list):
            printFormatted(
                "Expecting %s, %s" %
                (desired, ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" %
                           result) if len(result) == len(desired) and all(
                    r[2] in desired for r in result) else '[RED]FAIL %s\n' %
                    result))
        else:
            printFormatted(
                "Expecting %s, %s" %
                (desired,
                 ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" %
                  result) if len(result) == 1 and result[0][2] == desired else '[RED]FAIL %s\n' %
                    result))
    f = FeatureContainer()
    f.debug = True
    f.debugMsg = printFormattedDim
    f.addFeature('chrY', 1, 3, 'A', '+', '')
    f.addFeature('chrY', 5, 8, 'B', '+', '')
    f.sort()
    expect(f.findFeaturesAt('chrY', 0, '+'), None)
    expect(f.findFeaturesAt('chrY', 1, '+'), 'A')
    expect(f.findFeaturesAt('chrY', 2, '+'), 'A')
    expect(f.findFeaturesAt('chrY', 3, '+'), 'A')
    expect(f.findFeaturesAt('chrY', 4, '+'), None)
    expect(f.findFeaturesAt('chrY', 5, '+'), 'B')
    expect(f.findFeaturesAt('chrY', 6, '+'), 'B')
    expect(f.findFeaturesAt('chrY', 8, '+'), 'B')

    f = FeatureContainer()
    f.debug = True
    f.debugMsg = printFormattedDim
    f.addFeature('chrX', 10, 1000, 'parentB', '+', '')
    f.addFeature('chrX', 500, 900, 'nestedB', '+', '')
    f.addFeature('chrX', 10000, 12000, 'C', '+', '')
    f.addFeature('chrX', 100000, 120000, 'D', '+', '')
    f.sort()
    print(f.findFeaturesAt('chrX', 550, '+'))
    print(f.findFeaturesAt('chrX', 10000, '+'))
    print(f.findFeaturesAt('chrX', 100000000, '+'))
    expect(f.findFeaturesAt('chrX', 9, '+'), None)
    expect(f.findFeaturesAt('chrX', 12001, '+'), None)
    expect(f.findFeaturesAt('chrX', 12000, '+'), 'C')
    expect(f.findFeaturesAt('chrX', 120000, '+'), 'D')
    expect(f.findFeaturesAt('chrX', 10, '+'), 'parentB')

    f = FeatureContainer()
    f.debug = True
    f.debugMsg = printFormattedDim

    f.addFeature(
        'chr1',
        100,
        200,
        '1',
        '+',
        'A forward feature from 100 to 200 chr1')
    f.addFeature(
        'chr1',
        110,
        200,
        '2',
        '-',
        'A reverse feature from 110 to 200 chr1')
    f.addFeature(
        'chr2',
        100,
        200,
        '3',
        '+',
        'A forward feature from 100 to 200 chr2')
    f.addFeature(
        'chr2',
        100,
        110,
        '4',
        '+',
        'A forward feature from 100 to 110 chr2')
    f.addFeature(
        'chr2',
        100,
        150,
        '5',
        '-',
        'A reverse feature from 100 to 150 chr2')

    f.addFeature('chr3', 100, 150, '6', '-', 'feature 6')
    f.addFeature('chr3', 200, 250, '7', '-', 'feature 7')
    f.addFeature('chr3', 200, 450, '8', '-', 'feature 8')
    f.addFeature('chr3', 10, 15, '9', '-', 'feature 9')

    printFormatted("[BRIGHT]Test for reference presence:")
    result = f.getReferenceList()
    desired = ['chr1', 'chr2', 'chr3']
    printFormatted(
        "Expecting %s, %s" %
        (desired,
         ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" %
          result) if len(result) == len(desired) and all(
             r in desired for r in result) else '[RED]FAIL %s\n' %
            result))

    printFormatted("[BRIGHT]Test on leftmost start:")
    expect(f.findFeaturesAt('chr1', 100, '+'), '1')

    printFormatted("[BRIGHT]Test on rightmost end:")
    expect(f.findFeaturesAt('chr1', 200, '+'), '1')

    printFormatted("[BRIGHT]Test on random location within feature:")
    expect(f.findFeaturesAt('chr1', 120, '+'), '1')

    printFormatted("[BRIGHT]Test on random location within feature:")
    expect(f.findFeaturesAt('chr2', 120, '+'), '3')

    printFormatted("[BRIGHT]Test on limit location of feature:")
    expect(f.findFeaturesAt('chr2', 200, '+'), '3')

    printFormatted(
        "[BRIGHT]Test on non matching location (match available on other side, and one base left of coord):")
    expect(f.findFeaturesAt('chr2', 151, '-'), None)

    printFormatted(
        "[BRIGHT]Tests on double matching locations without strand spec")
    expect(f.findFeaturesAt('chr1', 120), ['1', '2'])
    expect(f.findFeaturesAt('chr2', 105), ['3', '4', '5'])

    printFormatted("[BRIGHT] ==== Range tests... ====")
    printFormatted(
        "[BRIGHT]Test for matching start and end coordinates overlapping one feature")
    expect(f.findFeaturesBetween('chr2', 102, 200, '+'), ['3', '4'])

    printFormatted("[BRIGHT]Test for matching all features on chromosome")
    expect(f.findFeaturesBetween('chr2', 0, 20000, None), ['3', '4', '5'])

    printFormatted(
        "[BRIGHT]Test for matching all but one features on chromosome")
    expect(f.findFeaturesBetween('chr3', 151, 20000, None), ['8', '7'])

    printFormatted(
        "[BRIGHT]Test for matching all but one features on chromosome")
    expect(f.findFeaturesBetween('chr3', 0, 195, None), ['6', '9'])

    printFormatted("[BRIGHT]Test for range finding non-existent feature")
    expect(f.findFeaturesBetween('chr2', 2001, 2000, '+'), None)

    printFormatted(
        "[BRIGHT]Test for finding non-existent feature LEFT NEXT to the point")
    expect(f.findNearestLeftFeature('chr2', 50, '+'), None)

    printFormatted(
        "[BRIGHT]Test for finding existent feature LEFT NEXT to the point")
    expect(f.findNearestLeftFeature('chr2', 250, None), '3')

    printFormatted(
        "[BRIGHT]Test for finding non-existent feature RIGHT NEXT to the point")
    expect(f.findNearestRightFeature('chr2', 250, '+'), None)

    printFormatted(
        "[BRIGHT]Test for finding existent feature RIGHT NEXT to the point")
    expect(f.findNearestRightFeature('chr2', 0, '-'), '5')

    printFormatted("[BRIGHT]Test for finding closest feature")
    print(f.findNearestFeature('chr1', 0, None))
    expect(f.findNearestFeature('chr1', 0, None), '1')

    # Sharing related stuff
    """from multiprocessing.managers import BaseManager
    import multiprocessing
    # Share the genome annotations over multiple processes:
    class bdbsSharedObjectManager(BaseManager): pass
    def Manager():
        m = bdbsSharedObjectManager()
        m.start()
        return m
    # initialisation #
    bdbsSharedObjectManager.register('FeatureContainer', FeatureContainer)
    sharedDataManager = Manager()
    f = sharedDataManager.FeatureContainer()


    def findFeaturesAt(args):
        featureContainer, chromosome, position, strand = args
        featureContainer.findFeaturesAt(chromosome, position, strand )

    print('Running with pool')
    pool = multiprocessing.Pool(8)
    print("With index:")
    bar = progressbar.ProgressBar(max_value=N)
    for q,result in enumerate(pool.imap( findFeaturesAt, (  (f, 'chr1', q, '+') for q in range(N) ),1000)):
         bar.update(q)
        #expect( f.findFeaturesAt('chr1',q,'+'), '1', True)
    bar.finish()
    exit()
    """
    #######################################
    import random
    import progressbar
    print('Creating random features')

    N = 1000000
    f.debug = False
    printFormatted(
        "[BRIGHT]Test with %s random features added to the chromosome" %
        N)
    expect(f.findFeaturesAt('chr1', 100, '+'), '1', True)
    for i in range(0, N):
        s = random.randint(0, 100_000_000)
        f.addFeature('chr1', s, s +
                     random.randint(1, 1000), 'art_%s' %
                     i, ['-', '+'][random.randint(0, 1)], 'A random feature')
        #f.findNearestFeature('chr1', random.randint(0,1000), '+')
    print('Constructing index...')
    f.sort()
    #######################################

    print("With index:")
    bar = progressbar.ProgressBar(max_value=N)
    for q in range(N):
        f.findFeaturesAt('chr1', q, '+')
        bar.update(q)
        #expect( f.findFeaturesAt('chr1',q,'+'), '1', True)
    bar.finish()

    print("Without index:")
    bar = progressbar.ProgressBar(max_value=N)
    for q in range(N):
        f.findFeaturesAt('chr1', q, '+', optim='nb')
        bar.update(q)
        #expect( f.findFeaturesAt('chr1',q,'+'), '1', True)
    bar.finish()
