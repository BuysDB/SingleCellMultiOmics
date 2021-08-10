#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import logging
from colorama import Fore
from colorama import Back
from colorama import Style
import os
import collections
import itertools
import gzip

#logging.basicConfig(level=logging.DEBUG)
# http://codereview.stackexchange.com/questions/88912/create-a-list-of-all-strings-within-hamming-distance-of-a-reference-string-with
def hamming_circle(s, n, alphabet):
    for positions in itertools.combinations(range(len(s)), n):
        for replacements in itertools.product(
                range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)


class BarcodeParser():

    def path_to_barcode_alias(self, path):
        barcodeFileAlias = os.path.splitext(
            os.path.basename(path))[0].replace('.gz','').replace('.bc','')
        return barcodeFileAlias

    def parse_pending_barcode_file_of_alias(self, alias):
        if not alias in self.pending_files:
            raise ValueError(f"trying to load {alias}, but barcode file not found")

        logging.info(f"Loading promised barcode file {alias}")
        self.parse_barcode_file( self.pending_files[alias] )
        logging.info(f"Performing hamming extension for {alias}")
        self.expand(self.hammingDistanceExpansion, alias=alias)
        del self.pending_files[alias]

    def parse_barcode_file(self, barcodeFile):

        barcodeFileAlias = self.path_to_barcode_alias(barcodeFile)
        # Decide the file type (index first or name first)
        indexNotFirst = False
        with gzip.open(barcodeFile,'rt') if barcodeFile.endswith('.gz') else open(barcodeFile) as f :
            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) == 1 and ' ' in line:
                    parts = line.strip().split(' ')
                if len(parts) == 1:
                    pass
                elif len(parts) == 2:
                    indexFirst = not all((c in 'ATCGNX') for c in parts[0])
                    if not indexFirst:
                        indexNotFirst = True
                    # print(parts[1],indexFirst)

        nospec=False
        i=0
        with gzip.open(barcodeFile,'rt') if barcodeFile.endswith('.gz') else open(barcodeFile) as f :

            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) == 1 and ' ' in line:
                    parts = line.strip().split(' ')
                if len(parts) == 1:
                    self.addBarcode(
                        barcodeFileAlias, barcode=parts[0], index=i+1)
                    #self.barcodes[barcodeFileAlias][parts[0]] = i
                    if not nospec: # only show this once:
                        logging.info(
                            f"\t{parts[0]}:{i} (No index specified in file)")
                        nospec=True
                elif len(parts) == 2:
                    if indexNotFirst:
                        barcode, index = parts
                    else:
                        index, barcode = parts
                    #self.barcodes[barcodeFileAlias][barcode] = index

                    # When the index is only digits, convert to integer
                    try:
                        if int(index)==int(str(int(index))):
                            index = int(index)
                        else:
                            pass
                    except Exception as e:
                        pass
                    self.addBarcode(
                        barcodeFileAlias, barcode=barcode, index=index)
                    if not nospec: # only show this once:
                        logging.info(
                            f"\t{barcode}:{index} (index was specified in file, {'index' if indexFirst else 'barcode'} on first column)")
                        nospec=True
                else:
                    e = f'The barcode file {barcodeFile} contains more than two columns. Failed to parse!'
                    logging.error(e)
                    raise ValueError(e)
        logging.info(f'done loading {i} barcodes')


    def __init__(
            self,
            barcodeDirectory='barcodes',
            hammingDistanceExpansion=0,
            spaceFill=False,
            lazyLoad=None # these aliases wiill not be loaded until requested, '*' matches all files
            ):

        self.hammingDistanceExpansion = hammingDistanceExpansion

        barcodeDirectory = os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)),
            barcodeDirectory)
        barcode_files = list(glob.glob(f'{barcodeDirectory}/*'))

        self.spaceFill = spaceFill
        self.hammingDistanceExpansion = hammingDistanceExpansion
        self.barcodes = collections.defaultdict(
            dict)  # alias -> barcode -> index
        # alias -> barcode -> (index, hammingDistance)
        self.extendedBarcodes = collections.defaultdict(dict)
        self.pending_files = dict()
        for barcodeFile in barcode_files:
            barcodeFileAlias = self.path_to_barcode_alias(barcodeFile)
            if lazyLoad is not None and (barcodeFileAlias in lazyLoad or lazyLoad =='*'):
                logging.info(f"Lazy loading {barcodeFile}, alias {barcodeFileAlias}")
                self.pending_files[barcodeFileAlias] = barcodeFile
                continue
            logging.info(f"Parsing {barcodeFile}, alias {barcodeFileAlias}")
            self.parse_barcode_file(barcodeFile)

            if hammingDistanceExpansion > 0:
                self.expand(hammingDistanceExpansion, alias=barcodeFileAlias)

    def getTargetCount(self, barcodeFileAlias):
        return(len(self.barcodes[barcodeFileAlias]), len(self.extendedBarcodes[barcodeFileAlias]))

    def expand(
            self,
            hammingDistanceExpansion,
            alias,
            reportCollisions=True,
            spaceFill=None):

        barcodes = self.barcodes[alias]
        # hammingBarcode -> ( ( distance,origin) )
        hammingSpace = collections.defaultdict(list)
        for barcode in barcodes:

            for hammingDistance in range(0, hammingDistanceExpansion + 1):
                for hammingInstance in hamming_circle(
                        barcode, hammingDistance, 'ACTGN'):
                    hammingSpace[hammingInstance].append(
                        (hammingDistance, barcode))
        # Resolve all
        for hammingBarcode in hammingSpace:

            # Check if there is a closest origin:
            sortedDistances = sorted(hammingSpace[hammingBarcode])
            if len(sortedDistances) > 1 and (
                    sortedDistances[0][0] == sortedDistances[1][0]):
                # We cannot resolve this, two or more origins are at the same distance:
                #print('Cannot resolve %s' % hammingBarcode)
                continue

            hammingDistance, origin = sortedDistances[0]

            self.addBarcode(
                alias,
                barcode=hammingBarcode,
                index=self.barcodes[alias][origin],
                hammingDistance=hammingDistance,
                originBarcode=origin)

    def addBarcode(
            self,
            barcodeFileAlias,
            barcode,
            index,
            hammingDistance=0,
            originBarcode=None):
        if hammingDistance == 0:
            self.barcodes[barcodeFileAlias][barcode] = index
        else:
            if originBarcode is None:
                raise ValueError()
            self.extendedBarcodes[barcodeFileAlias][barcode] = (
                index, originBarcode, hammingDistance)

    # get index and hamming distance to barcode  returns none if not Available
    def getIndexCorrectedBarcodeAndHammingDistance(self, barcode, alias, try_lazy_load_pending=True):

        # Check if the alias still needs to be loaded:


        if barcode in self.barcodes[alias]:
            return (self.barcodes[alias][barcode], barcode, 0)
        if barcode in self.extendedBarcodes[alias]:
            return self.extendedBarcodes[alias][barcode]

        if alias in self.pending_files:
            if not try_lazy_load_pending:
                raise RecursionError()
            self.parse_pending_barcode_file_of_alias(alias)
            return self.getIndexCorrectedBarcodeAndHammingDistance(barcode, alias, try_lazy_load_pending=False)

        return (None, None, None)

    def list(self, showBarcodes=5):  # showBarcodes=None show all
        for barcodeAlias, mapping in self.barcodes.items():
            print(
                f'{len(mapping)} barcodes{Style.DIM} obtained from {Style.RESET_ALL}{barcodeAlias}')
            if len(mapping):
                for bcId in list(mapping.keys())[:showBarcodes]:
                    try:
                        print(('%s%s%s%s%s%s' % (Fore.GREEN, bcId,
                                                 Fore.WHITE, '→', mapping[bcId], Style.RESET_ALL)))
                    except Exception as e:
                        print(('%s%s%s%s%s%s' % (Fore.GREEN, bcId,
                                                 Fore.WHITE, '->', mapping[bcId], Style.RESET_ALL)))
                if showBarcodes is not None and len(mapping) > showBarcodes:
                    print(
                        f'{Style.DIM} %s more ...\n{Style.RESET_ALL}' %
                        (len(mapping) - showBarcodes))

    def getBarcodeMapping(self):
        return self.barcodes
