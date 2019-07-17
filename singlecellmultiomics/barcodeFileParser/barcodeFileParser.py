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


#http://codereview.stackexchange.com/questions/88912/create-a-list-of-all-strings-within-hamming-distance-of-a-reference-string-with
def hamming_circle(s, n, alphabet):
    for positions in itertools.combinations(range(len(s)), n):
        for replacements in itertools.product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)


class BarcodeParser():

    def __init__(self, barcodeDirectory='barcodes', hammingDistanceExpansion=0,spaceFill=False ):


        barcodeDirectory= os.path.join( os.path.dirname(os.path.realpath(__file__)), barcodeDirectory)
        barcode_files = list(glob.glob(f'{barcodeDirectory}/*'))

        self.spaceFill = spaceFill
        self.hammingDistanceExpansion = hammingDistanceExpansion
        self.barcodes = collections.defaultdict(dict) # alias -> barcode -> index
        self.extendedBarcodes = collections.defaultdict(dict) # alias -> barcode -> (index, hammingDistance)

        for barcodeFile in barcode_files:
            barcodeFileAlias  = os.path.splitext(os.path.basename(barcodeFile))[0]
            logging.info(f"Parsing {barcodeFile}, alias {barcodeFileAlias}")

            # Decide the file type (index first or name first)
            indexNotFirst = False
            with open(barcodeFile) as f:
                for i, line in enumerate(f):
                    parts = line.strip().split()
                    if len(parts)==1 and ' ' in line:
                        parts = line.strip().split(' ')
                    if len(parts)==1:
                        pass
                    elif len(parts)==2:
                        indexFirst = not all((c in 'ATCGNX') for c in parts[0])
                        if not indexFirst:
                            indexNotFirst=True
                        #print(parts[1],indexFirst)

            with open(barcodeFile) as f:
                for i, line in enumerate(f):
                    parts = line.strip().split()
                    if len(parts)==1 and ' ' in line:
                        parts = line.strip().split(' ')
                    if len(parts)==1:
                        self.addBarcode( barcodeFileAlias, barcode=parts[0], index=i)
                        #self.barcodes[barcodeFileAlias][parts[0]] = i
                        logging.info(f"\t{parts[0]}:{i} (No index specified in file)")
                    elif len(parts)==2:
                        if indexNotFirst:
                            barcode, index = parts
                        else:
                            index,barcode  = parts
                        #self.barcodes[barcodeFileAlias][barcode] = index
                        self.addBarcode( barcodeFileAlias, barcode=barcode, index=index)
                        logging.info(f"\t{barcode}:{index} (index was specified in file, {'index' if indexFirst else 'barcode'} on first column)")
                    else:
                        e = f'The barcode file {barcodeFile} contains more than two columns. Failed to parse!'
                        logging.error(e)
                        raise ValueError(e)

        if hammingDistanceExpansion>0:
            for alias in list(self.barcodes.keys()):
                #print(alias)
                self.expand(hammingDistanceExpansion, alias=alias)

    def getTargetCount(self, barcodeFileAlias):
        return( len(self.barcodes[barcodeFileAlias] ), len(self.extendedBarcodes[barcodeFileAlias]) )


    def expand(self,hammingDistanceExpansion, alias, reportCollisions=True, spaceFill=None):

        barcodes = self.barcodes[alias]
        hammingSpace = collections.defaultdict(list) # hammingBarcode -> ( ( distance,origin) )
        for barcode in barcodes:

            for hammingDistance in range(0, hammingDistanceExpansion+1):
                for hammingInstance  in hamming_circle(barcode,hammingDistance,'ACTGN' ):
                    hammingSpace[hammingInstance].append((hammingDistance,barcode))
        # Resolve all
        for hammingBarcode in hammingSpace:

            # Check if there is a closest origin:
            sortedDistances = sorted(hammingSpace[hammingBarcode])
            if len(sortedDistances)>1 and ( sortedDistances[0][0] ==  sortedDistances[1][0]):
                # We cannot resolve this, two or more origins are at the same distance:
                #print('Cannot resolve %s' % hammingBarcode)
                continue

            hammingDistance, origin = sortedDistances[0]

            self.addBarcode( alias, barcode=hammingBarcode, index=self.barcodes[alias][origin], hammingDistance=hammingDistance, originBarcode=origin)


    def expand_old(self,hammingDistanceExpansion, alias, reportCollisions=True, spaceFill=None ): # Space fill fills all hamming instances, even if they are not resolvable

        if spaceFill is None:
            spaceFill=self.spaceFill
        #print("Expanding Hamming distance %s" % alias)
        hammingMatrix = {}
        collisions = collections.Counter()
        collisionsPerBarcode = {}
        barcodes = self.barcodes[alias]

        # Iterate all barcodes
        for barcode in barcodes:
            k = barcodes[barcode]
            # Perform iteration per hamming distance
            for hammingDistance in range(0, hammingDistanceExpansion+1):
                for hammingInstance  in hamming_circle(barcode,hammingDistance,'ACTGN' ):
                    if spaceFill:
                        self.addBarcode( alias, barcode=hammingInstance, index=k, hammingDistance=hammingDistance, originBarcode=barcode)
                        continue
                    # The hamming instance is a hammingDistance mutation of the current barcode
                    if not hammingInstance in hammingMatrix:
                        # The hamming Matrix contains [ mutatedBarcode ] = [ targetIndex, distance, [(originBarcode, targetIndex), ... ], originBarcode]
                        hammingMatrix[hammingInstance] = [k, hammingDistance, [(k, barcode)], barcode]
                    else:
                        if k != hammingMatrix[hammingInstance][0]: # The means another barcode was already assigned to this mutated version
                            collisions[hammingInstance]+=1
                            hammingMatrix[hammingInstance][2].append((k,barcode))
                        else: # There is no collision
                            raise ValueError('This should never happen')
                            continue
                        # Getting here means there is a collision
                        if hammingMatrix[hammingInstance][1] == hammingDistance: #Collision
                            hammingMatrix[hammingInstance][0] = None # set collision
                        # Maybe there is a collision, but one barcode is further away than the other:
                        elif hammingMatrix[hammingInstance][1] > hammingDistance:

                            hammingMatrix[hammingInstance][0] = k
                            hammingMatrix[hammingInstance][1] = hammingDistance
                            hammingMatrix[hammingInstance][3] = barcode

        if sum(collisions.values())>0 and reportCollisions:
            print("%s%s collisions for %s in Hamming space: %s" % (Fore.RED, sum(collisions.values()), alias ,Style.RESET_ALL ))
            #for hammingDistance in range(0, args.hd+1):
            #    print("%s %s collisions" % (hammingDistance,collisions[hammingDistance]))
            showCollisions = 20
            shown = 0
            totalCollisions = 0
            for idx,hammingInstance in enumerate(hammingMatrix):
                assignedSample, hammingDistance, collidingSamples, origin = hammingMatrix[hammingInstance]
                if len(collidingSamples)>1:
                    if shown<showCollisions:
                        print("%s\t%sat%s %s %sdistance%s %s" % (",".join(["%s%s[%s]%s" % (Fore.GREEN if x[0] is assignedSample else Fore.RED ,x[0], x[1], Style.RESET_ALL) for x in collidingSamples]), Style.DIM, Style.RESET_ALL, hammingInstance, Style.DIM, Style.RESET_ALL,hammingDistance))
                        shown+=1
                    totalCollisions+=1
            if totalCollisions>showCollisions:
                print(('%s more ...\n' % (totalCollisions-shown)))

        if not spaceFill:
            mapping = {} # @todo: we don't need this variable anymore
            for sequence in hammingMatrix:
                if hammingMatrix[sequence][0]!=None:
                    mapping[sequence] = hammingMatrix[sequence][0]
                    self.addBarcode( alias, barcode=sequence, index=hammingMatrix[sequence][0], hammingDistance=hammingMatrix[sequence][1], originBarcode=hammingMatrix[sequence][3])

        #Show a couple of barcodes
        if False:
            print(('\n%s hamming extended barcodes will be used for demultiplexing strategy %s' % (len(mapping), alias)))
            if len(mapping):
                showBarcodes = 7
                for bcId in list(mapping.keys())[:showBarcodes]:
                    print(('%s%s%s → %s%s'  % (Fore.GREEN, bcId, Fore.WHITE, mapping[bcId], Style.RESET_ALL )))
                if len(mapping)>showBarcodes:
                    print(('%s more ...\n' % (len(mapping)-showBarcodes)))


    def addBarcode(self, barcodeFileAlias, barcode, index, hammingDistance=0, originBarcode=None):
        if hammingDistance==0:
            self.barcodes[barcodeFileAlias][barcode] = index
        else:
            if originBarcode is None:
                raise ValueError()
            self.extendedBarcodes[barcodeFileAlias][barcode] = (index,originBarcode, hammingDistance)

    # get index and hamming distance to barcode  returns none if not Available
    def getIndexCorrectedBarcodeAndHammingDistance(self, barcode, alias):
        if barcode in self.barcodes[alias]:
            return (self.barcodes[alias][barcode], barcode, 0)
        if barcode in self.extendedBarcodes[alias]:
            return self.extendedBarcodes[alias][barcode]
        return (None,None,None)


    def list(self, showBarcodes=5): #showBarcodes=None show all
        for barcodeAlias,mapping in self.barcodes.items():
            print( f'{len(mapping)} barcodes{Style.DIM} obtained from {Style.RESET_ALL}{barcodeAlias}')
            if len(mapping):
                for bcId in list(mapping.keys())[:showBarcodes]:
                    try:
                        print(('%s%s%s%s%s%s'  % (Fore.GREEN, bcId, Fore.WHITE,  '→',  mapping[bcId], Style.RESET_ALL )))
                    except Exception as e:
                        print(('%s%s%s%s%s%s'  % (Fore.GREEN, bcId, Fore.WHITE,  '->',  mapping[bcId], Style.RESET_ALL )))
                if showBarcodes is not None and len(mapping)>showBarcodes:
                    print(f'{Style.DIM} %s more ...\n{Style.RESET_ALL}' % (len(mapping)-showBarcodes))


    def getBarcodeMapping(self):
        return self.barcodes
