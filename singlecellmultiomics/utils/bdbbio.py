from Bio import SeqIO
from Bio.Seq import Seq
import os
import uuid
from Bio import pairwise2
import math
import subprocess
from colorama import Fore #,Back, Style
from colorama import Back
from colorama import Style
from colorama import init
init(autoreset=True)
import numpy as np
import re
from itertools import takewhile,repeat
#import localInstall

if not os.name=='nt':
    import pysam
### custom path config: (not used when not available)

#localInstall.addToPath('/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1/')
#localInstall.addToPath('/hpc/hub_oudenaarden/bin/software/bwa-0.7.10/')
#localInstall.addToPath('/media/sf_externalTools/')

### end

class AlignmentVisualiser:
    def __init__(self, initialiser):

        if type(initialiser) is pysam.AlignedSegment:
            self.alignment = initialiser

    def getPhredScoreColor(self, phredScore):

        p = ord(phredScore)-33.0
        if p>32:
            return(Fore.GREEN + Style.BRIGHT)
        if p>25:
            return(Fore.GREEN)
        if p>20:
            return(Fore.GREEN + Style.DIM)
        if p>18:
            return(Fore.YELLOW)
        return(Fore.RED)

    def visualizeAlignment(self, annotateReference=None, annotateQuery=None, intronSummary=True):

        reference = []
        query = []
        alignmentSymbols = []
        qualitieSymbols = []
        prevIntron = 0
        for index,pair in  enumerate(self.alignment.get_aligned_pairs( matches_only=False, with_seq=True) ):
            #Pair structure = (queryPosition, referencePosition, referenceBase)

            if pair[0] is None and pair[2] is None:
                prevIntron+=1
                continue
            else:
                if prevIntron>0:
                    intronString = ' N%s ' % prevIntron
                    alignmentSymbols.append(Style.DIM + intronString + Style.RESET_ALL)
                    qualitieSymbols.append(' '*len(intronString))
                    reference.append(Style.DIM +  ('-'*len(intronString)) + Style.RESET_ALL)
                    query.append( Style.DIM +  (' '*len(intronString)) + Style.RESET_ALL)
                    prevIntron = 0

            queryPosition = pair[0]
            queryNucleotide= "%s%s-%s"% (Back.RED, Fore.WHITE, Style.RESET_ALL)
            referenceNucleotide= ' '# % (Fore.RED, Style.RESET_ALL) #if queryPosition>self.alignment.query_alignment_start and queryPosition<=self.alignment.query_alignment_end else ' '
            queryPhred=' '
            if queryPosition!=None:
                queryNucleotide = self.alignment.query_sequence[queryPosition].upper()
                queryPhred = self.alignment.qual[queryPosition]
                qualitieSymbols.append(self.getPhredScoreColor(queryPhred) + queryPhred + Style.RESET_ALL)
            else:
                qualitieSymbols.append(' ')

            #print("Q:%s  S:%s" % (queryPosition, self.alignment.query_alignment_start))

            if pair[2]!=None:
                referenceNucleotide =  pair[2].upper()
            elif queryPosition is None:
                referenceNucleotide = "%s%s-%s"% (Back.RED, Fore.WHITE, Style.RESET_ALL)
            elif queryPosition>self.alignment.query_alignment_start and queryPosition<=self.alignment.query_alignment_end:

                referenceNucleotide = "%s%s-%s"% (Back.RED, Fore.WHITE, Style.RESET_ALL)

            if pair[2]==queryNucleotide:
                alignmentSymbols.append('|')
            else:
                alignmentSymbols.append(' ')

            reference.append(referenceNucleotide)
            if queryPosition!=None and annotateQuery!=None:

                if ((queryPosition - self.alignment.query_length) in  annotateQuery and  (queryPosition - self.alignment.query_length<0 )):
                    query.append(getattr(Back,annotateQuery[(queryPosition - self.alignment.query_length)]) + Fore.BLACK + queryNucleotide + Style.RESET_ALL)
                elif queryPosition in annotateQuery:
                    query.append(getattr(Back,annotateQuery[queryPosition]) + Fore.WHITE + queryNucleotide + Style.RESET_ALL)
                else:
                    #DIM soft clipped things:
                    if queryPosition<self.alignment.query_alignment_start or queryPosition>self.alignment.query_alignment_end:
                        query.append(Style.DIM + Fore.YELLOW + queryNucleotide + Style.RESET_ALL)

                    else:
                        query.append(queryNucleotide)
            else:
                query.append(queryNucleotide)
        rStart = self.alignment.reference_end if  self.alignment.is_reverse else self.alignment.reference_start
        rEnd =  self.alignment.reference_start if  self.alignment.is_reverse else self.alignment.reference_end
        refStartStr = str(rStart)
        qStartStr = str(self.alignment.query_alignment_start)
        labelSize = max(len(refStartStr), len(qStartStr))


        annotatedCigar = []
        for operation,length in self.alignment.cigartuples:
            annotatedCigar.append( '%s%s%s' % (
                [ '%s%sM%s'% (Fore.GREEN, Style.BRIGHT, Style.NORMAL),
                 '%s%sI%s' % (Fore.YELLOW,Style.BRIGHT, Style.NORMAL),
                 '%s%sD%s' % (Fore.RED,Style.BRIGHT, Style.NORMAL),
                 '%sN' % Fore.MAGENTA,
                 '%s%sS' % (Style.DIM,Fore.RED),
                 '%s%sH' % (Style.DIM,Fore.RED),
                 '%s%sP' % (Style.DIM,Fore.RED),
                 '%s%s=' % (Style.DIM,Fore.GREEN),
                 '%s%sX' % (Style.DIM,Fore.BLUE)]
                [operation]
                , length, Style.RESET_ALL))

        multiMapping = ''
        try:
            multiMappingCount = int(self.alignment.get_tag('NH:i'))

            if multiMappingCount==1:
                multiMapping= Fore.GREEN + 'Unique map' + Style.RESET_ALL
            elif multiMappingCount>1:
                multiMapping= Fore.WHITE + Back.RED + ('Mapping to at least %s locations' % multiMappingCount) + Style.RESET_ALL
        except:
            pass

        print( ('%sReference:%s %s %s %s MQ:%s %s' % ( Style.BRIGHT, Style.RESET_ALL, self.alignment.reference_name, "".join(annotatedCigar), ("FWD" if not self.alignment.is_reverse else "REV"), self.alignment.mapping_quality, multiMapping))+ (' MULTIMAPPING' if self.alignment.is_secondary else '' )  )
        print(('%s%s %s %s (%s%s)' % (
            Style.DIM+refStartStr+Style.RESET_ALL,
            " "*(labelSize-len(refStartStr)),
            "".join(reference),
            (Style.DIM+str(rEnd)+Style.RESET_ALL),
            "" if self.alignment.is_reverse else "+",
            (Style.DIM+str(rEnd-rStart)+Style.RESET_ALL))))
        print(('%s%s' % (" "*(labelSize+1), "".join(alignmentSymbols))))
        print(('%s%s %s %s (%s)' % (
            " "*(labelSize-len(qStartStr)),
            Style.DIM+qStartStr+Style.RESET_ALL,
            "".join(query),
            (Style.DIM+str(self.alignment.query_alignment_end)+Style.RESET_ALL),
            (Style.DIM+str(self.alignment.query_alignment_end-self.alignment.query_alignment_start)+Style.RESET_ALL),
            )))
        print(('%s%s' % (" "*(labelSize+1), "".join(qualitieSymbols))))
        print(('%sQuery:%s %s' % (Style.BRIGHT, Style.RESET_ALL,self.alignment.query_name)))

def getLevenshteinDistance(source, target, maxDist=99999999999999):
    if len(source) < len(target):
        return getLevenshteinDistance(target, source)

    # So now we have len(source) >= len(target).
    if len(target) == 0:
        return len(source)

    # We call tuple() to force strings to be used as sequences
    # ('c', 'a', 't', 's') - numpy uses them as values by default.
    source = np.array(tuple(source))
    target = np.array(tuple(target))

    # We use a dynamic programming algorithm, but with the
    # added optimization that we only need the last two rows
    # of the matrix.
    previous_row = np.arange(target.size + 1)

    for s in source:
        # Insertion (target grows longer than source):
        current_row = previous_row + 1

        # Substitution or matching:
        # Target and source items are aligned, and either
        # are different (cost of 1), or are the same (cost of 0).
        current_row[1:] = np.minimum(
                current_row[1:],
                np.add(previous_row[:-1], target != s))

        # Deletion (target grows shorter than source):
        current_row[1:] = np.minimum(
                current_row[1:],
                current_row[0:-1] + 1)

        previous_row = current_row
        if( np.amin(previous_row)>maxDist):
            return(None)


    return previous_row[-1]


def cigarStringToDict(cigarString):
    if cigarString==None:
        cs = ''
    else:
        cs = cigarString
    cigarStringDict = {}
    for symbol in ['M','D','I','S']:
         # Extract the integer which is found before the symbol supplied
        # returns 0 when the symbol is not found as the sum of the list will be zero
        matches = re.findall('[0-9]*%s' % symbol, cs)
        cigarStringDict[symbol] = sum( [ int(x.replace(symbol,'')) for x in matches])

    return( cigarStringDict)


def distanceLineToPoint(x1, y1, x2, y2, x0, y0):

    return( float( abs( (y2-y1)*x0 - (x2 - x1 )*y0 + x2*y1 - y2*x1 ) ) / (math.sqrt( math.pow(y2-y1,2) + math.pow(x2-x1,2))) )


def fixGraphML(graphmlPath):
    with open(graphmlPath) as f:
        lines = []
        for line in f:
            lines.append( line.replace('"d3"','"fixedD3"') )
    if lines!=None:
        with open(graphmlPath,'w') as f:
            for line in lines:
                f.write(line)


def getHammingDistance(s1, s2,maxDist=999999):
    d = 0
    for index,base in enumerate(s1):
        d+=(base!=s2[index] and s2[index]!="N" and base!="N")
        if d>maxDist:
            return(maxDist+1)
    return(d)


def getHammingIndices( seqA, seqB, maxIndices=999999 ):

    indices = []
    for index,base in enumerate(seqA):
        if len(seqB)<index:
            break
        if base!=seqB[index]:
            indices.append(index)
            if indices>=maxIndices:
                return(indices)

    return(indices)

def humanReadable(value, targetDigits=2,fp=0):

        if value == 0:
                return('0')

        baseId = int(math.floor( math.log10(float(value))/3.0 ))
        suffix = ""
        if baseId==0:
                sVal =  str(round(value,targetDigits))
                if len(sVal)>targetDigits and sVal.find('.'):
                        sVal = sVal.split('.')[0]

        elif baseId>0:

                sStrD = max(0,targetDigits-len(str( '{:.0f}'.format((value/(math.pow(10,baseId*3)))) )))


                sVal = ('{:.%sf}' % min(fp, sStrD)).format((value/(math.pow(10,baseId*3))))
                suffix = 'kMGTYZ'[baseId-1]
        else:

                sStrD = max(0,targetDigits-len(str( '{:.0f}'.format((value*(math.pow(10,-baseId*3)))) )))
                sVal = ('{:.%sf}' %  min(fp, sStrD)).format((value*(math.pow(10,-baseId*3))))
                suffix = 'mnpf'[-baseId-1]

                if len(sVal)+1>targetDigits:
                        # :(
                        sVal = str(round(value,fp))[1:]
                        suffix = ''


        return('%s%s' % (sVal,suffix))




##Pathtools
def decodePath( path, spacers={'flag': ',', 'param':'_', 'kv':'=', 'invalid':['=','_']}):


        invalidFileNameCharacters = spacers['invalid']
        flagSpacer = spacers['flag']
        paramSpacer = spacers['param']
        keyValueSpacer = spacers['kv']

        decodedPath = {
            'runId':None,
            'step':-1,
            'name':'unknown',
            'parameters':{},
            'flags' : [],
            'extension':''
        }

        parts = os.path.splitext(path)
        basePath = parts[0]
        decodedPath['extension'] = parts[1]

        #Check if the anlysis base directory is in the path:
        parts = basePath.split('/')
        if 'analysis_' in parts[0]:
            decodedPath['runId'] = parts[0]
            toDecode = parts[1]
        else:
            toDecode = path

        keyValuePairs = basePath.split(paramSpacer)
        for pair in keyValuePairs:
            keyValue = pair.split(keyValueSpacer)
            if len(keyValue)!=2:
                print(('Parameter decoding failed: %s, this probably means input files contain invalid characters such as %s' % (pair, ", ".join(invalidFileNameCharacters))))
            else:
                key = keyValue[0]
                value = keyValue[1]
                if key!='flags' and key!='name':
                    decodedPath['parameters'][key] = value
                elif key=='flags':
                    decodedPath['flags'] = value.split(flagSpacer)
                elif key=='name':
                    decodedPath['name'] = value

        return(decodedPath)


def invBase(nucleotide):

        if nucleotide=='A':
            return(['T','C','G'])
        if nucleotide=='C':
            return(['T','A','G'])
        if nucleotide=='T':
            return(['C','A','G'])
        if nucleotide=='G':
            return(['C','A','T'])
        return([''])

def encodePath(step, name, params, extension, spacers={'flag': ',', 'param':'_', 'kv':'=', 'invalid':['=','_']}):

    invalidFileNameCharacters = spacers['invalid']
    flagSpacer = spacers['flag']
    paramSpacer = spacers['param']
    keyValueSpacer = spacers['kv']

    #Todo get keyValueSpacer in here
    encodedPath = 'step={:0>2d}'.format(step)
    encodedPath += '_name=%s' % name
    for parname in params:
        if parname.lower()!="flags":
            if not( any(i in parname for i in invalidFileNameCharacters) or any(i in params[parname] for i in invalidFileNameCharacters) ):
                encodedPath+='%s%s%s%s' % (paramSpacer, parname, keyValueSpacer, params[parname])
            else:
                raise ValueError('Parameter encoding failed: %s %s contains at least one invalid character' % (parname, params[parname]) )

    #Todo check flags for mistakes
    if 'flags' in params:
        encodedPath+='%sflags%s%s' % (paramSpacer, keyValueSpacer,flagSpacer.join(params['flags']))

    if extension!=None and extension!=False:
        encodedPath = '%s.%s' % (encodedPath, extension)

    return(encodedPath)


def locateIndel( seqA, seqB, verbose=False ):

    alignment = pairwise2.align.globalmx(seqA, seqB, 2, -1, one_alignment_only=True)[0]

    #pairwise2.format_alignment(*alignment)
    indelLocA = -1
    indelLocB = -1

    if alignment:
        align1 = alignment[0]
        align2 = alignment[1]

        localAlign1 = align1.strip('-')
        localAlign2 = align2.strip('-')

        dels = min( localAlign1.count('-'), localAlign2.count('-') )
        ins =  max( localAlign1.count('-'), localAlign2.count('-') )
        if (ins==1) and dels==0:

            seqNHasIndel = localAlign2.count('-')==0 #bool: 0: seqA has indel, 1: seqB has indel

            #Find indel location:

            if seqNHasIndel:
                indelLocA = align1.find('-')
                if verbose:
                    print(('seqA:%s\nseqB:%s'% (align1,align2)))
                    print((' '*(indelLocA+5)+'^'))
            else:
                indelLocB = align2.find('-')
                if verbose:
                    print(('seqA:%s\nseqB:%s'% (align1,align2)))
                    print((' '*(indelLocB+5)+'^'))


            #print(alignment[2]) is score

        return(indelLocA,indelLocB)


class Histogram():
    def __init__(self):
        self.data = {}

    def addCount(self,idx, count=1):
        if not idx in self.data:
            self.data[idx] = 0
        self.data[idx]+=count

    def write(self, path):
        with open(path,'w') as f:
            for index in self.data:
                f.write('%s;%s\n' % (index, self.data[index]))



globalWritingUpdates = False
class InternalTool():
    def __init__(self, screenName = "unnamed"):
        self.verbose = True
        self.screenName = screenName

    def setVerbosity(self, verbosity):
        self.verbose = verbosity

    def __priorRunningPrint(self):
        if self.isUpdating():
            print('')
            self.setUpdating(False)

    def isUpdating(self):
        global globalWritingUpdates
        if globalWritingUpdates:
            return(True)

    def setUpdating(self, boolean):

        global globalWritingUpdates
        globalWritingUpdates = boolean


    def log(self, message):
        if self.verbose:
            self.__priorRunningPrint()
            print((Fore.GREEN + self.screenName + ', INFO: ' + Style.RESET_ALL + message))

    def warn(self, message):
        self.__priorRunningPrint()
        print((Fore.RED + self.screenName + ', WARN: ' + Style.RESET_ALL + message))

    def softWarn(self, message):
        self.__priorRunningPrint()
        print((Fore.YELLOW + self.screenName + ', SOFT WARN: ' + Style.RESET_ALL + message))



    def update(self, message):
        print('\r' + Fore.GREEN + self.screenName + ', UPDATE: ' + Style.RESET_ALL +message, end=' ')
        self.setUpdating(True)


    def info(self, message):
        self.log(message)


    def fatal(self, message):
        self.__priorRunningPrint()
        print((Fore.RED + self.screenName + ', FATAL ERROR: ' + Style.RESET_ALL + message))
        exit()


class ExternalTool(InternalTool):
    def __init__(self, screenName = "unnamed"):
        InternalTool.__init__(self)
        #Locations to look for binaries when they are not found on the PATH


    def setExecutablePath(self, path):
        self.executable = path

    #Only works in unix and cygwin
    def isAvailable(self):
        available = subprocess.call("type " + self.executable.split()[0], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
        return(available or os.path.isfile(self.executable))


class FastqToFasta(ExternalTool):
    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('fastq_to_fasta')
        self.screenName = "FASTQ_TO_FASTA"

    def execute(self, fqPath, fastaPath):
        os.system('%s -n -i %s -o %s'  % (self.executable, fqPath, fastaPath ))


class FastQQualityFilter(ExternalTool):

    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('fastq_quality_filter')
        self.minBaseQuality = 36
        self.sequenceCountSuffix = 'fastqQualityResultCount'
        self.screenName = "FASTQ_QUALITY_FILTER"

    def setMinBaseQuality(self, minBaseQuality):
        self.minBaseQuality = minBaseQuality

    def getCountPath(self,filteredTargetPath):
        return( '%s_%s' % (filteredTargetPath,self.sequenceCountSuffix) )


    def execute(self, fqPath, filteredTargetPath):
        os.system('%s -q %s -i %s -o %s -v -sam > %s' % (self.executable, self.minBaseQuality, fqPath, filteredTargetPath, self.getCountPath()))
        r = 0
        with open(self.getCountPath(), 'w') as f:
            r = int(f.readline())

        return(r)

class ART(ExternalTool):

    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('./art_bin_ChocolateCherryCake/art_illumina')
        self.screenName = "ART"

    def simulateAmpliconReads(self,sequenceAbundanceCounter, prefix, readSize=76, coverage=10, reverseComplementAmplicon=False):

        sequences = []
        print(( sequenceAbundanceCounter.most_common(1)))
        (v,hCount) = sequenceAbundanceCounter.most_common(1)[0]
        for iteration in range(0,hCount):

            for sequence,count in sequenceAbundanceCounter.most_common():
                if count>=hCount:
                    sequences.append(sequence)
                else:
                    break
            hCount-=1


        bpythonSeqs=[]
        for index,sequence in enumerate(sequences):
            if reverseComplementAmplicon :
                bpythonSeqs.append(SeqIO.SeqRecord(Seq(sequence).reverse_complement(), '%s-%s' % (str(index),str(sequence))))
            else:
                bpythonSeqs.append(SeqIO.SeqRecord(Seq(sequence), '%s-%s' % (str(index),str(sequence))))

        fastaPath = getTempFileName('art')+'.fa'
        SeqIO.write(bpythonSeqs, fastaPath, "fasta")

        os.system('%s -amp -i %s -l %s -f %s -o %s -sam -ss HS25' %(self.executable, fastaPath, readSize, coverage, prefix))



    def setMinBaseQuality(self, minBaseQuality):
        self.minBaseQuality = minBaseQuality

    def getCountPath(self,filteredTargetPath):
        return( '%s_%s' % (filteredTargetPath,self.sequenceCountSuffix) )


    def execute(self, fqPath, filteredTargetPath):
        os.system('%s -q %s -i %s -o %s -v > %s' % (self.executable, self.minBaseQuality, fqPath, filteredTargetPath, self.getCountPath()))
        r = 0
        with open(self.getCountPath(), 'w') as f:
            r = int(f.readline())

        return(r)


class NeedleAll(ExternalTool):

    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('needleall')
        self.screenName = "NEEDLE_ALL"


    def execute(self, fqPathA, fqPathB, outputFilePath):
        command = '%s -auto -stdout -asequence "%s" -bsequence "%s" > %s' % (self.executable, fqPathA, fqPathB, outputFilePath)
        print(command)
        os.system(command)
        r = 0

        #This is an inverted distance matrix
        invScores = DistanceMatrix()

        #find max score:
        maxScore = 0
        # To perform T-SNE a normalised matrix is required: the total sum should be 1
        totalScore = 0
        cells = 0
        with open(outputFilePath, 'r') as f:
           for line in f:
                #print(line)
                parts = line.strip().replace('(','').replace(')','').split(' ')
                if len(parts)==4:
                    value = float(parts[3])
                    if value>maxScore:
                        maxScore=value
                    totalScore+=value
                    cells+=1

        #Normalise the scores


        matrixSum = 0
        with open(outputFilePath, 'r') as f:
           for line in f:
                #print(line)
                parts = line.strip().replace('(','').replace(')','').split(' ')
                #print('%s %s %s %s' % (parts[0], parts[1],parts[2],parts[3]))
                if len(parts)==4:
                    if float(parts[3])>0:
                        value = (maxScore-float(parts[3]))
                        matrixSum+= value
                        invScores.setDistance(str(parts[0]), str(parts[1]), value)

        print(('Sum of distance matrix is %s' % matrixSum))
        return(invScores)



class ClustalOmega(ExternalTool):

    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('clustalo')
        self.screenName = "CLUSTAL-OMEGA"

    def alignFasta(self, fastaPath,outputFastaFilePath):
        self.log('aligning %s to %s' % (fastaPath, outputFastaFilePath))
        os.system('%s -i %s -o %s --force' % (self.executable, fastaPath, outputFastaFilePath))


    def alignFastqs(self, fqPath, outputFastaFilePath):
        #self.log('converting %s to fasta' % fqpath)
        f = FastqToFasta()
        tmpFasta = getTempFileName()
        f.execute(fqPath, tmpFasta)
        self.alignFasta(tmpFasta,outputFastaFilePath)
        os.remove(tmpFasta)




class SamTools(ExternalTool):
    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('samtools')
        self.screenName = "SAMTOOLS"

    def samToSortedBam(self, samPath, bamPath, deleteSam=False):
        self.log('Converting SAM to sorted BAM')
        os.system('%s view -bS %s | %s sort - -o %s.bam' % (self.executable, samPath, self.executable, bamPath))
        self.log('Completed conversion')
        if deleteSam:
            os.remove(samPath)

    def bamToSortedBam(self, bamPath, sortedBamPath):
        self.log('Converting BAM to sorted BAM')
        os.system('%s sort -o %s.bam %s' % (self.executable, sortedBamPath, bamPath))
        self.log('Completed conversion')


    def index(self, bamPath):
        self.log('Indexing BAM file')
        os.system('%s index %s' % (self.executable, bamPath))
        self.log('Completed indexing')

    def readBamToDictOld(self,bamFilePath, columns={'cigar':5, 'left':3,'scarSequence':9}, index='scarSequence', appendTo={}): # columns:{ targetName:columnIndex, ... }
        retDict = appendTo
        #os.system("samtools view %s | cut -f1,6 | grep -oEi '^.*([0-9]+D){1,1}.*' ")
        #tempfile = './temp/temp_' + str(uuid.uuid1())
        tempfile = getTempFileName()
        self.log("samtools view %s > %s" %(bamFilePath,tempfile))
        os.system("samtools view %s > %s" %(bamFilePath,tempfile))
        with open(tempfile) as f:
            #Go over all lines in the BAM file
            for line in f:
                #Extract all columns:
                parts = line.replace('\n','').split()
                #Extract the primary key value:
                indexKey = parts[columns[index]]
                #Initialise the existence of the primary key value in the return dictionary
                if not indexKey in retDict:
                    retDict[indexKey] = {}
                #Iterate over the keysm except the index
                for key in columns:
                    if key is not index:
                        #Get the value of the key on this line in the BAM file
                        v = parts[columns[key]]
                        if not v in retDict[indexKey]:
                            retDict[indexKey][v] = 1
                        else:
                            retDict[indexKey][v] += 1
        os.remove(tempfile)
        return(retDict)

    def readBamToDict(self,bamFilePath, columns={'cigar':5, 'id':0,'left':3,'scarSequence':9}, index='scarSequence', appendTo={}): # columns:{ targetName:columnIndex, ... }
        retDict = appendTo
        #os.system("samtools view %s | cut -f1,6 | grep -oEi '^.*([0-9]+D){1,1}.*' ")
        tempfile = './temp/temp_' + str(uuid.uuid1())
        self.log("samtools view %s > %s" %(bamFilePath,tempfile))
        os.system("samtools view %s > %s" %(bamFilePath,tempfile))
        with open(tempfile) as f:
            #Go over all lines in the BAM file
            for line in f:
                #Extract all columns:
                parts = line.replace('\n','').split()
                #Extract the primary key value:
                indexKey = parts[columns[index]]
                #Initialise the existence of the primary key value in the return dictionary
                if not indexKey in retDict:
                    retDict[indexKey] = {}
                #Iterate over the keysm except the index
                for key in columns:
                    if key is not index:
                        #Get the value of the key on this line in the BAM file
                        if columns[key]<len(parts):
                            v = parts[columns[key]]
                            retDict[indexKey][key] = v

        try:
            os.remove(tempfile)
        except:
            print(('Temporary file could not be removed %s' % tempfile))
            pass
        return(retDict)


# Make sure to index the reference genome!
class BWA_MEM(ExternalTool):
    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('bwa mem')
        self.screenName = "BWA MEM"
        self.fixedFlags = []

    def setReferencePath(self, referenceGenomePath):
        self.referenceGenomePath = referenceGenomePath

    def execute(self, fqPath, samPath, readGroupInfo=None, secondSplitHits=False, getSecondaryAlignments=False):
        if not isinstance(fqPath, str):
            fqPath = " ".join(fqPath)

        if self.referenceGenomePath!=None:
            self.log('Aligning %s to %s' %(fqPath, self.referenceGenomePath))

            flags=[]+self.fixedFlags
            if secondSplitHits:
                flags.append("-M")
            if readGroupInfo!=None:
                flags.append("-R '%s'" % readGroupInfo)
            if getSecondaryAlignments:
                flags.append("-a")

            self.log("Flags: %s" % (" ".join(flags)))
            cmd = '%s %s %s %s > %s' % (self.executable, " ".join(flags),  self.referenceGenomePath, fqPath, samPath)
            self.log(cmd)
            os.system(cmd)
        else:
            self.warn('Warning; no reference genome specified for BWA MEM, canceled alignment.')

# Make sure to index the reference genome!
class STAR(ExternalTool):
    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('STAR')
        self.screenName = "STAR"
        self.threads = 4

    def index(self):

        if not os.path.exists(self.refDirectory):
            try:
                os.mkdir(self.refDirectory )
            except:
                self.fatal('Could not create index directory on %s' % self.refDirectory)

        os.system('%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %s' %
                  (self.executable, self.refDirectory,  self.referenceGenomePath, self.threads))



    def setReferencePath(self, referenceGenomePath):
        self.referenceGenomePath = referenceGenomePath
        self.refDirectory = '%s_STAR_INDEX' % os.path.basename(self.referenceGenomePath)
        if not os.path.exists(self.refDirectory):
            self.info("Reference index %s does not exist." % referenceGenomePath)
            self.index()


    def execute(self, fqPath, outputPath):
        os.system('%s --genomeDir %s --readFilesIn %s --runThreadN %s --outFileNamePrefix %s' %
                (self.executable, self.refDirectory, fqPath, self.threads, outputPath))


    #!!!sjdbOverhang  should match the length of the reads minus one!!!
    def pass2alignment(self, SJouttabPath, sjdbOverhang=74):

        base = os.path.basename(SJouttabPath)

        self.p2refDirectory = '%s_P2_STAR_INDEX' % os.path.basename(SJouttabPath)
        if not os.path.exists(self.p2refDirectory):
            try:
                os.mkdir(self.p2refDirectory )
            except:
                self.fatal('Could not create index directory on %s' % self.p2refDirectory)


    def mapSingleEnd(self, fqPath, samPath):
        if not isinstance(fqPath, str):
            fqPath = " ".join(fqPath)

        if self.referenceGenomePath!=None:
            self.log('Aligning %s to %s' %(fqPath, self.referenceGenomePath))
            os.system('%s --genomeDir %s --readFilesIn %s --runThreadN %s' %
                (self.executable, self.refDirectory, fqPath, self.threads))
        else:
            self.warn('Warning; no reference genome specified for STAR, canceled alignment.')

class GSNAP(ExternalTool):
    def __init__(self):
        ExternalTool.__init__(self)
        self.setExecutablePath('gsnap')
        self.screenName = "GSNAP"
        #self.cmd=  'sudo gsnap -d lintrace analysis_bdefa5bc-6b5c-11e5-a487-0800279bd60a_fullScarList_Protocol1-DNA-seq.fq -t 4 -A sam > analysis_bdefa5bc-6b5c-11e5-a487-0800279bd60a_fullScarList_Protocol1-DNA-seq.sam'

    #Reference genome path should be the name of the reference build folder
    def setReferencePath(self, referenceGenomePath):
        self.referenceGenomePath = referenceGenomePath

    def execute(self, fqPath, samPath):
        if self.referenceGenomePath!=None:
            self.log('Aligning %s to %s' %(fqPath, self.referenceGenomePath))
            os.system('%s -d %s %s --format=sam > %s' % (self.executable, self.referenceGenomePath, fqPath, samPath))
        else:
            self.warn('Warning; no reference genome specified for GSNAP, canceled alignment.')


def getTempFileName(centerFix=""):
    if not os.path.exists('./temp'):
        os.mkdir('./temp')
    return('./temp/temp_' + centerFix + str(uuid.uuid1()))

def mapReads(readsFile, baseName, mapper, samTools=None,readGroupInfo=None, secondSplitHits=False):

    if samTools==None:
        samTools = SamTools()

    basePath = '%s/%s' % ( '.', baseName)
    samPath = '%s.sam' %( basePath)
    bamPath = '%s.bam' %( basePath)


    if readGroupInfo==None and secondSplitHits==False:
        mapper.execute(readsFile, samPath)
    else:
        mapper.execute(readsFile, samPath, readGroupInfo=None, secondSplitHits=False)

    samTools.samToSortedBam(samPath, basePath, True)
    samTools.index(bamPath)
    return({'samPath':samPath, 'bamPath':bamPath})

def mapSequencesToDict(mapper, sequences):

    bpythonSeqs=[]
    for index,sequence in enumerate(sequences):
        bpythonSeqs.append(SeqIO.SeqRecord(Seq(sequence), str(index)))

    fastaPath = getTempFileName('bwa')+'.fa'
    SeqIO.write(bpythonSeqs, fastaPath, "fasta")
    bPath = getTempFileName('samtools')
    samPath = bPath+'.sam'
    bamPath = bPath+'.bam'
    mapper.execute(fastaPath, samPath)
    samtools = SamTools()
    samtools.samToSortedBam(samPath, bPath, True)
    samtools.index(bamPath)
    d = samtools.readBamToDict(bamPath)
    return(d)

def mapDictOfSequences(mapper, sequenceDict, bPath=None):

    bpythonSeqs=[]
    for seqId in sequenceDict:
        bpythonSeqs.append(SeqIO.SeqRecord(Seq(sequenceDict[seqId]), str(seqId)))

    fastaPath = getTempFileName('bwa')+'.fa'
    SeqIO.write(bpythonSeqs, fastaPath, "fasta")
    if bPath==None:
        bPath = getTempFileName('samtools')
    samPath = bPath+'.sam'
    bamPath = bPath+'.bam'
    mapper.execute(fastaPath, samPath)
    samtools = SamTools()
    samtools.samToSortedBam(samPath, bPath, False)
    samtools.index(bamPath)
    return(bamPath)

def mapReadsToDict(mapper, fqPath, index='id',basePath=None, readGroupInfo=None, secondSplitHits=False):
    if basePath==None:
        bPath = getTempFileName('samtools')
    else:
        bPath=basePath
    samPath = bPath+'.sam'
    bamPath = bPath+'.bam'
    mapper.execute(fqPath, samPath,  readGroupInfo, secondSplitHits)
    samtools = SamTools()
    samtools.samToSortedBam(samPath, bPath, True)
    samtools.index(bamPath)
    r = {}
    d = samtools.readBamToDict(bamPath, index=index, appendTo=r)
    return(d)

# Reference to this function: http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
# This method is the quickest way to count lines in big files. (Quicker than wc..)
def fileLineCount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )


# Deprecated: should use Numpy sparse matrix
class DistanceMatrix():
    def __init__(self, maxDist=1000000000):
        self.data = {}
        self.keys = set()
        self.maxDist = maxDist

        self.finalKeying = True # Do not store keys on every update

    def setDistances(self,fromList, toList, distanceList ):
        # Add keys:
        if not self.finalKeying:
            self.keys.update(fromList+toList)

        for index,a in enumerate(fromList):
            self._update(a,toList[index], distanceList[index])

    def getNumpyArray(self):
        #Make sure the keys are properly indexed
        self._performKeying()

        #Initialise empty matrix:
        narray = np.zeros((len(self.keys), len(self.keys)))

        for aIndex,keyA in enumerate(self.keys):
            for bIndex,keyB in enumerate(self.keys):
                narray[aIndex,bIndex] = self.getDistance(keyA,keyB)

        return(narray)


    def _update(self, a,b,distance):
        if distance<=self.maxDist:
            if (not a in self.data) and (not b in self.data):
                self.data[a] = {}
                self.data[a][b] = distance
                return(True)
            if a in self.data and b in self.data[a]:
                self.data[a][b] = distance
                return(True)
            if b in self.data and a in self.data[b]:
                self.data[b][a] = distance
                return(True)


        # A already exists but B does not: (And b also does not exists as main entry)
        if a in self.data and not b in self.data and not b in self.data[a]:
            self.data[a][b] = distance
            return(True)
        elif b in self.data:
            if not a in self.data[b]:
                self.data[b][a] = distance
                return(True)
        return(False)


    def setDistance(self, a,b,distance ):
        if not self.finalKeying:
            if a not in self.keys and b not in self.keys:
                self.keys.update([a,b])
            else:
                if a not in self.keys:
                    self.keys.add(a)
                if b not in self.keys:
                    self.keys.add(b)

        return( self._update(a,b, distance ))


    def getDistance(self, a,b ):
        if a in self.data and b in self.data[a]:
            return(self.data[a][b])
        if b in self.data and a in self.data[b]:
            return(self.data[b][a])
        return(self.maxDist)


    def _performKeying(self):
        #Retrieve all keys
        self.keys = set()
        addKeys = []
        self.keys.update(list(self.data.keys()))
        for a in self.data:
            for b in self.data[a]:
                if b not in self.keys:
                    addKeys.append(b)
        self.keys.update(addKeys)

    def writeToFile(self, path):

        if self.finalKeying:
            self._performKeying()

        separator = ';'
        with open(path, 'w') as f:
            values= ['name']
            for keyA in self.keys:
                values.append(keyA)
            f.write('%s\n' % (separator.join(values)))
            for keyA in self.keys:
                values = []
                for keyB in self.keys:
                    values.append(str(self.getDistance(keyA, keyB)))
                f.write('%s%s%s\n' % (keyA, separator,separator.join(values)))


    def loadFromFile(self, path,separator=None):
        with open(path) as f:
            idx = 0
            header = {}
            for line in f:
                parts = line.rstrip().split(separator)
                if idx==0:
                    for keyId,keyName in enumerate(parts):
                        header[keyId] = keyName
                else:

                    for partIndex, part in parts:
                        if partIndex==0:
                            keyA = part
                        else:
                            keyB=header[partIndex]
                            self.setDistance(keyA,keyB, float(part))

                idx+=1
