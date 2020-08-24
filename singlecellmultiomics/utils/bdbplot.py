from lxml import etree
import math
import collections
from collections import Counter
from collections import OrderedDict
import numpy as np
from singlecellmultiomics.utils import bdbbio
import os
import matplotlib.cm
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from colorama import Fore #,Back, Style
from colorama import Back
from colorama import Style
from colorama import init
import scipy
import scipy.cluster
import time
import itertools

init(autoreset=True)



#Convert a nested dictionary to a matrix
# ({'A':{'1':2}, 'B':{'1':3, '2':4}}) will become
#(array([[  2.,  nan],
#        [  3.,   4.]]), ['A', 'B'], ['1', '2'])


def interpolateBezier( points, steps=10, t=None):
    if len(points)==3:
        mapper = lambda t,p: (1-t)**2 * p[0] + 2*(1-t)*t*p[1] + t**2*p[2]
    elif len(points)==4:
        mapper = lambda t,p: (np.power( (1-t),3)*p[0] +\
         3* np.power((1-t),2) *t *p[1] +\
         3*(1-t)*np.power(t,2)*p[2] +\
         np.power(t,3)*p[3])

    if t is not None:
        return   mapper(t, [q[0] for q in points]), mapper(t, [q[1] for q in points])
    xGen = ( mapper(t, [q[0] for q in points]) for t in np.linspace(0, 1, steps) )
    yGen = ( mapper(t, [q[1] for q in points]) for t in np.linspace(0, 1, steps) )

    return zip(xGen, yGen)

def interpolateBezierAngle(points, t, ds=0.001):
    x0, y0 = interpolateBezier(points, t=t-ds)
    x1, y1 = interpolateBezier(points, t=t+ds)
    return np.arctan2( y1-y0, x0-x1)


def initMatrix(rowNames,columnNames, mtype="obj"):
    if mtype=="obj":
        matrix = np.empty( (len(rowNames), len(columnNames)), dtype=object)
    elif mtype=="npzeros":
        matrix = np.zeros( (len(rowNames), len(columnNames)))
    return(matrix)

def nestedDictionaryToNumpyMatrix( nestedDictionary, setNan=True, mtype="obj", transpose=False, indicateProgress=False):

    rowNames = sorted( list(nestedDictionary.keys()), key=int )
    columnNames = set()
    for key in nestedDictionary:
        columnNames.update( set(nestedDictionary[key].keys() ))
    columnNames = sorted( list(columnNames) )

    keys = list(columnNames)

    if ':' in keys[0]:
        sargs = np.argsort( [ int(k.split(':')[1]) for k in keys] )
        print(sargs)

        columnNames = [ keys[index] for index in sargs ]
        print(columnNames)

    matrix = initMatrix(rowNames,columnNames, mtype)

    if setNan:
        matrix[:] = np.nan

    prevTime = time.time()
    for rowIndex,rowName in enumerate(rowNames):
        if indicateProgress and (time.time()-prevTime)>1:
            prevTime = time.time()
            print("Matrix creation progress: %.2f%%" % (100.0*rowIndex/len(rowNames)))

        for colIndex,colName in enumerate(columnNames):
            try:
                matrix[rowIndex,colIndex] = nestedDictionary[rowName][colName]
            except:
                pass

    if transpose:
        matrix = matrix.transpose()
        columnNames, rowNames =  rowNames, columnNames

    if indicateProgress:
        print("Matrix finished")
    return( (matrix, rowNames, columnNames) )

def pruneNonUniqueColumnsFromMatrix(matrix,rows,columns, minInstances=1,minOccurence=1):
    colsToKeep = []
    for columnIndex in range(matrix.shape[1]):
        if len( set(np.unique( matrix[:,columnIndex].astype(str) ))-set( ["nan"]) )>minInstances:

            cnts = Counter( list(matrix[:,columnIndex].astype(str)) )

            counts = Counter({k: cnts[k] for k in  cnts if cnts[k] >= minOccurence})
            del counts['nan']

            if len(counts.values())>minInstances:
                colsToKeep.append(columnIndex)

    matrix = matrix[:,colsToKeep]
    columns = np.array(columns)[colsToKeep]
    return(matrix, rows, columns)



# Convert dictionary of tuples to a numpy matrix
def tupleAnnotationsToNumpyMatrix( originRowNames, originColNames, tuples, setNan=True, mtype="obj" ):

    m = initMatrix(originRowNames,originColNames, mtype) #np.zeros( (len(originRowNames), len(originColNames)))
    if setNan:
        m[:] = np.nan
    for tup in tuples:
        value = tuples[tup]
        if tup[0] in originRowNames and tup[1] in originColNames:
         m[originRowNames.index(tup[0]), originColNames.index(tup[1])  ] = value
    return(m)

def dictAnnotationsToNumpyMatrix( originRowNames, originColNames, dictionary, mtype="obj" ):
    tuples = {}
    for rowKey in dictionary:
        for columnKey in dictionary[rowKey]:
            tuples[ (rowKey, columnKey) ] = dictionary[rowKey][columnKey]
    return(tupleAnnotationsToNumpyMatrix(originRowNames, originColNames, tuples, mtype=mtype))


def getSomeColors(n):
    return( plt.cm.Set1(np.linspace(0, 1, n)) )



def _ipol(a, b, first, last,  interpolateValue):
    #Due to floating point rounding errors the interpolate value can be very close to last,
    # it is ok to return last in those cases
    if last>first and interpolateValue>=last:
        return(b)
    if last<first and interpolateValue>=first:
        return(a)

    y_interp = scipy.interpolate.interp1d([first, last], [a,b])
    return( y_interp(interpolateValue) )

def interpolate(interpolateValue,  colorScaleKeys, nodeColorMapping):
        #Seek positions around value to interpolate
        first = colorScaleKeys[0]
        index = 0
        last = first
        for value in colorScaleKeys:
            if value>=interpolateValue:
                last = value
                break
            else:
                first = value
            index+=1
        if value==interpolateValue:
            return(nodeColorMapping[value])

        #Do interpolation
        colorA = nodeColorMapping[first]
        colorB = nodeColorMapping[last]
        dx = last-first

        # Check out of bounds condition
        if interpolateValue< first:
            return(colorA)
        if interpolateValue>last:
            return(colorB)



        return( _ipol(colorA[0], colorB[0], first, last, interpolateValue), _ipol(colorA[1], colorB[1], first, last, interpolateValue), _ipol(colorA[2], colorB[2], first, last, interpolateValue))



def plotFeatureSpace(features, classLabels, featureNames, path=None, bins = 50, title=None):

    print(Fore.GREEN + "Feature space plotter:")
    #1d mode:
    print(features.shape)
    classAbundance = Counter(classLabels)
    print(classAbundance)

    classes = set(classLabels)
    classColors = ['#FF6A00','#0066FF','#FF33FF','#666666']
    classColors = [ tuple(float(int(hexColor.replace('#','')[i:i+2], 16))/255.0 for i in (0, 2 ,4)) for hexColor in classColors]

    if features.shape[1]==1:
        print("Performing 1-D density plot of %s samples" % ( features.shape[0]))

        plt.close('all')

        histStart = features.min()
        histEnd = features.max()
        if histStart == histEnd:
            histEnd += 1
            histStart -= 1
        precision = (histEnd-histStart)/bins


        print("Histogram will be plotted from %s to %s " % (histStart, histEnd))
        fig, ax = plt.subplots() #figsize=(120, 10))

        print(classColors)
        for classIndex,className in enumerate(list(classes)):
            boolList = np.array(classLabels)==np.array(className)
            classSize = classAbundance[className]
            ax.hist(
                features[boolList,0],
                np.arange(histStart,histEnd+precision,precision),
                normed=False, fc= (classColors[classIndex]+ (0.5,)),
                ec=classColors[classIndex],
                lw=1.5, histtype='stepfilled',
                label='%s[%s]' % (className,classSize)
                )

        plt.ylabel("Density")
        plt.xlabel(featureNames[0])
        if title is not None:
            plt.title(title)
        ax.legend(loc='upper right')
        #plt.yscale('log', nonposy='clip')
        if path is None:
            plt.show()
        else:
            plt.savefig(path)
        return(True)


    if features.shape[1]==2:
        #2d
        print("Performing 2-D density plot of %s samples" % ( features.shape[0]))
        fig, ax = plt.subplots()

        for classIndex,className in enumerate(list(classes)):
            #print(np.where( classLabels==className, features ))
            boolList = np.array(classLabels)==np.array(className)
            plt.plot( features[boolList,0], features[boolList,1], ".",label='%s[%s]' % (className,sum(boolList)), c= (classColors[classIndex] + (0.5,)))

        plt.xlabel(featureNames[0])
        plt.ylabel(featureNames[1])
        plt.legend(loc="lower right")
        plt.tight_layout()
        try:
            plt.savefig(path)

        except Exception as e:
            print(e)
        return(True)
    print(Fore.RED + "Invalid amount of dimensions for feature space plotting (%s)" % features.shape[1])

def matplotHeatmap( D, YC, figsize=(10,10), clust=True, xLab=None, yLab=None, show=True, colormap=plt.cm.YlGnBu_r, colorbarLabel=None ):
    plt.rcParams["axes.grid"] = False
    import scipy
    import scipy.cluster.hierarchy as sch
    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=figsize)
    if not clust:
        idx1 = range(0, len(YC))
        idx1 = range(0, len(YC))

    if clust:
        ax1 = fig.add_axes([0.09,0.1,0.2,0.6])

        L = sch.linkage(D, method='centroid')
        Z1 = sch.dendrogram(L, orientation='right')
        ax1.set_xticks([])
        ax1.set_yticks([])

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
        Z2 = sch.dendrogram(L)
        ax2.set_xticks([])
        ax2.set_yticks([])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        D = D[idx1,:]
        D = D[:,idx1]

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])

    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=colormap)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    ###

    if xLab is None:
        axmatrix.set_xticks(range(len(YC)))
        axmatrix.set_xticklabels(YC[idx1])
    else:
        axmatrix.set_xticks(range(len(xLab)))
        axmatrix.set_xticklabels(xLab)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()

    plt.xticks(rotation=-90)

    if yLab is None:
        axmatrix.set_yticks(range(len(YC)))
        axmatrix.set_yticklabels(YC[idx1], minor=False)
    else:
        axmatrix.set_yticks(range(len(yLab)))
        axmatrix.set_yticklabels(yLab, minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()


    # Plot colorbar.
    axcolor = fig.add_axes([1.05,0.1,0.02,0.6])
    cbar = fig.colorbar(im, cax=axcolor)
    if colorbarLabel is not None:
        cbar.set_label(colorbarLabel, rotation=270,  labelpad=15)
    if show:
        fig.savefig('dendrogram.png')

def tsnePlot(data, labels=None, components=2, perplexity=30.0, iterations=1000):
    from sklearn.manifold import TSNE
    #from MulticoreTSNE import MulticoreTSNE as TSNE
    model = TSNE(n_components=components, perplexity=perplexity, n_iter=iterations ) #random_state=0, n_jobs=8,
    transformedPoints = model.fit_transform(data.astype(np.float64))


    classes = list(set(labels))
    classColors = getSomeColors(len(classes))
    color = np.array([ classColors[classes.index(label)] for label in labels ])
    nplabels = np.array(labels)
    print("TSNE input:")
    print(data.shape)
    print("TSNE plotting for %s classes " %  len(classes))

    #print("Color mapping is %s" % ",".join(color))
    #Plot the data:
    if components==2:
        fig = plt.figure()
        #plt.style.use('ggplot')
        ax = fig.add_subplot(111)

        print(color)
        #plt.scatter(transformedPoints[:, 0], transformedPoints[:, 1], c=color, cmap=plt.cm.Spectral, s=1, alpha=0.5) #labels=labels,
        for classIndex, className in enumerate(classes):
            classColor = classColors[classIndex]
            print(className)
            print(classColor)
            plt.scatter(transformedPoints[className==nplabels, 0], transformedPoints[className==nplabels, 1],   s=3, alpha=0.9, label=className)  #c=classColor,, cmap=plt.cm.Spectral,

        #plt.axis('tight')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.show()


    elif components==3:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for classIndex,className in enumerate(classes):
            samplesForClass = (className==nplabels)
            ax.scatter(transformedPoints[samplesForClass,0], transformedPoints[samplesForClass,1], transformedPoints[samplesForClass,2], c=classColors[classIndex], label=classIndex)
        ax.legend()
        plt.show()


class BDBcolor():

    def __init__(self, r=0, g=0, b=0, a=1.0 ):

        if str(r)[0]=='#':
            #parse hex colour:
            #parts = r.replace('#','').replace('(','').replace(')','').split(',')
            cleaned = r.replace('#','')

            #parts = [int(i) for i in parts]
            r =  int(cleaned[0:2], 16)
            g =  int(cleaned[2:4], 16)
            b =  int(cleaned[4:6], 16)

        self.r = max(0,min(255,r))
        self.g = max(0,min(255,g))
        self.b = max(0,min(255,b))
        self.a = max(0,min(1.0,a))

    def getRGBStr(self):
        return('rgb(%s,%s,%s)' % (int(self.r), int(self.g), int(self.b)))

    def getRGBAStr(self):
        return('rgba(%s,%s,%s,%s)' % (self.r, self.g, self.b, self.a))


    def getReadableInverted(self):

        hsv = self.getHSV()
        hsv['v'] = (  255.0-( hsv['v'] ) )
        rgb = self.HSVtoRGB(0,0,hsv['v'])
        return( BDBcolor( rgb['r'],rgb['g'],rgb['b'] ))

    def getHSV(self):
        h=0
        s=0
        v=0
        minV = min( self.r, self.g, self.b )
        maxV = max( self.r, self.g, self.b )
        v = maxV
        delta = maxV - minV
        if maxV != 0:
            s = delta / float(maxV)
        else:
            s = 0
            h = -1
            return({'h':h, 's':s, 'v':v})

        if delta==0:
                h = 255
        else:
                if self.r == maxV:
                    h = ( self.g - self.b ) / float(delta)

                else:
                        if self.g == maxV:

                                h = 2.0 + float( self.b - self.r ) / float(delta)
                        else:
                                h = 4.0 + float( self.r - self.g ) / float(delta)


        h *= 60
        if h < 0:
            h += 360

        return({'h':h, 's':s, 'v':v})


    def HSVtoRGB( self, h,s,v ):


        i=0
        f=0
        p=0
        q=0
        t = 0

        if s == 0:
            #Grey
            r = g = b = v
            return({'r':round(r), 'g':round(g), 'b':round(b)})

        h /= 60         # sector 0 to 5
        i = math.floor( h )
        f = h - i           # factorial part of h
        p = v * ( 1 - s )
        q = v * ( 1 - s * f )
        t = v * ( 1 - s * ( 1 - f ) )

        if i==0:
            r = v
            g = t
            b = p
        elif i==1:
            r = q
            g = v
            b = p
        elif i==2:
            r = p
            g = v
            b = t
        elif i==3:
            r = p
            g = q
            b = v
        elif i==4:
            r = t
            g = p
            b = v
        else:
            r = v
            g = p
            b = q

        return({'r':round(r), 'g':round(g), 'b':round(b)})


class BDBPlot():


    def __init__(self):

        # We need to declare the xlink namespace, to create references to things in our own file
        self.xlink =  'http://www.w3.org/1999/xlink'

        NSMAP = {'xlink':self.xlink }

        self.svgTree = etree.Element('svg',nsmap = NSMAP)
        self.svgTree.set('xmlns','http://www.w3.org/2000/svg')
        self.svgTree.set('version','1.2')



        self.root = self.svgTree.getroottree()

        self.nextFilterId = 0
        self.nextObjId = 0
        self.nextTspanId = 0
        #Create definition element
        self.svgTree.append( self.getDefinitionBlock() )
        self.debug = 2 # 2 all

        self.xMin = 0
        self.xMax = 10
        self.yMin = 0
        self.yMax = 10
        self.plotStartX = 30
        self.plotStartY = 30
        self.plotHeight = 400
        self.plotWidth = 600

        self.setWidth(800)
        self.setHeight(1000)
        self.script = ""

    def clear(self):

        toRm = []
        for child in self.svgTree:
            toRm.append(child)
        for child in toRm:
            self.svgTree.remove(child)

        self.svgTree.append( self.getDefinitionBlock() )
        self.nextFilterId = 0
        self.nextObjId = 0

    def getGroup(self, identifier, zIndex=0 ):
        g = etree.Element('g')
        g.set('id', str(identifier))
        self.svgTree.append(g)
        return(g)


    def getTspan(self):
        tspan = etree.Element('tspan')
        tspan.set('id', str(self.nextTspanId))
        self.nextTspanId+=1
        return(tspan)

    def addLegend(self,colorMapping):

        y = 0
        c = self.getYLabelCoord(y)
        yp = c[1]+10 + 80
        for color in colorMapping:


            text = self.getText(str(colorMapping[color]), c[0]+self.plotStartX, yp,BDBcolor(80,80,80,1))
            text.set('text-anchor','begin')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            text.set('fill','%s' % color)
            self.svgTree.append( text )
            yp+=20



    def getGroupColors(self, n):

        if n==1:
            return(['#3770C4'])
        if n==2:
            return(['#3770C4','#66A43E'])
        if n==3:
            return(['#3770C4','#66A43E','#F6853A'])
        if n==4:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A'])
        if n==5:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2'])
        if n==6:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2','#AAD400'])
        if n==7:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2','#AAD400','#9DAC93'])
        if n==8:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2','#AAD400','#9DAC93','#7FCADF'])
        if n==9:
            return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2','#AAD400','#9DAC93','#7FCADF','#D1AC17'])



        return(['#3770C4','#A43E3E','#66A43E','#F6853A','#A33DA2','#AAD400','#9DAC93','#7FCADF','#D1AC17','#000080','#FF0066','#6C5D53'] + ['#333333']*n)


    #Macro to set a title quickly
    def setTitle(self, string, x=None, y=10, size=25, fill='#333333' ):

        centerX = x is None
        if x is None:
            x,_= self.getPlottingCoord(self.xMin + 0.5*(self.xMax - self.xMin), 0)


        text = self.getText(str(string), x,y, fill=fill)
        if centerX:
            text.set('text-anchor','middle')
        else:
            text.set('text-anchor','begin')
        text.set('dominant-baseline','central')
        text.set('font-family','Gill Sans MT')
        text.set('font-size', str(size))
        self.svgTree.append(text)
        return(self)

    def setSubtitle(self, string, x=None, y=40, size=15, fill='#222222'):
        self.setTitle(string,x,y,size,fill)

    def setWidth(self, width):
        self.width = width
        self.svgTree.set('width','%s' % width)

    def setHeight(self, height):
        self.height = height
        self.svgTree.set('height','%s' % height)


    def getDx(self):
        return( float(self.plotWidth )/float((self.xMax - self.xMin)))


    def getDy(self):
        return( float(self.plotHeight )/float((self.yMax - self.yMin)))


    def getPlottingCoord(self, x,y,z=0):

        return( (self.plotStartX + (float(x)/(self.xMax - self.xMin))*self.plotWidth, self.plotHeight+self.plotStartY -(( float(y)/(self.yMax - self.yMin)))*self.plotHeight))


    def getXLabelCoord(self, x):
        return( (self.plotStartX + (float(x)/(self.xMax - self.xMin))*self.plotWidth, self.plotHeight+self.plotStartY+2 ))


    def getYLabelCoord(self, y):
        return( (self.plotStartX, self.plotHeight+self.plotStartY -(( float(y)/(self.yMax - self.yMin)))*self.plotHeight ))

    def getNextObjId(self):
        self.nextObjId+=1
        return(str(self.nextObjId))


    def getDefinitionBlock(self):

        self.defs = etree.Element('defs')
        self.defs.set('id','defs0')
        return(self.defs)


    def filter(self):
        filterDef = etree.Element('filter')
        filterDef.set('id','filter_%s' % (self.nextFilterId))
        self.nextFilterId+=1
        return(filterDef)

    def getAxis(self,hv=0):

        if hv==1:
            p = self.getPath(self.getPathDefinition([self.getPlottingCoord(self.xMin, self.yMin), self.getPlottingCoord(self.xMax, self.yMin)]))
        elif hv==2:
            p = self.getPath(self.getPathDefinition([self.getPlottingCoord(self.xMin, self.yMax),self.getPlottingCoord(self.xMin, self.yMin)]))
        else:
            p = self.getPath(self.getPathDefinition([self.getPlottingCoord(self.xMin, self.yMax),self.getPlottingCoord(self.xMin, self.yMin), self.getPlottingCoord(self.xMax, self.yMin)]))
        return(p)

    def getPathDefinition(self, coordinates, preventAliasing=False ):

        definition = []
        for idx,coordinateTuple in enumerate(coordinates):

            if preventAliasing:
                coordinateTuple = ( round(coordinateTuple[0])+0.5, round(coordinateTuple[1])+0.5 )

            if idx==0:
                definition.append('M%s,%s' % (coordinateTuple[0],coordinateTuple[1]))
            else:
                definition.append('L%s,%s' % (coordinateTuple[0],coordinateTuple[1]))
        return(' '.join(definition))


    def getLinearGradientDefinition(self, tuplesWithStops): #Format: %x, color

        definitionElement =  etree.Element('linearGradient')
        definitionElement.set('id', self.getNextObjId()) #not really needed
        definitionElement.set('x1', tuplesWithStops[0][0])
        definitionElement.set('y1', tuplesWithStops[0][0])
        definitionElement.set('x2', tuplesWithStops[-1][0])
        definitionElement.set('y2', tuplesWithStops[-1][0])

        for i, tup in enumerate(tuplesWithStops):

            stop = etree.SubElement(definitionElement, 'stop')
            stop.set('offset',tup[0])
            stop.set('stop-color',tup[1])
            if i==0:
                stop.set('class','start')
        stop.set('class','stop')
        return(definitionElement)

    def shadow(self, dy=2, dx=2, gaussStd=2,color='rgb(0,0,0)', floodOpacity=0.9,
     width=None, # Width: set a pixel region around the filter to prevent clipping (https://stackoverflow.com/questions/17883655/svg-shadow-cut-off)
     height=None):
        f = self.filter()

        if width is not None:
            f.set('width', str(width))
            f.set('x', '-%s' % (width*0.5))
        if height is not None:
            f.set('height', str(height))
            f.set('y', '-%s' % (height*0.5))

        f.set('color-interpolation-filters','sRGB')
        self.nextFilterId+=1
        #Flood
        flood = etree.SubElement(f, 'feFlood')
        flood.set('result','flood')
        flood.set('flood-color',color)
        flood.set('flood-opacity','%s' % floodOpacity)

        #Composite filter
        composite1 = etree.SubElement(f, 'feComposite')
        composite1.set('in2','SourceGraphic')
        composite1.set('operator','in')
        composite1.set('in','flood')
        composite1.set('result','composite1')

        #Gaussian blur
        gauss = etree.SubElement(f, 'feGaussianBlur')
        gauss.set('result','blur')
        gauss.set('stdDeviation','%s' % gaussStd)

        #Shadow offset
        offset = etree.SubElement(f, 'feOffset')
        offset.set('result','offset')
        offset.set('dy','%s' % dy)
        offset.set('dx','%s' % dx)

        #Final composite filter
        composite2 = etree.SubElement(f, 'feComposite')
        composite2.set('in2','offset')
        composite2.set('operator','over')
        composite2.set('in','SourceGraphic')
        composite2.set('result','composite2')
        return(f)


    def makeInnerShadow(self, shadow ):

        i = etree.SubElement(shadow, 'feComposite')
        i.set('operator','in')
        i.set('in2','SourceGraphic')



    def addDef(self, filterDef, defId=None):
        if defId is not None:
            filterDef.set('id',defId)
        self.defs.append(filterDef)
        return( filterDef.get('id') )

    def hasDef(self, defId):
        return( len(self.defs.findall(".//*[@id='%s']" % defId))>0 )

    def getDef(self,defId):
        if not self.hasDef(defId):
            print(('Definition %s was not found' % defId))
            exit()
        return( self.defs.findall(".//*[@id='%s']" % defId)[0] )

    def warn(self, msg):
        print(('[WARN] %s' % msg))

    def getRectangle(self, x,y, width, height):
        rectangle = etree.Element('rect')
        rectangle.set('id', self.getNextObjId())
        rectangle.set('x',str(x))
        rectangle.set('y',str(y))
        rectangle.set('width',str(width))
        rectangle.set('height',str(height))
        rectangle.set('style',"fill:rgba(100,100,100,1);stroke:#1b1b1b")
        return(rectangle)

    def getImage( self, path, x=None, y=None, width=None, height=None, preserveAspectRatio=None):
        image = etree.Element('image')
        image.set('id', self.getNextObjId())
        if x is not None:
            image.set('x',str(x))
        if y is not None:
            image.set('y',str(y))
        if width is not None:
            image.set('width',str(width))

        if preserveAspectRatio is not None:
            image.set('preserveAspectRatio',str(preserveAspectRatio))

        if height is not None:
            image.set('height',str(height))

        image.set('{%s}href'% self.xlink ,str(path))

        return(image)



    #Modify attribute in style
    def modifyStyleString(self, style, setAttr={},remove=[]):

        attributes = {}

        if style is not None and style.strip()!='':
            parts = style.split(';')

            for part in parts:
                kvPair = part.split(':')
                if len(kvPair)==2:
                    key = kvPair[0]
                    value = kvPair[1]

                    if key not in remove:
                        attributes[kvPair[0]] = kvPair[1]
                else:
                    self.warn('Style parsing %s failed (ignoring)' % part)



            if self.debug>=3:
                print('Style decomposition')
                for key in attributes:
                    print(('%s\t:\t%s' % (key, attributes[key]) ))

        #Roll changes:
        for attribute in setAttr:
            attributes[attribute] = setAttr[attribute]

        #Create new style string
        newStyle = []
        for attr in attributes:
            newStyle.append('%s:%s' % (attr, attributes[attr]))

        return(';'.join(newStyle))

    def modifyStyle(self, element,setAttr={},remove=[]):

        if not 'style' in element.attrib:
            element.set('style','')

        element.set('style', self.modifyStyleString(element.get('style'), setAttr, remove ))


    def setTextRotation(self, element, angle ):
        element.set('transform','rotate(%s, %s, %s)'%(angle,element.get('x'), element.get('y')))



    def humanReadable(self, value, targetDigits=2,fp=0):

        #Float:
        if value<1 and value>0:
            return('%.2f' % value )

        if value == 0.0:
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



    def getText(self, text, x=0, y=0, fill=BDBcolor(), pathId=None):


        textElement =  etree.Element('text')
        textElement.set('id', self.getNextObjId())

        if pathId != None:
            tp =  etree.Element('textPath')
            tp.text = text
            tp.set('{%s}href'%self.xlink,  '#%s' % (pathId))
            tp.set("startOffset", "50%")

            textElement.append(tp)
        else:
            textElement.text = str(text)

        textElement.set('x',str(x))
        textElement.set('y',str(y))

        if type(fill) is str:
            textElement.set('fill',fill)
        else:
            textElement.set('fill',fill.getRGBStr())
        textElement.set('shape-rendering','crispEdges')


        return(textElement)

    def getCenteredText(self,text, x,y, bold=False, fontSize=14, fill='rgba(50,50,50,1)', **kwargs):
        text = self.getText(text,x,y,**kwargs)
        text.set('text-anchor','middle')
        text.set('dominant-baseline','middle')
        text.set('font-family','Helvetica')
        #text.set('font-family','Cambria Math')
        if bold:
            text.set('font-weight', 'bold')
        text.set('font-size', str(fontSize))
        text.set('fill', fill)
        return text

    def addTspan(self, textObject, text=None):
        tspan = self.getTspan()
        textObject.append(tspan)

        #Fill with text if supplied:
        if text is not None:
            tspan.text = text

        return(tspan)
    #superscript
    def addSuper(self, text, superText, offset=-10):

        superElement =  etree.Element('tspan')
        superElement.text = superText
        superElement.set('dy', '%s' % offset)
        text.append(superElement)




    def polarToCartesian(self, centerX, centerY, radius, angleInDegrees = None, angleInRadians = None):

        angleInRadians = (angleInDegrees-90) * math.pi / 180.0 if angleInDegrees is not None else angleInRadians
        return({
          'x': centerX + (radius * math.cos(angleInRadians)),
          'y': centerY + (radius * math.sin(angleInRadians))
        })

    def describeArc(self, x, y, radius, startAngle, endAngle, sweep=0, largeArcFlag=None):

        start = self.polarToCartesian(x, y, radius, endAngle)
        end = self.polarToCartesian(x, y, radius, startAngle)

        if largeArcFlag==None:
            if endAngle - startAngle <= 180:
                largeArcFlag  = "0"
            else:
                largeArcFlag  = "1"

        d = " ".join([str(x) for x in [
            "M", start['x'], start['y'],
            "A", radius, radius, 0, largeArcFlag, sweep, end['x'], end['y']
        ]])

        return(d)

    def describeArcRad(self, x, y, radius, startAngle, endAngle, sweep=0, largeArcFlag=None):

        start = self.polarToCartesian(x, y, radius, angleInRadians=startAngle )
        end = self.polarToCartesian(x, y, radius, angleInRadians=endAngle  )

        if largeArcFlag==None:
            if endAngle - startAngle <= math.pi:
                largeArcFlag  = "0"
            else:
                largeArcFlag  = "1"

        d = " ".join([str(x) for x in [
            "M", start['x'], start['y'],
            "A", radius, radius, 0, largeArcFlag, sweep, end['x'], end['y']
        ]])

        return(d)


    def getCircle(self, centerX, centerY, radius):
        circle = etree.Element('circle')
        circle.set('id', self.getNextObjId())
        circle.set('cx',str(centerX))
        circle.set('cy',str(centerY))
        circle.set('r',str(radius))
        circle.set('style',"fill:none;stroke:#1b1b1b;stroke-width:1.29999995;stroke-linecap:round;stroke-miterlimit:4;stroke-opacity:1;stroke-dasharray:5.2, 5.2;stroke-dashoffset:0")
        return(circle)


    def getPath(self, pathDef):
        path = etree.Element('path')
        path.set('id', self.getNextObjId())
        path.set('d', pathDef)
        path.set('style',"fill:none;stroke:#1b1b1b;stroke-width:1;stroke-linecap:round")
        return(path)


    def dump(self):
        print(( etree.toString(self.root, pretty_print=True)))




    def write(self, path, pretty=False, htmlCallback=None, bodyCallback=None):

        try:
            os.makedirs(os.path.dirname(path),exist_ok=True)
        except:
            pass

        if len(self.script)>0:
            html = etree.Element('html')


            head = etree.Element('head')
            body = etree.Element('body')

            s = etree.Element('script')
            s.set('type','text/javascript')
            s.text = self.script

            jquery = etree.Element('script')
            jquery.set('src','https://ajax.googleapis.com/ajax/libs/jquery/2.2.3/jquery.min.js')
            jquery.text = ' '
            head.append(jquery)


            jqueryUi = etree.Element('link')
            jqueryUi.set('rel', 'stylesheet')
            jqueryUi.set('href','https://ajax.googleapis.com/ajax/libs/jqueryui/1.12.1/themes/smoothness/jquery-ui.css')
            jqueryUi.text = ' '
            head.append(jqueryUi)

            jqueryUi = etree.Element('script')
            jqueryUi.set('src','https://ajax.googleapis.com/ajax/libs/jqueryui/1.12.1/jquery-ui.min.js')
            jqueryUi.text = ' '
            head.append(jqueryUi)

            jqueryColor = etree.Element('script')
            jqueryColor.set('src','http://code.jquery.com/color/jquery.color-2.1.2.js')
            jqueryColor.text = ' '
            head.append(jqueryColor)

            body.append(self.svgTree)
            html.append(head)

            body.append(s)

            if htmlCallback is not None:
                htmlCallback(html)

            if bodyCallback is not None:
                bodyCallback(body)
            html.append(body)

            import html as pyhtml
            if pretty:
                import xml.dom.minidom as minidom
                with open(path, 'w') as f:
                    f.write( pyhtml.unescape(minidom.parseString(etree.tostring(html.getroottree())).toprettyxml(indent=" ").decode('utf-8') ))
            else:

                #html.getroottree().write(path)
                with open(path, 'w') as f:
                    f.write( pyhtml.unescape( etree.tostring(html.getroottree()).decode('utf-8') ) )
        else:
            try:
                self.svgTree.getroottree().write(path)
            except:
                print("failed saving %s" % path)
        return(path)

    def SVGtoPNG(self,svgPath, pngPath, width=None, inkscapePath="C:\Program Files (x86)\Inkscape\inkscape.exe"):

        if width is not None:
            pass
        else:
            width = self.width


        #cmd = '"%(INKSCAPE_PATH)s" -z --verb=org.ekips.filter.embedimage --verb=FileSave --verb=FileClose -f %(source_svg)s -w %(width)s -j -e %(dest_png)s' %  {'INKSCAPE_PATH':inkscapePath, 'source_svg':svgPath, 'dest_png':pngPath, 'width':width}
        cmd = '"%(INKSCAPE_PATH)s" -z  --verb=FileSave --verb=FileClose -f %(source_svg)s -w %(width)s -e %(dest_png)s' %  {'INKSCAPE_PATH':inkscapePath, 'source_svg':svgPath, 'dest_png':pngPath, 'width':width}
        os.system('%s' % cmd)
        os.system('%s' % cmd)





#circle = bdbplot.getCircle(200,200,50)
#bdbplot.modifyStyle(circle, {'filter':'url(#%s)'%shadow.get('id')})
#bdbplot.svgTree.append( circle )

#circle = bdbplot.getCircle(250,250,50)
#bdbplot.modifyStyle(circle, {'filter':'url(#%s)'%shadow.get('id')})
#bdbplot.svgTree.append( circle )

#path = bdbplot.getPath('M100 100 L300 100 L300 300 L300 100')
#bdbplot.svgTree.append( path )


##
# Spaghettogram
##



##
# Histogram
##

#
# dictionary of read abundace->freq
def readCountHistogram(abundanceFreqDict, logAbundance=True):

    #We expect a distribution which is very steep.
    lookAhead = 3


    if logAbundance:
        f = abundanceFreqDict
        abundanceFreqDict = Counter({})
        for a in f:
            #print(("%s %s" % (a, f[a])))
            try:
                logA = int(round( math.log10(int(a)*100), 0))
            except Exception as e:
                logA = 0
                print(e)
            abundanceFreqDict[logA] += f[a]
           # print(("%s %s -> %s %s" % (a, f[a], logA, abundanceFreqDict[logA] )))

    #Find the highest abundant read:
    hfreq = max([n for n in abundanceFreqDict])

    #Find closed distribution:
    #Mapping from abundance value to plotting X coordinate
    xxMapping = {}

    closedEnd = 1
    perBin = 100
    for c in range(1,perBin):

        if abundanceFreqDict[c]==0 and 0==sum(abundanceFreqDict[q] for q in range(c,c+lookAhead+1)):
            closedEnd = c-1
            break

        else:
            xxMapping[c] = c-0.5

    #check how many extra blocks we need

    extraBlocks = 0
    prevX = closedEnd-0.5
    maxX = 1
    prevAbundance = closedEnd
    extraBlockCoords = []
    extraBlockContinuous = {}
    for abundance in sorted(abundanceFreqDict.keys()):
        if abundance > closedEnd:


            if (abundance-prevAbundance) > 1:
                xxMapping[abundance] = prevX+2
                extraBlockContinuous[extraBlocks] = True
            else:
                xxMapping[abundance] = prevX+1
                extraBlockContinuous[extraBlocks] = False

            prevAbundance=abundance
            extraBlockCoords.append(abundance)
            prevX = xxMapping[abundance]
            maxX= xxMapping[abundance]+1
            extraBlocks+=1
    #print(xxMapping)
    bdbplot = BDBPlot()
    bdbplot.plotStartX = 100
    bdbplot.plotStartY = 150

    bdbplot.plotHeight =400
    bdbplot.plotWidth = max(600, (maxX) * 25)

    bdbplot.setWidth(bdbplot.plotWidth+bdbplot.plotStartX+10)
    bdbplot.setHeight(800)


    bdbplot.xMax = max(1,maxX) # prevent 0 (breaks everything, 0 divisions and such)
    bdbplot.yMax = max(1,int(math.log10(hfreq)))

    axis = bdbplot.getAxis(2)
    bdbplot.svgTree.append( axis )

    #Draw specialised x-axis

    p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, bdbplot.yMin),bdbplot.getPlottingCoord(xxMapping[closedEnd]+1, bdbplot.yMin)]))
    bdbplot.svgTree.append( p )

    for extraBlock in range(0,extraBlocks):
        x = xxMapping[extraBlockCoords[extraBlock]]

        #p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(x, bdbplot.yMin),bdbplot.getPlottingCoord(x+1, bdbplot.yMin)]))

        if  extraBlockContinuous[extraBlock]:
            d = 0.15
            p = bdbplot.getPath(bdbplot.getPathDefinition([
                bdbplot.getPlottingCoord(x-1, bdbplot.yMin),
                bdbplot.getPlottingCoord(x-0.75, bdbplot.yMin+d),
                bdbplot.getPlottingCoord(x-0.25, bdbplot.yMin-d),
                bdbplot.getPlottingCoord(x, bdbplot.yMin)
                ]))

            bdbplot.modifyStyle(p, {'stroke-width':'1', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )


        p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(x, bdbplot.yMin),bdbplot.getPlottingCoord(x+1, bdbplot.yMin)]))
        bdbplot.svgTree.append( p )



        #p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(x+2, bdbplot.yMin),bdbplot.getPlottingCoord(x+3, bdbplot.yMin)]))
        #bdbplot.modifyStyle(p, {'stroke-width':'1', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
        #bdbplot.svgTree.append( p )


    #Draw fine grid
    for y in range(1,bdbplot.yMax+1):
        p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)]))

        bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
        bdbplot.svgTree.append( p )


    ### Block plotting ###
    rectangles = []
    barShadow = bdbplot.shadow(1,1)
    bdbplot.addDef(barShadow)

    for abundance in range(1,int(hfreq)+1):

        if abundanceFreqDict[abundance]>0:
            plotX = xxMapping[abundance]
            frequency = abundanceFreqDict[abundance]

            if frequency==1:
                value=0.20
            else:
                value = math.log10(frequency)

            c = bdbplot.getPlottingCoord(plotX,value)
            origin = bdbplot.getPlottingCoord(plotX,0)

            barWidth = float(bdbplot.plotWidth)/(bdbplot.xMax+1)

            rectangleParams = (c[0], c[1], barWidth,  (float(value)/bdbplot.yMax) * bdbplot.plotHeight-3)
            rectangles.append(rectangleParams)
            bar = bdbplot.getRectangle( *rectangleParams )
            bdbplot.modifyStyle(bar, {'filter':'url(#%s)'%barShadow.get('id'),'fill':'rgba(255,255,255,1)'})
            bdbplot.svgTree.append( bar )

            text = bdbplot.getText(str( bdbplot.humanReadable(frequency,1 ) ), c[0]+0.5*barWidth, c[1]-10,BDBcolor(0,0,0,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Gill Sans MT')
            #text.set('font-family','Cambria Math')
            text.set('font-size', '14')
            text.set('font-weight', 'bold')
            text.set('fill', 'rgba(50,50,50,1)')
            bdbplot.svgTree.append( text )

            #AXIS LABEL
            c = bdbplot.getXLabelCoord(plotX+0.5)
            text = bdbplot.getText(str(bdbplot.humanReadable(abundance)), c[0], c[1]+15,BDBcolor(80,80,80,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.svgTree.append( text )




    for rect in rectangles:
        bar = bdbplot.getRectangle( *rect )
        bdbplot.modifyStyle(bar, {'fill':bdbplot.getGroupColors(1)[0], 'stroke':'#FFFFFF','stroke-width':'1.5'})
        bdbplot.svgTree.append( bar )



    #Y axis label
    for y in range(1,bdbplot.yMax+1):
        c = bdbplot.getYLabelCoord(y)

        value = math.pow(10,y)

        text = bdbplot.getText(str(10), c[0]-10, c[1],BDBcolor(80,80,80,1))
        text.set('text-anchor','end')
        text.set('dominant-baseline','middle')
        text.set('font-family','Cambria Math')
        bdbplot.addSuper(text,str(y))
        bdbplot.svgTree.append( text )



    c = bdbplot.getYLabelCoord( (bdbplot.yMax/2))
    text = bdbplot.getText('Frequency', c[0]-60, c[1]-30,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    #bdbplot.modifyStyle(text, {'font-size': '20'})
    bdbplot.setTextRotation(text,270)
    bdbplot.svgTree.append( text )

    c = bdbplot.getXLabelCoord(bdbplot.xMax/2)
    text = bdbplot.getText('Read abundance', c[0], c[1]+50,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    bdbplot.svgTree.append( text )


    return(bdbplot)




class subdividedHistClass():

    def __init__(self, name, dataPoints, logTransform=False, offset=0):
        self.logTransform = logTransform
        self.totalValue = 0
        self.barSpacerWidth = 5
        self.barWidth = 40
        self.maxValue = 0
        self.bars = []
        self.name = name
        self.startX = offset
        x = offset
        self.width = 0
        for dName,count in dataPoints.most_common():
            self.totalValue+=count
            x+=self.barSpacerWidth
            self.bars.append({'x':x, 'y':count,'name':dName})
            x+=self.barWidth

            self.maxValue = max(self.maxValue, count)

        self.width = x + self.barSpacerWidth - offset

    def plot(self, bdbplot, scarAliases,subClassColors ):
        #Draw full class rectangle:
        value = 1
        if self.logTransform and self.totalValue>0:
            value = math.log(self.totalValue)

        if not self.logTransform:
            value = self.totalValue

        c = bdbplot.getPlottingCoord(self.startX,value)
        origin = bdbplot.getPlottingCoord(self.startX,0)

        rectangleParams = (c[0], c[1], self.width,  (float(value)/bdbplot.yMax) * bdbplot.plotHeight)
        bar = bdbplot.getRectangle( *rectangleParams )
        bdbplot.modifyStyle(bar, {'fill':'rgba(150,150,150,0.8)','stroke-width':'0'})
        bdbplot.svgTree.append( bar )


        text = bdbplot.getText(self.name, c[0]+0.5*self.width, c[1]- 10,BDBcolor(0,0,0,1))
        text.set('text-anchor','middle')
        text.set('dominant-baseline','middle')
        text.set('font-family','Gill Sans MT')
        #text.set('font-family','Cambria Math')
        text.set('font-size', '14')
        text.set('font-weight', 'bold')
        text.set('fill', 'rgba(50,50,50,1)')
        bdbplot.svgTree.append( text )

        barShadow = bdbplot.shadow(1,1)
        bdbplot.addDef(barShadow)
        rectangles= []
        for bar in self.bars:
            #Add bar:
            plotX = bar['x']
            frequency= bar['y']
            className = scarAliases[ bar['name'] ]
            barColor = subClassColors[bar['name']]

            #Add class label to X axis:
            c = bdbplot.getXLabelCoord(plotX+self.barWidth*0.5)
            text = bdbplot.getText(className, c[0], c[1]+15,BDBcolor(80,80,80,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.svgTree.append( text )

            if self.logTransform:
                if frequency==1:
                    value=0.20
                else:
                    value = math.log10(frequency)
            else:
                value = frequency

            c = bdbplot.getPlottingCoord(plotX,value)
            origin = bdbplot.getPlottingCoord(plotX,0)

            #barWidth = float(bdbplot.plotWidth)/(bdbplot.xMax+1)
            barWidth = self.barWidth
            rectangleParams = (c[0], c[1], barWidth,  (float(value)/bdbplot.yMax) * bdbplot.plotHeight)
            rectangles.append(rectangleParams)
            bar = bdbplot.getRectangle( *rectangleParams )
            bdbplot.modifyStyle(bar, {'stroke-width':'0','filter':'url(#%s)'%barShadow.get('id'),'fill':barColor})
            bdbplot.svgTree.append( bar )

            text = bdbplot.getText(str( bdbplot.humanReadable(frequency,3,2 ) ), c[0]+0.5*barWidth, c[1]-10,BDBcolor(0,0,0,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Gill Sans MT')
            #text.set('font-family','Cambria Math')
            text.set('font-size', '14')
            text.set('font-weight', 'bold')
            text.set('fill', 'rgba(50,50,50,1)')

            bdbplot.svgTree.append( text )

            #Percentile:
            text = bdbplot.getText(str( bdbplot.humanReadable(100*(float(frequency)/self.totalValue),3,2 )+'%' ), c[0]+0.5*barWidth, c[1]+10,BDBcolor(255,255,255,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Gill Sans MT')
            #text.set('font-family','Cambria Math')
            text.set('font-size', '14')
            text.set('font-weight', 'bold')
            text.set('fill', 'rgba(255,255,255,1)')

            bdbplot.svgTree.append( text )



def subdividedClassHistogram( classes, logTransform= False, scarAliases={}):

    classSpacerWidth = 20

    currentX = classSpacerWidth
    classIndex = 0
    maxValue = 0

    subHistClassList = []

    for className in classes:
        shc = subdividedHistClass(className,classes[className], logTransform, currentX)
        currentX += shc.width + classSpacerWidth

        if logTransform:
            if shc.totalValue>0:
                maxValue = max(maxValue,math.log(shc.totalValue))
        else:
            maxValue = max(maxValue,shc.totalValue)

        subHistClassList.append( shc )


    bdbplot = BDBPlot()

    ## color list:
    subClassColors = {}
    idx = 0
    #print(len(scarAliases))
    gc = bdbplot.getGroupColors( len(scarAliases) )
    for scar in scarAliases:
        subClassColors[scar] = gc[idx]
        print(('%s -> %s' %(scar, gc[idx])))
        idx+=1

    ## Plot area preparation
    bdbplot.plotStartX = 100
    bdbplot.plotStartY = 100

    bdbplot.plotHeight =400
    bdbplot.plotWidth = currentX

    bdbplot.setWidth(bdbplot.plotWidth+bdbplot.plotStartX+10)
    bdbplot.setHeight(800)
    bdbplot.xMax = max(1,currentX) # prevent 0 (breaks everything, 0 divisions and such)
    if logTransform:
        bdbplot.yMax = max(1,int(maxValue)+1)
    else:
        bdbplot.yMax = max(1,int(maxValue+1))

    axis = bdbplot.getAxis()
    bdbplot.svgTree.append( axis )
    if logTransform:
        for y in range(1,bdbplot.yMax+1):
            p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)], True))

            bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )

    #Draw fine grid
    if logTransform:
        for y in range(1,bdbplot.yMax+1):
            p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)], True))

            bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )

    else:
        stepSize = 50000
        if bdbplot.yMax<101:
            stepSize=10
        for y in range(0,bdbplot.yMax+1,stepSize):
            p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)], True))

            bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )

            c = bdbplot.getYLabelCoord(y)

            text = bdbplot.getText(bdbplot.humanReadable(y,2,3 ), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.svgTree.append( text )


    if logTransform:
        for y in range(1,bdbplot.yMax+1):
            c = bdbplot.getYLabelCoord(y)


            value = math.pow(10,y)

            text = bdbplot.getText(str(10), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.addSuper(text,str(y))
            bdbplot.svgTree.append( text )


    for shc in subHistClassList:
        shc.plot(bdbplot,scarAliases,subClassColors)


    c = bdbplot.getYLabelCoord( (bdbplot.yMax/2))
    text = bdbplot.getText('Reads', c[0]-60, c[1]-30,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    #bdbplot.modifyStyle(text, {'font-size': '20'})
    bdbplot.setTextRotation(text,270)
    bdbplot.svgTree.append( text )

    c = bdbplot.getXLabelCoord(bdbplot.xMax/2)
    text = bdbplot.getText('Sample', c[0], c[1]+50,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    bdbplot.svgTree.append( text )


    return({'plot':bdbplot,'colorMapping':subClassColors})



def classHistogram( classCountMapping, logTransform = False, classColors=None,  placeLeft=None,placeRight=None, reverseOrder=True, height=400, xLabel='Sample', yLabel='Reads', yStepper = 50000, classWidth=50, freqFontSize=14, xLabelFontSize=10,  rotateClassLabels=0,zebraFillMode=False, defaultFillColor='#404040', barSpacing=5, xLabelOffset=50, showZeros=False, axisLabelFontSize=25, drawFreqLabels=True, title=None, freqMethod='humanReadable' ): # freqmethod 'humanReadable' ,'float'

    amountOfClasses = len(classCountMapping)
    #classWidth = 50
    maxValue = 0
    for className in classCountMapping:
        maxValue = max(classCountMapping[className],maxValue)


    bdbplot = BDBPlot()
    bdbplot.plotStartX = 100
    bdbplot.plotStartY = 100

    bdbplot.plotHeight = height
    bdbplot.plotWidth = max(600, (amountOfClasses) * classWidth)

    bdbplot.setWidth(bdbplot.plotWidth+bdbplot.plotStartX+10)
    bdbplot.setHeight(700)

    if title is not None:
        text = bdbplot.getText(title, 10,20, fill='#666666')
        text.set('text-anchor','begin')
        text.set('dominant-baseline','central')
        text.set('font-family','Gill Sans MT')
        text.set('font-size', '25')
        bdbplot.svgTree.append(text)


    bdbplot.xMax = max(1,amountOfClasses) # prevent 0 (breaks everything, 0 divisions and such)

    if logTransform:
        bdbplot.yMax = max(1,int(math.log10(maxValue)+1))
    else:
        bdbplot.yMax = max(1,int(maxValue+1))

    axis = bdbplot.getAxis(2)
    bdbplot.svgTree.append( axis )

    classIndex = 0

    rectangles = []
    barShadow = bdbplot.shadow(1,1)
    bdbplot.addDef(barShadow)

    whiteShadow =  bdbplot.shadow(1,1,3,'rgb(255,255,255)')
    bdbplot.addDef(whiteShadow)

    #Draw fine grid
    if logTransform:
        for y in range(1,bdbplot.yMax+1):
            p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)], True))

            bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )

    else:
        for y in np.arange(0,bdbplot.yMax+yStepper,yStepper):
            p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)], True))

            bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
            bdbplot.svgTree.append( p )

            c = bdbplot.getYLabelCoord(y)
            if(freqMethod=='humanReadable'):
                text = bdbplot.getText(bdbplot.humanReadable(y,2,3 ), c[0]-10, c[1],BDBcolor(80,80,80,1))
            else:
                text = bdbplot.getText(float(y), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')

            bdbplot.svgTree.append( text )



    if isinstance(classCountMapping, collections.OrderedDict):

        if reverseOrder:
            classOrderKeys = list(reversed(list(classCountMapping.keys())))
        else:
            classOrderKeys = list(classCountMapping.keys())

        classOrder = [
            (key, classCountMapping[key]) for key in classOrderKeys
        ]

    else:
        if reverseOrder:
            classOrder = list(reversed(classCountMapping.most_common()))
        else:
            classOrder = list(classCountMapping.most_common())
    classOrderKeys = {}
    for className,freq in classOrder:
        classOrderKeys[className] = (className, freq)
    #Prepend left desired classes to the left
    if placeLeft is not None:
        for className in placeLeft:
            #print(className)
            tup = classOrderKeys[className]
            classOrder.remove(tup)
            classOrder.insert(0, tup)
    if placeRight is not None:
        for className in placeRight:
            #print(className)
            tup = classOrderKeys[className]
            classOrder.remove(tup)
            classOrder.append( tup)

    barWidth = (float(bdbplot.plotWidth)/(bdbplot.xMax+1))
    for className, frequency in classOrder:

        #Add class label to X axis:
        c = bdbplot.getXLabelCoord(classIndex)
        text = bdbplot.getText(str(className), c[0] + ( 0.5*(barWidth)), c[1]+15,BDBcolor(80,80,80,1))
        text.set('text-anchor','middle')
        text.set('dominant-baseline','middle')
        text.set('font-family','Cambria Math')
        text.set('font-size',str(xLabelFontSize))

        if rotateClassLabels!=0:
            bdbplot.setTextRotation(text, rotateClassLabels)
            text.set('text-anchor','start')
        if rotateClassLabels>90:
            bdbplot.setTextRotation(text, rotateClassLabels)
            text.set('text-anchor','end')


        bdbplot.svgTree.append( text )

        #Add bar:
        plotX = classIndex

        if logTransform:


            if frequency==1:
                value=0.20
            elif frequency==0:
                value= 0
            else:
                value = math.log10(frequency)
        else:
            value = frequency

        c = bdbplot.getPlottingCoord(plotX,value)
        origin = bdbplot.getPlottingCoord(plotX,0)




        #barWidth = float(bdbplot.plotWidth)/(bdbplot.xMax+1)


        rectangleParams = (c[0]+0.5*barSpacing, c[1], barWidth-barSpacing,  (float(value)/bdbplot.yMax) * bdbplot.plotHeight-3)
        rectangles.append(rectangleParams)
        bar = bdbplot.getRectangle( *rectangleParams )
        bdbplot.modifyStyle(bar, {'filter':'url(#%s)'%barShadow.get('id'),'fill':'rgba(255,255,255,1)'})
        bdbplot.svgTree.append( bar )

        if (showZeros or frequency>0) and drawFreqLabels:
            if(freqMethod=='humanReadable'):
                text = bdbplot.getText(str( bdbplot.humanReadable(frequency,2,1 ) ), c[0]+0.5*barWidth, c[1]-10,BDBcolor(0,0,0,1))
            else:
                text = bdbplot.getText('%.2f' % frequency  , c[0]+0.5*barWidth, c[1]-10,BDBcolor(0,0,0,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Gill Sans MT')
            #text.set('font-family','Cambria Math')
            text.set('font-size', str(freqFontSize))
            text.set('font-weight', 'bold')
            text.set('fill', 'rgba(50,50,50,1)')
            bdbplot.modifyStyle(text, {'filter':'url(#%s)'%whiteShadow.get('id')})

            bdbplot.svgTree.append( text )

        classIndex += 1


    if logTransform:
        for y in range(1,bdbplot.yMax+1):
            c = bdbplot.getYLabelCoord(y)


            value = math.pow(10,y)

            text = bdbplot.getText(str(10), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.addSuper(text,str(y))
            bdbplot.svgTree.append( text )


    for idx,rect in enumerate(rectangles):
        bar = bdbplot.getRectangle( *rect )

        if zebraFillMode:
            fillColor = bdbplot.getGroupColors(3)[idx%2]
        else:
            fillColor = defaultFillColor
        if classColors is not None:
            if classOrder[idx][0] in classColors:
                fillColor = classColors[classOrder[idx][0]]
            else:
                #print(('Setting %s to default color' % classOrder[idx][0]))
                pass

        bdbplot.modifyStyle(bar, {'fill':fillColor, 'stroke':'#FFFFFF','stroke-width':'1.75'})
        bdbplot.svgTree.append( bar )




    c = bdbplot.getYLabelCoord( (bdbplot.yMax/2))
    text = bdbplot.getText(yLabel, c[0]-60, c[1]-30,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', str(axisLabelFontSize))
    #bdbplot.modifyStyle(text, {'font-size': '20'})
    bdbplot.setTextRotation(text,270)
    bdbplot.svgTree.append( text )

    c = bdbplot.getXLabelCoord(bdbplot.xMax/2)
    text = bdbplot.getText(xLabel, c[0], c[1]+xLabelOffset,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', str(axisLabelFontSize))
    bdbplot.svgTree.append( text )


    return(bdbplot)



class Heatmap(object):

    def __init__(self, npMatrix, colorMatrix=None, rowNames=None,  rowColors=None, columnColors=None, cellFormat=None, cellIdentifiers=None, cellRotations=None, cellStrings=None, columnNames=None, cellSize=25, title=None, subtitle=None, nominalColoring=False, rotateColumnLabels=90, cellAnnot=None, cellAnnotFormat=None, metaDataMatrix=None, cluster=False, groupSize=10 ):

        self.nominalColoring = nominalColoring
        self.nominalColoringMapping = None
        self.zeroColor = None
        self.colormap = matplotlib.cm.get_cmap('plasma')
        self.NanColor = (0.9,0.9,0.9,1)
        self.rotateColumnLabels = rotateColumnLabels
        self.footerHeight = 400
        self.groupSpacerSize = 10
        self.cellFont = 'Cambria'
        self.labelFont = 'Cambria'

        print("Plotting %s by %s matrix" % npMatrix.shape)
        if rowNames is not None:
            print("Supplied %s rownames" % len(rowNames))
            rowNames = [ "%s: %s"%t for t in enumerate(rowNames) ]

        if columnNames is not None:
            print("Supplied %s column names" % len(columnNames))

        if cluster:
            print("Clustering")
            self.matrix = npMatrix
            clusterMatrix = np.zeros( npMatrix.shape )
            for (column,row), value in np.ndenumerate(npMatrix):
                l = self.nominalIndex(value)
                clusterMatrix[column,row] = l

            distances = scipy.spatial.distance.pdist( np.nan_to_num(clusterMatrix.transpose()), 'cityblock' )
            mdistMatrix = scipy.spatial.distance.squareform(distances)
            clustering = scipy.cluster.hierarchy.linkage( mdistMatrix, 'ward' )
            leavesList = list( scipy.cluster.hierarchy.leaves_list(clustering) )
            #npMatrix  = clusterMatrix
            print(leavesList)
        else:
            leavesList = list(range(npMatrix.shape[1]))
        print("Ranging %s" % len(leavesList))

        self.matrix = npMatrix[:,leavesList]
        if colorMatrix is not None:
            self.colorMatrix = colorMatrix[:,leavesList] #Values between zero and one
        else:
            self.colorMatrix = None
        self.rowColors = np.array(rowColors)[leavesList] if rowColors is not None else []
        self.columnColors = columnColors if columnColors is not None else []
        self.cellIdentifiers = cellIdentifiers if cellIdentifiers is not None else None
        self.cellAnnot = cellAnnot
        self.rowNames =np.array(rowNames)[leavesList] if rowNames is not None else []
        self.columnNames = columnNames if columnNames is not None else []
        self.cellSize = cellSize
        self.cellFormat = cellFormat if cellFormat is not None else lambda x: x
        self.cellAnnotFormat  = cellAnnotFormat if cellAnnotFormat is not None else lambda x: x
        self.cellStrings = cellStrings[:,leavesList] if cellStrings is not None else []
        self.cellRotations = cellRotations[:,leavesList] if cellRotations is not None else None

        self.metaDataMatrix = metaDataMatrix[:,leavesList] if metaDataMatrix is not None else None
        self.title = title
        self.subtitle = subtitle
        self.leftMargin = 80
        self.labelWidth = 20
        self.topMargin = 150
        self.cellSpacing = 1
        self.cellFontSize = 10
        self.labelFontSize = 15
        self.groupSize = groupSize
        #self.colormap = matplotlib.cm.get_cmap('inferno')
        print(self.rowNames)


    def getRowName(self, index):
        try:
            return(str(self.rowNames[index]))
        except:
            return('')
    def getColName(self, index):
        try:
            return(str(self.columnNames[index]))
        except:
            return('')

    def getRowColor(self,index):
        try:
            c = self.rowColors[index]

            return( c )
        except:
            return( BDBcolor( 50,50,50, 1) )

    def getColumnColor(self,index):
        try:
            c = self.columnColors[index]
            return( c )
        except:
            return( BDBcolor( 50,50,50, 1) )

    def getCellString(self, row,column):
        try:
            return( self.cellFormat(self.cellStrings[column,row]))
        except:
            return('')

    def getCellAnnotString(self, row,column):
        try:
            return( str(self.cellAnnotFormat(self.cellAnnot[column, row])))
        except:
            return('')

    def getCellId(self, row,column):
        try:
            return( str(self.cellIdentifiers[column, row]) )
        except:
            return(None)
    def getColor(self, value):

        theValueIsNan = self.isnan(value)

        if theValueIsNan:
            r,g,b,a  = self.NanColor
        elif self.zeroColor is not None and value==0:
            r,g,b,a  = self.zeroColor
        else:

            if self.nominalColoring:
                r,g,b,a = self.nominalColor(value)
            else:
                try:
                    r,g,b,a = self.colormap(value)
                except:

                    print("Reverted to nominal coloring mode")
                    self.nominalColoring = True
                    r,g,b,a = self.nominalColor(value)
        return(r,g,b,a)


    def nominalIndex(self,value):
        #Force build of nominal matrix
        r,g,b,a = self.nominalColor(value)
        #get index:
        try:
            idx = list(self.nominalColoringMapping.keys()).index(value)
        except:
            idx=0.01
        return(idx / len(self.nominalColoringMapping.keys()))


    def getCellMetaData(self, row,column):
        try:
            return( str(self.metaDataMatrix[column, row]) )
        except:
            return(None)

    def addTitle(self, plot):
        if self.title is not None:
            text = plot.getText(self.title, 10,20, fill='#666666')
            text.set('text-anchor','begin')
            text.set('dominant-baseline','central')
            text.set('font-family','Gill Sans MT')
            text.set('font-size', '25')
            plot.svgTree.append(text)

    def getCellRotation(self, row, column):
        try:
            if self.cellRotations[column, row] == np.nan:
                return(None)

            return( str(self.cellRotations[column, row]) )
        except:
            return(None)


    def addSubtitle(self, plot):
        if self.subtitle is not None:
            text = plot.getText(self.subtitle, 10,40, fill='#222222')
            text.set('text-anchor','begin')
            text.set('dominant-baseline','central')
            text.set('font-family','Cambria')
            text.set('font-size', '15')
            plot.svgTree.append(text)

    def addGroupTitle(self, plot, title, indexStart, indexEnd):

        matrixGroup = plot.svgTree.findall(".//g[@id='matrix']")[0]
        #Calculate x-starting coordinate
        xStart = self.columnIndexToXCoordinate(indexStart)
        #Calculate x-ending coordinate
        xEnd = self.columnIndexToXCoordinate(indexEnd)+self.cellSize
        yStart = self.rowIndexToYCoordinate(0)-self.cellSize
        yEnd = yStart - self.cellSize*0.5

        p = plot.getPath(plot.getPathDefinition([(xStart, yStart),(xStart,yEnd),(xEnd, yEnd), (xEnd, yStart)]))
        plot.modifyStyle(p, {'stroke-width':'1', 'stroke-linecap':'round','stroke-dashoffset':'0', 'stroke':'#333333'} )
        matrixGroup.append( p )

        text = plot.getText(title, xStart+0.5*( (xEnd-xStart)), yEnd - 8, fill='#000000')
        text.set('text-anchor','middle')
        text.set('dominant-baseline','central')
        text.set('font-family','Gill Sans MT')
        #stext.set('font-family','Cambria')
        matrixGroup.append(text)


    def columnIndexToXCoordinate(self, columnIndex, objSize=None):
        objSize = self.cellSize if objSize is None else objSize
        return(  (self.cellSize + (self.cellSize-objSize)*0.5)  + (columnIndex-1) * self.cellSize+ self.cellSpacing*columnIndex + int(columnIndex/self.groupSize)*self.groupSpacerSize)
        #return(columnIndex * (self.cellSize+self.cellSpacing) + int(columnIndex/4)*self.groupSpacerSize)


    def rowIndexToYCoordinate(self, rowIndex, objSize=None):
        objSize = self.cellSize if objSize is None else objSize
        return( (self.cellSize + (self.cellSize-objSize)*0.5) + (rowIndex-1) * (self.cellSize) + self.cellSpacing*rowIndex + int(rowIndex/self.groupSize)*self.groupSpacerSize)


    def nominalColor(self, value):
        if self.nominalColoringMapping == None:
            #Find all unique values in the matrix:

            uniqueValues = list(set( v for _,v in np.ndenumerate( self.matrix.astype(str ) ) ))
            try:
                uniqueValues = sorted(uniqueValues)
            except:
                pass

            if len(uniqueValues)!=0:
                self.nominalColoringMapping = {val:self.colormap(float(index)/float(len(uniqueValues))) for index,val in enumerate(uniqueValues)}
            else:
                return(self.NanColor)
            print(self.nominalColoringMapping)


        return( self.nominalColoringMapping.get(value, self.NanColor) )

    def isnan(self,value):

        theValueIsNan = False
        if value is None:
            theValueIsNan = True
        else:
            try:
                theValueIsNan = np.isnan(value)
            except:
                theValueIsNan = False
        return(theValueIsNan)

    def getPlot(self):

        plot = BDBPlot()

        matrixGroup = plot.getGroup('matrix')
        matrixGroup.set('transform', 'translate(%s, %s)' % (self.leftMargin,self.topMargin ))
        columnCount, rowCount = self.matrix.shape

        plotWidth = self.leftMargin + (columnCount*(self.cellSize+self.cellSpacing )) + self.groupSpacerSize*(columnCount/4)
        plotHeight = self.topMargin + self.rowIndexToYCoordinate(rowCount) +self.cellSize+self.cellSpacing  + self.footerHeight
        plot.setWidth(plotWidth)
        plot.setHeight(plotHeight)
        ySlack = self.cellSize/2

        cellShadow = plot.shadow(0.5,0.5,1, 'rgb(0,0,0)', 0.98)

        plot.addDef(cellShadow)


        foreGroundTilesGroup = plot.getGroup('foreGroundTiles')
        foreGroundTilesGroup.set('style', 'filter:url(#%s);'%cellShadow.get('id') )
        matrixGroup.append(foreGroundTilesGroup)
        for (column,row), value in np.ndenumerate(self.matrix):

            zIndex = 0
            cellSize = self.cellSize*0.75 if self.getCellRotation(row, column)!=None and self.getCellRotation(row, column)!=0 and self.getCellRotation(row, column)!='nan' else self.cellSize
            x = self.columnIndexToXCoordinate(column,cellSize)
            y = self.rowIndexToYCoordinate(row,cellSize)
            rect = plot.getRectangle(x,y,cellSize, cellSize)

            try:
                cVal = self.colorMatrix[column,row]
            except:
                cVal = value


            theValueIsNan = self.isnan(cVal)
            r,g,b,a = self.getColor(value)

            r = int(r*255.0)
            g = int(g*255.0)
            b = int(b*255.0)
            #print('rgb(%s,%s,%s)' %  (r,g,b))


            rect.set( 'fill','rgb(%s,%s,%s)' %  (r,g,b))
            plot.modifyStyle(rect, {'fill': 'rgba(%s,%s,%s,1)' % (r,g,b), 'stroke-width':'0', 'stroke':'rgba(%s,%s,%s,1)' % (0,0,0)})

            if self.getCellId(row,column) is not None:
                rect.set('cell_id',self.getCellId(row,column))
            if self.getCellMetaData(row,column) is not None:
                rect.set('meta',self.getCellMetaData(row,column))

            if self.getCellRotation(row, column)!=None and self.getCellRotation(row, column)!=0 and self.getCellRotation(row, column)!='nan':
                rect.set('transform','rotate(%s, %s, %s)'%(self.getCellRotation(row, column),x+cellSize*0.5, y+cellSize*0.5))
                zIndex = len(matrixGroup)-1

            if theValueIsNan:
                matrixGroup.insert( 0, rect)
            else:
                foreGroundTilesGroup.insert( zIndex, rect)

            brightness = 255 - ((r+g+b)/3.0)
            if ((r+g+b)/3.0)<100:
                c = BDBcolor( brightness, brightness, brightness, 1)
            else:
                c = BDBcolor( 0, 0, 0, 1)

            if theValueIsNan==False:
                if self.cellAnnot is None:
                    text = plot.getText(str(self.getCellString(row,column)), x+0.5*cellSize, y+0.5*cellSize, fill=c.getRGBStr())
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','central')
                    #text.set('font-family','Gill Sans MT')
                    text.set('font-family',self.cellFont)
                    text.set('font-size', str(self.cellFontSize))
                    matrixGroup.append(text)
                else:
                    text = plot.getText(str(self.getCellString(row,column)), x+0.5*cellSize, y+0.3*cellSize, fill=c.getRGBStr())
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','central')
                    #text.set('font-family','Gill Sans MT')
                    text.set('font-family',self.cellFont)
                    text.set('font-size', str(self.cellFontSize))
                    #text.set('font-weight','bold')
                    matrixGroup.append(text)

                    text = plot.getText(str(self.getCellAnnotString(row,column)), x+0.5*cellSize, y+0.75*cellSize, fill=c.getRGBStr())
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','central')
                    #text.set('font-family','Gill Sans MT')
                    text.set('font-family',self.cellFont)
                    text.set('font-size', str(self.cellFontSize*0.90))
                    matrixGroup.append(text)


            if column==0:
                #x -= self.labelWidth

                self.leftMargin = max( len( self.getRowName(row) )*8, self.leftMargin )

                text = plot.getText(self.getRowName(row), self.columnIndexToXCoordinate(column)-self.labelWidth,  self.rowIndexToYCoordinate(row)+ySlack, fill=self.getRowColor(row))
                text.set('text-anchor','end')
                text.set('dominant-baseline','central')
                #text.set('font-family','Gill Sans MT')
                text.set('font-family',self.labelFont)
                text.set('font-size', str(self.labelFontSize))

                matrixGroup.append(text)
            if row==(rowCount-1) or row==0:
                offset = self.labelWidth+self.cellSize if row==(rowCount-1) else -self.cellSize*0.5
                text = plot.getText(self.getColName(column), self.columnIndexToXCoordinate(column)+self.cellSize*0.5,  self.rowIndexToYCoordinate(row)+offset, fill=self.getColumnColor(column))
                text.set('text-anchor','middle')
                text.set('dominant-baseline','central')
                #text.set('font-family','Gill Sans MT')
                text.set('font-family',self.labelFont)
                text.set('font-size', str(self.labelFontSize))
                if self.rotateColumnLabels!=0:
                    plot.setTextRotation(text, self.rotateColumnLabels)
                    if self.rotateColumnLabels>90:
                        text.set('text-anchor','end' if row==(rowCount-1) else 'start')
                    else:
                        text.set('text-anchor','start' if row==(rowCount-1) else 'end')
                matrixGroup.append(text)

        plot.svgTree.append(matrixGroup)
        self.addTitle(plot)
        self.addSubtitle(plot)

        matrixGroup.set('transform', 'translate(%s, %s)' % (self.leftMargin,self.topMargin ))
        plotWidth = self.leftMargin + (columnCount*(self.cellSize+self.cellSpacing )) + self.groupSpacerSize*(columnCount/4)
        plotHeight = self.topMargin + self.rowIndexToYCoordinate(rowCount) +self.cellSize+self.cellSpacing  + self.footerHeight
        plot.setWidth(plotWidth)
        plot.setHeight(plotHeight)
        return(plot)



def testHeatmap():
    xmin = 0.0
    xmax = 10.0
    dx = 1.0
    ymin=0.0
    ymax=5.0
    dy = 1.0
    x,y = np.meshgrid(np.arange(xmin,xmax,dx),np.arange(ymin,ymax,dy))
    npMat = (x*y)
    m = npMat.max()
    print(m)
    npMat /= m
    print(npMat)
    xlabels = ["X"+str(x) for x in range(npMat.shape[1])]
    ylabels = ["Y"+str(y) for y in range(npMat.shape[0])]
    h = Heatmap(npMat, npMat, rowNames=xlabels, columnNames=ylabels)
    p = h.getPlot()
    p.write('test.svg')


def histogram(values = [1,7,3,2,1,0,0,0,1], rebin=False, binCount=9, reScale=False, logScale=False, logScaleData=False):

    if rebin:
        newBars = {}
        bars = [0]*(binCount+1)
        frequencies = dict()
        for v in values:
            if not v in frequencies:
                frequencies[v]=1
            else:
                frequencies[v]+=1

        #Take log of frequencies
        if logScale:
            for q in frequencies:

                if frequencies[q]>0:
                    frequencies[q] = math.log10(frequencies[q])
                else:
                    frequencies[q] = -1


        minValue = min([float(x) for x in list(frequencies.keys())])
        maxValue = max([float(x) for x in list(frequencies.keys())])
        binSize = (maxValue - minValue)/binCount

        sampleTotal = sum(frequencies.values())

        for binIndex in range(0,binCount+1):
            binStart = minValue+ binIndex*binSize
            binEnd = binStart+binSize
            binTotal = 0
            for d in frequencies.irange(binStart, binEnd, (True,False)):
                binTotal+=frequencies[d]
                if reScale:
                    newBars[binStart+0.5*binSize] =  float(binTotal)/float(sampleTotal)
                    bars[binIndex] = float(binTotal)/float(sampleTotal)
                else:
                    newBars[binStart+0.5*binSize] =  binTotal
                    bars[binIndex] = binTotal


    else:
        bars=values


    bdbplot = BDBPlot()


    bdbplot.plotStartX = 100
    bdbplot.plotStartY = 100
    bdbplot.xMax = max(1,int(math.ceil(len(bars)))) # prevent 0 (breaks everything, 0 divisions and such)
    bdbplot.yMax = max(1,int(math.ceil(max(bars))))


    #plotWall = bdbplot.getRectangle(bdbplot.plotStartX,bdbplot.plotStartY,bdbplot.plotWidth,bdbplot.plotHeight)
    #bdbplot.svgTree.append( plotWall )
    #bdbplot.modifyStyle(plotWall, {'filter':'url(#%s)'%shadow.get('id'),'fill':'rgba(255,255,255,1)'})
    #bdbplot.modifyStyle(plotWall, {'fill':'rgba(255,255,255,1)'})

    axis = bdbplot.getAxis()
    bdbplot.svgTree.append( axis )

    #Draw fine grid
    for y in range(1,bdbplot.yMax+1):
        p = bdbplot.getPath(bdbplot.getPathDefinition([bdbplot.getPlottingCoord(bdbplot.xMin, y),bdbplot.getPlottingCoord(bdbplot.xMax, y)]))

        bdbplot.modifyStyle(p, {'stroke-width':'0.5', 'stroke-linecap':'round', 'stroke-dasharray':'2 2','stroke-dashoffset':'0'} )
        bdbplot.svgTree.append( p )


    # Add axis labels:
    if logScaleData:
        for x in range(0,bdbplot.xMax+1):
            c = bdbplot.getXLabelCoord(x)

            text = bdbplot.getText(str(10), c[0], c[1]+15,BDBcolor(80,80,80,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.addSuper(text,str(x))
            bdbplot.svgTree.append( text )

    else:
        for x in range(0,bdbplot.xMax+1):
            c = bdbplot.getXLabelCoord(x)
            text = bdbplot.getText(str(x), c[0], c[1]+15,BDBcolor(80,80,80,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.svgTree.append( text )
    if logScale:
        for y in range(1,bdbplot.yMax+1):
            c = bdbplot.getYLabelCoord(y)

            value = math.pow(10,y)

            text = bdbplot.getText(str(10), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.addSuper(text,str(y))
            bdbplot.svgTree.append( text )

    else:
        for y in range(0,bdbplot.yMax+1):
            c = bdbplot.getYLabelCoord(y)
            text = bdbplot.getText(str(y), c[0]-10, c[1],BDBcolor(80,80,80,1))
            text.set('text-anchor','end')
            text.set('dominant-baseline','middle')
            text.set('font-family','Cambria Math')
            bdbplot.svgTree.append( text )



    c = bdbplot.getYLabelCoord( (bdbplot.yMax/2))
    text = bdbplot.getText('Frequency', c[0]-60, c[1]-30,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    #bdbplot.modifyStyle(text, {'font-size': '20'})
    bdbplot.setTextRotation(text,270)
    bdbplot.svgTree.append( text )

    c = bdbplot.getXLabelCoord(bdbplot.xMax/2)
    text = bdbplot.getText('Read abundance', c[0], c[1]+50,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    text.set('font-size', '25')
    bdbplot.svgTree.append( text )


    barShadow = bdbplot.shadow(1,1)
    bdbplot.addDef(barShadow)

    for barIndex,barValue in enumerate(bars):
        c = bdbplot.getPlottingCoord(barIndex,barValue)
        origin = bdbplot.getPlottingCoord(barIndex,0)
        bar = bdbplot.getRectangle(  c[0], c[1], float(bdbplot.plotWidth)/(len(bars)+1),  (float(barValue)/bdbplot.yMax) * bdbplot.plotHeight )

        bdbplot.modifyStyle(bar, {'filter':'url(#%s)'%barShadow.get('id'),'fill':'#FFFFFF'})
        bdbplot.svgTree.append( bar )

    for barIndex,barValue in enumerate(bars):
        c = bdbplot.getPlottingCoord(barIndex,barValue)
        origin = bdbplot.getPlottingCoord(barIndex,0)

        barWidth = float(bdbplot.plotWidth)/(len(bars)+1)

        bar = bdbplot.getRectangle(  c[0], c[1], barWidth,  (float(barValue)/bdbplot.yMax) * bdbplot.plotHeight )
        bdbplot.modifyStyle(bar, {'fill':bdbplot.getGroupColors(1)[0], 'stroke':'#FFFFFF','stroke-width':'1.5'})
        bdbplot.svgTree.append( bar )

        if barValue!=-1:
            text = bdbplot.getText(str( bdbplot.humanReadable( int(math.pow(10,barValue)) ) ), c[0]+0.5*barWidth, c[1]-10,BDBcolor(0,0,0,1))
            text.set('text-anchor','middle')
            text.set('dominant-baseline','middle')
            text.set('font-family','Gill Sans MT')
            #text.set('font-family','Cambria Math')
            text.set('font-size', '14')
            text.set('font-weight', 'bold')
            text.set('fill', 'rgba(50,50,50,1)')
            bdbplot.svgTree.append( text )

    return(bdbplot)

def densityXY(scatterData, plotPath, xlabel='x', ylabel='y', logX=False, forceShow=False, logY=False):


    from scipy.stats import gaussian_kde

    if len(scatterData['x'])==0:
        #self.warn('No datapoints left for comparison')
        return(False)

    #@todo: cast this earlier.
    scatterData['x'] = np.array(scatterData['x'])
    scatterData['y'] = np.array(scatterData['y'])

    xy = np.vstack([scatterData['x'],scatterData['y']])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    scatterData['x'] = scatterData['x'][idx]
    scatterData['y'] = scatterData['y'][idx]
    z = z[idx]

    plt.close('all')
    fig, ax = plt.subplots()
    ax.scatter(scatterData['x'], scatterData['y'], c=z, s=50, edgecolor='')
    if logX:
        ax.set_xscale('log')
    if logY:
        ax.set_yscale('log')

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if plotPath and forceShow:
        plt.show()
    if plotPath is None:
        plt.show()
    else:
        plt.savefig(plotPath, bbox_inches='tight')
        plt.close('all')


##
# SIMPLE X Y PLOT
##
def simpleXY():

    bdbplot = BDBPlot()

    shadow = bdbplot.shadow()
    bdbplot.addDef(shadow)
    bdbplot.plotStartX = 100
    bdbplot.plotStartY = 100
    bdbplot.xMax = 10
    bdbplot.yMax = 10

    #plotWall = bdbplot.getRectangle(bdbplot.plotStartX,bdbplot.plotStartY,bdbplot.plotWidth,bdbplot.plotHeight)
    #bdbplot.svgTree.append( plotWall )
    #bdbplot.modifyStyle(plotWall, {'filter':'url(#%s)'%shadow.get('id'),'fill':'rgba(255,255,255,1)'})
    #bdbplot.modifyStyle(plotWall, {'fill':'rgba(255,255,255,1)'})

    axis = bdbplot.getAxis()
    bdbplot.svgTree.append( axis )

    # Add axis labels:
    for x in range(0,11):
        c = bdbplot.getXLabelCoord(x)
        text = bdbplot.getText(str(x), c[0], c[1]+15,BDBcolor(80,80,80,1))
        text.set('text-anchor','middle')
        text.set('dominant-baseline','middle')
        text.set('font-family','Cambria Math')
        bdbplot.addSuper(text,'y')
        bdbplot.svgTree.append( text )

    for y in range(1,11):
        c = bdbplot.getYLabelCoord(y)
        text = bdbplot.getText(str(y), c[0]-15, c[1],BDBcolor(80,80,80,1))
        text.set('text-anchor','middle')
        text.set('dominant-baseline','middle')
        text.set('font-family','Cambria Math')
        bdbplot.svgTree.append( text )



    c = bdbplot.getYLabelCoord(5)
    text = bdbplot.getText('Y Axis', c[0]-40, c[1],BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    bdbplot.setTextRotation(text,270)
    bdbplot.svgTree.append( text )

    c = bdbplot.getXLabelCoord(5)
    text = bdbplot.getText('X Axis', c[0], c[1]+40,BDBcolor(0,0,0,1))
    text.set('text-anchor','middle')
    text.set('dominant-baseline','middle')
    text.set('font-family','Gill Sans MT')
    bdbplot.svgTree.append( text )




    for x in range(0,10):

        c = bdbplot.getPlottingCoord(x,x)
        circle = bdbplot.getCircle(c[0],c[1],2)
        #bdbplot.modifyStyle(circle, {'filter':'url(#%s)'%shadow.get('id')})
        bdbplot.svgTree.append( circle )

    #
    #for i in range(0,10):
    #
    #   a = ((math.pi*2)/10.0) * i
    #   text = bdbplot.getText('%s' % i,250 + 50*math.cos(a),250 + 50*math.sin(a),BDBcolor(80,80,80,1))
    #   text.set('text-anchor','middle')
    #   text.set('dominant-baseline','middle')
    #   text.set('font-family','Cambria Math')
    #   bdbplot.svgTree.append( text )


    bdbplot.dump()
    bdbplot.write('test.svg')


#vals = [0,11]
#
#import random
#for i in range(0,16):
#   vals += [i]* int(math.ceil(math.exp(i/2+1)))
#
#
##print('c(%s)' % ','.join(str(i) for i in vals))
#
#plot = histogram(vals, True,15, False, True)
#
#
#plot.write('test.svg')

#
#d = Counter({1:100,2:50,3:10,5:5,6:10,  1000:1, 500:2, 10000:3})
#plot = readCountHistogram(d)
#text = plot.getText('Embryo 1, component 1',10,40)
#text.set('font-family','Gill Sans MT')
#text.set('font-size', '42')
#plot.svgTree.append( text )
#
#totalCount = sum( [v*d[v]for v in d] )
#text = plot.getText( '%s reads total' % plot.humanReadable(totalCount),10,75, BDBcolor(77,77,77,1))
#text.set('font-family','Gill Sans MT')
#text.set('font-size', '23')
#plot.svgTree.append( text )
#
#
#plot.write('test.svg')

import networkx as nx
import scipy.interpolate
from Bio import pairwise2


class GraphRenderer():

    def interpolate(self, interpolateValue,  colorScaleKeys, nodeColorMapping):

        #Seek positions around value to interpolate
        first = colorScaleKeys[0]
        index = 0
        last = first
        for value in colorScaleKeys:

            if value>=interpolateValue:
                last = value
                break
            else:
                first = value
            index+=1
        if value==interpolateValue:
            return(nodeColorMapping[value])

        #Do interpolation
        colorA = nodeColorMapping[first]
        colorB = nodeColorMapping[last]
        dx = last-first

        return( self._ipol(colorA[0], colorB[0], first, last, interpolateValue), self._ipol(colorA[1], colorB[1], first, last, interpolateValue), self._ipol(colorA[2], colorB[2], first, last, interpolateValue))


    def _ipol(self,a, b, first, last,  interpolateValue):
        #Due to floating point rounding errors the interpolate value can be very close to last,
        # it is ok to return last in those cases
        if last>first and interpolateValue>=last:
            return(b)
        if last<first and interpolateValue>=first:
            return(a)

        y_interp = scipy.interpolate.interp1d([first, last], [a,b])
        return( y_interp(interpolateValue) )


    def sortByIndexAndBase(self, value):

        parts = value.split('_')
        pos = ['A','T','C','G','N'].index(parts[1])
        if pos==None:
            pos=0
        else:
            pos+=1
        return( int(parts[0]) + pos*0.1 )


    def __init__(self, nxGraph, coloringMode = 'nodeRGB', coloringAttribute='confidence', performDistanceMeasure=True, performFrequencyMeasure=True, alias='none'):
        self.g = nxGraph
        self.undirectedG = self.g.to_undirected()
        self.plot = BDBPlot()

        self.nodeShadow = self.plot.shadow(0.5,0.5,1, 'rgb(0,0,0)', 0.98)
        self.plot.addDef(self.nodeShadow)

        minX = 0
        maxX = 0
        minY = 0
        maxY = 0
        for nodeName in self.g:
            node = self.g.node[nodeName]
            if 'x' in node and 'y' in node and 'size' in node :
                if (node['x']-node['size'])<minX:
                    minX = node['x']-node['size']
                if (node['x']+node['size'])>maxX:
                    maxX = node['x']+node['size']

                if (node['y']-node['size'])<minY:
                    minY = node['y']-node['size']
                if (node['y']+node['size'])>maxY:
                    maxY = node['y']+node['size']



        #Estimate color scale: #############
        colorScaleY = 0
        colorScaleHeight = 0
        colorScaleGraphSpacing = 0
        createdColorScale = False
        if coloringAttribute is not None and coloringMode == 'nodeRGB':
            createdColorScale= True
            colorScaleParts = 10
            colorScaleWidth = 400
            colorScaleHeight = 35
            colorScaleSpacing = 5
            colorScaleShadow = self.plot.shadow(1,1,1)
            colorScaleDeltaX = float(colorScaleWidth-colorScaleSpacing)/colorScaleParts
            self.plot.addDef(colorScaleShadow)
            labelHeight = 15
            colorScaleX = 5
            colorScaleGraphSpacing = 10
            colorScaleY = (maxY-minY)+colorScaleGraphSpacing

            nodeColorMapping = {}
            abundanceMapping = Counter({})
            lowestValue = 100000
            highestValue = -lowestValue
            for nodeName in self.g:
                node = self.g.node[nodeName]
                if coloringAttribute in node and 'r' in node:
                    if node[coloringAttribute]==1:
                        value = 0
                    else:
                        value = -math.log( 1.0 - node[coloringAttribute],10 )

                    abundanceMapping[value]+= node['abundance']
                    lowestValue = min(lowestValue, value)
                    highestValue = max(highestValue, value)
                    nodeColorMapping[value] = (node['r'], node['g'], node['b'])

            #Create color scale; ##########
            lowestValue = 2.5
            highestValue = 3.4

            #lowestValue = 3.0
            #highestValue = 3.6
            print((lowestValue, highestValue))
            colorScale = self.plot.getGroup('colorScale')

            self.plot.svgTree.append(colorScale)
            r = self.plot.getRectangle(colorScaleX,colorScaleY,colorScaleWidth,colorScaleHeight+labelHeight+colorScaleSpacing )
            colorScale.append(r)

            r.set('style', 'fill:#FFFFFF' )
            #r.set('style', 'filter:url(#%s);fill:#FFFFFF'%colorScaleShadow.get('id') )
            colorScaleKeys = sorted(OrderedDict(sorted(nodeColorMapping.items())) )
            #print(nodeColorMapping)
            deltaValue = (highestValue-lowestValue) / (colorScaleParts-1)
            currentValue = lowestValue
            x = colorScaleSpacing+colorScaleX
            y = colorScaleY + colorScaleSpacing

            idx = 0
            while idx<colorScaleParts:
                print(('Interpolating for %s' % currentValue ))
                c = self.interpolate(currentValue, colorScaleKeys, nodeColorMapping)
                r = self.plot.getRectangle(x,y,colorScaleDeltaX-colorScaleSpacing,colorScaleHeight-2*colorScaleSpacing )
                r.set('fill',  'rgb(%s, %s, %s)' % c)
                r.set('style', 'filter:url(#%s);'%  self.nodeShadow.get('id') )
                r.set('stroke',  'None')
                colorScale.append( r )

                text = self.plot.getText('%.1f' %  (10.0*currentValue) ,x+0.5*(colorScaleDeltaX-colorScaleSpacing),y+colorScaleHeight+colorScaleSpacing, BDBcolor(0,0,0,1))
                #self.plot.setTextRotation(text, 90)
                text.set('text-anchor','middle')
                text.set('dominant-baseline','central')
                text.set('font-family','Gill Sans MT')
                #text.set('font-family','Cambria Math')
                text.set('font-size', '10')
                #text.set('font-weight', 'bold')
                text.set('fill', 'rgba(0,0,0,1)')
                self.plot.svgTree.append( text )


                currentValue+=deltaValue
                x += colorScaleDeltaX
                idx+=1


        #######################


        #Colorize nodes by distance to center


        distancesFound = Counter({})
        distancesReads = Counter({})
        for node in self.g:
            if 'idStr' in self.g.node[node] and  'Wt' == self.g.node[node]['idStr']:
                centerNode = node
                break

        ldistThreshold = 4
        maxHamming = 8
        classColors = {'H0':'#0000DD','H1': '#66A43E', 'H2': '#0D40DB', 'H3': '#3970DD', 'H4': '#769AE0', 'H5': '#A6B9DD', 'H6': '#D7DAE0', 'H7': '#FFFFFF', 'N1':'#FFCC00','N2':'#FF6600', 'N3':'#C83737','N4':'#800000'}
        if performDistanceMeasure:

            print('Estimating all distances...')
            weirdSequences = []
            for nodeIndex,node in enumerate(self.undirectedG):

                if nodeIndex%100==0:
                    completion = 100.0*(float(nodeIndex)/len(self.undirectedG))
                print('\rcompletion %s    ' % completion, end=' ')

                abundance = self.g.node[node]['abundance']
                hammingDistance = bdbbio.getHammingDistance(node, centerNode)
                unformattedAlignments = pairwise2.align.localxx(node,centerNode) #bdbbio.getLevenshteinDistance(node,centerNode)
                ldist = len(node)

                self.g.node[node]['exactHammingDistance'] =  hammingDistance
                self.g.node[node]['exactNWDistance'] =  ldist
                if len(unformattedAlignments)>0:
                    ldist = int(round(len(node)-float(unformattedAlignments[0][2])))
                    #pairwise2.format_alignment(*unformattedAlignments[0])

                if hammingDistance<maxHamming and ldist<=hammingDistance:
                    distancesReads['H%s'%hammingDistance]+=abundance

                    if hammingDistance==0:
                        self.g.node[node]['color'] = '#66A43E'
                        distancesFound['H0']+=1
                    elif hammingDistance==1:
                        self.g.node[node]['color'] = '#0D40DB'
                        distancesFound['H1']+=1
                    elif hammingDistance==2:
                        self.g.node[node]['color'] = '#2362E0'
                        distancesFound['H2']+=1
                    elif hammingDistance==3:
                        self.g.node[node]['color'] = '#769AE0'
                        distancesFound['H3']+=1
                    elif hammingDistance==4:
                        self.g.node[node]['color'] = '#A6B9DD'
                        distancesFound['H4']+=1
                    elif hammingDistance==5:
                        self.g.node[node]['color'] = '#D7DAE0'
                        distancesFound['H5']+=1
                    elif hammingDistance==6:
                        self.g.node[node]['color'] = '#FFFFFF'
                        distancesFound['H6']+=1
                else:


                    if ldist<=ldistThreshold:

                        distancesFound['N%s'%ldist]+=1
                        distancesReads['N%s'%ldist]+=abundance
                        if 'N%s'%ldist not in classColors:

                            brightness = 255-int(round((float(ldist)/ldistThreshold)*100))
                            self.g.node[node]['color'] = 'rgb(%s,%s,%s)' % (brightness,0,0)
                            classColors['N%s'%ldist] = self.g.node[node]['color']
                        self.g.node[node]['color'] = classColors['N%s'%ldist]
                    else:
                        self.g.node[node]['color'] = '#404040'
                        #distancesFound['l%s'%ldist]+=1
                        #distancesReads['l%s'%ldist]+=abundance
                        distancesFound['N>%s'%ldistThreshold]+=1
                        distancesReads['N>%s'%ldistThreshold]+=abundance
                        weirdSequences.append(SeqIO.SeqRecord(Seq(node), 'NW%s-a%s-%s' % (str(ldist),abundance,str(nodeIndex))))

            #classHistogram( classCountMapping, logTransform = False, classColors=None, placeLeft=None ):

            classHistogram(distancesFound, True, classColors, None,['N>%s'%ldistThreshold], False).write('%s_distancesByNodes.svg' % alias)
            classHistogram(distancesReads, True, classColors, None,['N>%s'%ldistThreshold], False).write('%s_distancesByCountB.svg' % alias)
            nx.write_graphml( self.g, './%s-distanceAnnotated.graphml' % alias)
            fastaPath = './%s-weirdSequences.fa' % alias
            SeqIO.write(weirdSequences, fastaPath, "fasta")




        if performFrequencyMeasure:

            sequenceColors = {'rest':'#404040'}
            sequenceFrequencies = Counter({})
            frequencyTable = []
            frequencyCounter = Counter({})
            for nodeIndex,node in enumerate(self.undirectedG):
                abundance = self.g.node[node]['abundance']
                nodeName = nodeIndex
                if 'idStr' in self.g.node[node]:
                    nodeName = self.g.node[node]['idStr']
                #if 'r' in self.g.node[node]:
                #   sequenceColors[str(nodeName)] = 'rgb(%s, %s, %s)' % (self.g.node[node]['r'],self.g.node[node]['g'],self.g.node[node]['b'])
                if 'color' in self.g.node[node]:
                    sequenceColors[str(nodeName)] = self.g.node[node]['color']

                if abundance>3:
                    sequenceFrequencies[str(nodeName)] += abundance

                else:
                    #sequenceFrequencies['rest'] += abundance
                    pass
                for index,base in enumerate(str(node)):

                    if index>=len(frequencyTable):
                        frequencyTable.append(Counter({}))
                        for b in ['A','T','C','G']:
                            frequencyCounter['%s_%s' % (index, b)] = 0
                    frequencyCounter['%s_%s' % (index, base)]+=abundance
                    frequencyTable[index][base]+= abundance


            freqKeys = sorted(list(frequencyCounter.keys()), key=lambda x: self.sortByIndexAndBase(x))
            colors = {}
            for key in freqKeys:
                parts = key.split('_')
                i = ['N','A','T','C','G'].index(parts[1])
                if i==None:
                    i = 0
                colors[key] = ['#404040', '#336bbd','#ff6600','#aa0000','#5aa02c'][i]


            classHistogram(sequenceFrequencies, True, sequenceColors, None,[], False, classWidth=40, rotateClassLabels=90).write('%s_sequenceFrequenciesLog.svg' % alias)
            classHistogram(sequenceFrequencies, False, sequenceColors, None,[], False, classWidth=40, rotateClassLabels=90).write('%s_sequenceFrequenciesLin.svg' % alias)
            classHistogram(frequencyCounter, True, colors, None,freqKeys, False, classWidth=30, rotateClassLabels=90).write('%s_baseFrequencies.svg' % alias)



        print(('Offsetting %s %s' % (minX, minY)))
        self.plot.setWidth( maxX-minX )
        if createdColorScale:
            self.plot.setHeight( ( (colorScaleY+colorScaleHeight+colorScaleGraphSpacing)-minY ))
        else:
            self.plot.setHeight( maxY - minY + 10)
        h = maxY-minY
        ## plotting
        edgeGroup = self.plot.getGroup('edges')
        nodeGroup = self.plot.getGroup('nodes')
        labelGroup = self.plot.getGroup('labels')

        print('Adding edges')

        for fromNode,toNode,data in self.g.edges(data=True):

            if 'hdist' in data and data['hdist']==1:
                fn = self.g.node[fromNode]
                tn = self.g.node[toNode]
                p = self.plot.getPath( self.plot.getPathDefinition([ (int(round(fn['x']-minX)), h-int(round(fn['y']-minY))), (int(round(tn['x']-minX)), h-int(round(tn['y']-minY))) ]) )

                etree.strip_attributes(p,'style')
                if coloringMode == 'nodeRGB':
                    if 'r' in tn:
                        p.set('stroke', 'rgb(%s, %s, %s)' % (tn['r'],tn['g'],tn['b']))
                else:
                    if 'color' in tn:
                        p.set('stroke', '%s' % (tn['color']))
                p.set('stroke-width', '0.75')
                p.set('stroke-opacity', '0.8')
                #self.plot.modifyStyle(p, { 'stroke-width':'0.5', 'stroke-opacity':"0.5"  } ) #,'stroke-width':'0.5' 'stroke':'rgba(80,80,80,0.8)'})
                #self.plot.modifyStyle(p, { 'stroke':'rgba(80,80,80,0.8)' } ) #,'stroke-width':'0.5' 'stroke':'rgba(80,80,80,0.8)'})

                edgeGroup.append(p)

        #self.plot.modifyStyle(edgeGroup, {'stroke-width':'0.5', 'stroke':'rgba(80,80,80,0.8)'})
        self.plot.svgTree.append(edgeGroup)
        self.plot.svgTree.append(nodeGroup)
        self.plot.svgTree.append(labelGroup)
        print('Adding nodes')





        for componentIndex,connectedComponent in enumerate(nx.connected_component_subgraphs(self.undirectedG)):
            componentGroup = self.plot.getGroup('component_%s' % componentIndex)
            smallNodes = self.plot.getGroup('component_%s_smallNodes' % componentIndex)
            bigNodes = self.plot.getGroup('component_%s_bigNodes' % componentIndex)
            componentGroup.append(smallNodes)
            componentGroup.append(bigNodes)
            nodeGroup.append(componentGroup)
            for nodeName in connectedComponent:
                node = self.g.node[nodeName]
                if 'x' in node and 'y' in node and 'size' in node :


                    if node['size']>0:
                        circle = self.plot.getCircle(int(round(node['x']-minX)), h-int(round(node['y']-minY)), int(round(node['size'])))
                        circle.set('style','fill:none')
                        if coloringMode == 'nodeRGB' and 'r' in node:
                            self.plot.modifyStyle(circle, {'fill':'rgb(%s,%s,%s)' % (node['r'], node['g'], node['b']), 'stroke':'rgba(0,0,247,0.8)'})
                        else:
                            if 'color' in node:
                                self.plot.modifyStyle(circle, {'fill':'%s' % (node['color']), 'stroke':'rgba(0,0,247,0.8)'})

                        if node['size']>0.0:
                            bigNodes.append(circle)
                            bigNodes.append(circle)
                            if node['abundance']>500:
                                if 'idStr' in node:
                                    label = node['idStr']
                                else:
                                    label= ''

                                text = self.plot.getText('%s %s' % (label,node['abundance']) ,int(round(node['x']-minX)), h-int(round(node['y']-minY))+4,BDBcolor(80,80,80,1))
                                text.set('text-anchor','middle')
                                text.set('dominant-baseline','central')
                                text.set('font-family','Gill Sans MT')
                                #text.set('font-family','Cambria Math')
                                text.set('font-size', '14')
                                #text.set('font-weight', 'bold')
                                text.set('fill', 'rgba(50,50,50,1)')
                                labelGroup.append( text )


                        else:
                            smallNodes.append(circle)
                            smallNodes.append(circle)

                else:
                    print('Skipped a node; could not find coordinates')

            bigNodes.set('style', 'filter:url(#%s);'%self.nodeShadow.get('id') )
            #self.plot.modifyStyle(bigNodes, {'filter':'url(#%s)'%self.nodeShadow.get('id')})








def testGraphRenderer():
    p = GraphRenderer( nx.read_graphml("C:\\Users\BuysDB\Desktop\Control1-g3.graphml"),'distances',alias='controlDist')
    p.plot.write('control_spaghettogramRenderDistancesB.svg')
    p.plot.SVGtoPNG('control_spaghettogramRenderDistancesB.svg', 'control_spaghettogramRenderDistancesB.png',2048)

    p = GraphRenderer( nx.read_graphml("C:\\Users\BuysDB\Desktop\Control1-g3.graphml"),'nodeRGB',alias='controlConf')
    p.plot.write('control_spaghettogramRenderConfidenceB.svg')
    p.plot.SVGtoPNG('control_spaghettogramRenderConfidenceB.svg', 'control_spaghettogramRenderConfidenceB.png',2048)


def artGraphRenderer():
    p = GraphRenderer( nx.read_graphml("C:\\Users\BuysDB\Desktop\ArtSimulatedGraphHC26k.graphml"),'nodeRGB','confidence',False, True, 'simulated')
    p.plot.write('ArtSimulatedGraphHC26k.svg')
    p.plot.SVGtoPNG('ArtSimulatedGraphHC26k.svg', 'ArtSimulatedGraphHC26k.png',2048)

def embryoGraphRenderer():
    p = GraphRenderer( nx.read_graphml("C:\\Users\BuysDB\Desktop\embryo7remapped2.graphml"),'nodeRGB','confidence',False, True, 'embryo7')
    p.plot.write('embryo7remapped3.svg')
    p.plot.SVGtoPNG('embryo7remapped3.svg', 'embryo7remapped3.png',2048)

def midbrainGraphRenderer():
    p = GraphRenderer( nx.read_graphml("C:\\Users\BuysDB\Desktop\midbrain.graphml"),'nodeRGB','confidence',False, True, 'brain')
    p.plot.write('midbrain.svg')
    p.plot.SVGtoPNG('midbrain.svg', 'midbrain.png',2048)

class SequenceBin():

    def __init__(self, sequence, abundance):
        self.sequence = sequence
        self.abundance = abundance
        self.confidences = []
        self.diffIndices = []
        self.x = 0
        self.y = 0




class HammingBin():

    def __init__(self, index, hammingDistance=0):
        self.hammingDistance = hammingDistance
        self.index = index
        self.sequences = []

    def addSequence(self, sequence, abundance):
        self.sequences.append(sequence)



#SVG table renderer
class SVGTable(BDBPlot):
    #@param datamatrix list of lists containing values
    #@param header list containing column names
    def __init__(self, dataMatrix, header):
        BDBPlot.__init__(self)
        self.data = dataMatrix
        self.header = header
        self.cellPointer = 0

    #Render a cell in the matrix
    def cellRenderFunction(self,x,y):
        self.svgTree.getGroup()



class SpaghettoPlot():

    def __init__(self, networkxGraph):
        self.g = networkxGraph
        self.minRadius = 1
        self.nucleotideColours = {'A':'#FF2222','T':'#22FF22', 'G':'#2222FF','C':'#FFFF22'}


    def getSurfaceBasedNodeRadius(self, abundance, maxRadius):

        #return(math.log(abundance+1)+1)
        maxO = math.pi*math.pow(float(maxRadius),2)
        return( max(self.minRadius,math.sqrt( float(abundance*maxO)/math.pi) ) )



    def layout(self):


        maxRadius = 100
        minDistance = 15 #Distance between the nodes

        components = []
        self.readLen = 0

        toRemove = []
        for nodeA,nodeB,d in self.g.edges_iter(data='hdist'):
            if d!=1:
                toRemove.append( (nodeA, nodeB) )
        self.g.remove_edges_from(toRemove)

        for componentIndex,connectedComponent in enumerate(nx.connected_component_subgraphs(self.g)):
            if len(connectedComponent)>3:
                #Find center(s)
                plot = BDBPlot()

                centerNode = None
                centerAbundance = 0
                for node in connectedComponent:
                    a = connectedComponent.node[node]['abundance']
                    if a > centerAbundance:
                        centerNode = node
                        centerAbundance = a
                    #Todo: compat for more centers
                    self.readLen = max(self.readLen, len(node))
                #Estimate radial size of connected component:
                #longest path from center to member
                radialSize = 1

                distanceMap = {0:[centerNode]}
                for targetNode in connectedComponent:
                    if targetNode!=centerNode:

                        pathLen = nx.shortest_path_length(self.g,source=centerNode,target=targetNode)
                        radialSize = max( radialSize, pathLen)
                        if not pathLen in distanceMap:
                            distanceMap[pathLen] = []
                        distanceMap[pathLen].append(targetNode)


                #print(distanceMap)

                #Read index to angle mapping
                phis = []
                for index in range(0,self.readLen):
                    phis.append( -math.pi*0.5 + float(math.pi*2) * (float(index)/self.readLen) )



                currentRadius = self.getSurfaceBasedNodeRadius(float(self.g.node[centerNode]['abundance'])/centerAbundance,maxRadius)
                centerRadius = currentRadius
                # Construct coordinates
                coordinates = {}
                coordinates[centerNode] = {'x':0,'y':0,'r':currentRadius}

                hammingRadials = []
                for distance in range(1,radialSize+1):

                    print(('Distance %s radius: %s' % (distance, currentRadius)))
                    if distance in distanceMap:
                        #Find the radius of this circle
                        maxNodeRadius = 0
                        for node in distanceMap[distance]:
                            maxNodeRadius = max(maxNodeRadius, self.getSurfaceBasedNodeRadius(float(self.g.node[node]['abundance'])/centerAbundance,maxRadius))

                        currentRadius += (maxNodeRadius+ minDistance) #added *.5!

                        #Calculate x and y coordinates:

                        for node in distanceMap[distance]:
                            r = self.getSurfaceBasedNodeRadius(float(self.g.node[node]['abundance'])/centerAbundance,maxRadius)

                            #Retrieve the hamming index of the node
                            indices = bdbbio.getHammingIndices(node,centerNode) #[ int(x) for x in self.g.node[node]['indices'].split(',') ]

                            #Take the first index as base
                            for index in indices:
                                angle = phis[ index ]

                                coordinates[node] = {'x':math.cos(angle)*currentRadius, 'y':math.sin(angle)*currentRadius, 'r':r, 'a':angle }
                            #print(coordinates[node])
                    hammingRadials.append(currentRadius)

                translateX = -currentRadius
                translateY = -currentRadius
                #for nodeName in coordinates:
                #   translateX = min(translateX, coordinates[nodeName]['x']-coordinates[nodeName]['r'] )
                #   translateY = min(translateY, coordinates[nodeName]['y']-coordinates[nodeName]['r'] )
                #Draw component

                for i,r in enumerate(hammingRadials):

                    circle = plot.getCircle( -translateX,  -translateY, r)
                    plot.modifyStyle(circle, {'stroke-width':'0.3', 'stroke-miterlimit':'4', 'stroke-dasharray':'2.08,2.08'})
                    plot.svgTree.append( circle )

                    text = plot.getText('h'+str(i+1),-translateX, -r-translateY+5, fill=BDBcolor(30,30,30,0))
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','central')
                    text.set('font-family','Cambria')
                    text.set('font-size','6')
                    plot.svgTree.append(text)


                    if i==0:
                        r=centerRadius*2 + 10
                        for index, angle in enumerate(phis):
                            text = plot.getText(centerNode[index], math.cos(angle)*r*0.5 - translateX, math.sin(angle)*r*0.5-translateY, fill=BDBcolor(30,30,30,0))
                            text.set('text-anchor','middle')
                            text.set('dominant-baseline','central')
                            text.set('font-family','Gill Sans MT')
                            text.set('font-size','10')
                            plot.svgTree.append(text)

                            p = r - 35
                            text = plot.getText(str(index), math.cos(angle)*p*0.5 - translateX, math.sin(angle)*p*0.5-translateY, fill=BDBcolor(30,30,30,0))
                            text.set('text-anchor','middle')
                            text.set('dominant-baseline','central')
                            text.set('font-family','Cambria')
                            text.set('font-size','8')
                            plot.svgTree.append(text)



                for nodeName in coordinates:

                    circle = plot.getCircle( coordinates[nodeName]['x']-translateX,  coordinates[nodeName]['y']-translateY, coordinates[nodeName]['r'])
                    plot.modifyStyle(circle, {'stroke-width':'0', 'fill':'#FF6655', 'fill-opacity':'0.30'})
                    plot.svgTree.append( circle )

                plot.write('./components/component%s.svg' % componentIndex )
