#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from singlecellmultiomics.utils import bdbplot,organoidTools
from copy import deepcopy
from importlib import reload
import numpy as np

class StatePlotter():

    def __init__(self, plot=None):


        self.canvas = bdbplot.BDBPlot() if plot is None else plot
        self.heightPerState = 10
        self.stateMargin = 1
        self.chromosomeMargin = 2
        self.pixelsPerBase = 1/10_000_000
        self.headerHeight=30
        self.bigFontSize = 10
        self.smallFontSize = 6
        self.gray = 220

        self.gainColor = '#B36466'
        self.gainTwoColor = '#FF0087'
        self.normalColor = '#FFFFFF'
        self.lossColor =  '#6562B5'
        self.totalLossColor  = '#004BD8'
        self.missingColor = '#AAAAAA'

    def plotStates(self, df, offset=(0,0), **kwargs):
        self.offset = offset
        stateGroup = self.canvas.getGroup('stateGrid')
        self.stateGroup = stateGroup
        stateGroup.set('transform' ,f"translate({offset[0]},{offset[1]})" )
        self.canvas.svgTree.append(stateGroup)
        self.states = sorted(list(set(df['cluster'])))
        self.chromosomeOrder =  sorted(sorted(list(set(df[ 'chromosome']))), key=organoidTools.chrom_sort_human)

        x = 10

        for i,state in enumerate(self.states):
            g = self.canvas.getGroup(f'state_{state}')
            stateGroup.append(g)

            result = self.plotState(self.canvas, g, df[ df['cluster']==state ], x,
                           y= (self.heightPerState*i + (i*self.stateMargin)), label=state ,
                                   plotChromosomeLabels=(i==0), **kwargs
                                   )

        return self.canvas


    def plotState(self, plot, g, row, x, y, label, plotChromosomeLabels=False,
                  logScale=False, chromosomeSizes=None, logRepeats=False):

        currentX = x


        shownAlleles = False
        for chromosome in self.chromosomeOrder:
            print(chromosome)
            isAllelic = chromosome.endswith('_A') or chromosome.endswith('_B')
            #Obtain how many bases this chromosome has
            chrom = chromosome.split('_')[0]
            if chromosomeSizes is not None and chrom in chromosomeSizes:
                chromosomeSize = chromosomeSizes[chrom]
            else:

                chromosomeSize = row[ row['chromosome']==chromosome ]['endCoordinate'].max()
                

            chromosomePixelWidth = chromosomeSize * self.pixelsPerBase

            rect =  self.canvas.getRectangle(currentX, y, chromosomePixelWidth, self.heightPerState)
            self.canvas.modifyStyle( rect, {'fill':f'rgb({self.gray},{self.gray},{self.gray})', 'stroke':'none'})
            g.append(rect)

            if plotChromosomeLabels:

                if  not chromosome.endswith('_B'): # allelic
                    text = plot.getText(chromosome.replace('chr','').replace('_A',''), currentX+chromosomePixelWidth*0.5 + (chromosomePixelWidth*0.5 if chromosome.endswith('_A') else 0), y - self.heightPerState)
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','middle')
                    text.set('font-family','Helvetica')
                    text.set('font-size', str(self.smallFontSize))

                    g.append(text)

                if isAllelic:
                    alleleDescriptor = 'allele ' if not shownAlleles else ''
                    allele = 'A' if chromosome.endswith('_A') else 'B'
                    text = self.canvas.getText( f'{alleleDescriptor}{allele}', currentX+chromosomePixelWidth*0.5 , y - 0.4*self.heightPerState)
                    text.set('text-anchor','middle')
                    text.set('dominant-baseline','middle')
                    text.set('font-family','Helvetica')
                    text.set('font-size', str(self.smallFontSize*0.8))
                    g.append(text)
                    if allele=='B':
                        shownAlleles=True

            if isAllelic:
                offset = self.chromosomeMargin*0.8
                r = self.canvas.getRectangle(currentX-offset, y-self.headerHeight*0.5,
                                      chromosomePixelWidth+offset*2, self.heightPerState+offset*2+self.headerHeight*0.5)
                r.set('z-index', '0')
                self.canvas.modifyStyle( r, {'fill':'#DCFFCE', 'stroke':'none'})
                self.stateGroup.insert(0,r)


            withinX = currentX

            binSizes = {}
            if logScale:
                currCoord = 0
                minBinSize = 1_000

                for binIndex in sorted(list(row[ row['chromosome']==chromosome ]['binIndex'])):

                    d = row[ row['chromosome']==chromosome ]['binIndex']==binIndex
                    dat = row[ row['chromosome']==chromosome ][d].iloc[0,:]
                    cn =dat['copyNumber']

                    space = dat['startCoordinate'] - currCoord
                    if space > minBinSize: # add intermediate bin:
                        binSizes[(binIndex, 'spacer')] = np.log(space)/10 if logRepeats else space/chromosomePixelWidth


                    size = np.log((dat['endCoordinate'] - dat['startCoordinate'] ))
                    binSizes[binIndex] = size
                    currCoord=dat['endCoordinate']

                space = chromosomeSize - currCoord
                if space > minBinSize: # add intermediate bin:
                    binSizes[(binIndex, 'spacerfinal')] = np.log(space)/10 if logRepeats else space/chromosomePixelWidth


            for binIndex in sorted(list(row[ row['chromosome']==chromosome ]['binIndex'])):
                d = row[ row['chromosome']==chromosome ]['binIndex']==binIndex
                dat = row[ row['chromosome']==chromosome ][d].iloc[0,:]
                cn =dat['copyNumber']

                if logScale:
                    if (binIndex, 'spacer') in binSizes:
                        withinX+=( binSizes[(binIndex,'spacer')] / sum( binSizes.values() ) ) * chromosomePixelWidth

                    size = ( binSizes[binIndex] / sum( binSizes.values() ) ) * chromosomePixelWidth


                else:
                    size = (dat['endCoordinate'] - dat['startCoordinate'] )*self.pixelsPerBase

                r = self.canvas.getRectangle(withinX, y, size, self.heightPerState)
                self.canvas.modifyStyle(r, {'stroke-width':'0.2'})
                fillAttr='fill'

                if chromosome.endswith('_A') or chromosome.endswith('_B'): # allelic
                    cn+=1 # (diploid color == white == 2 )


                if cn >= 4:
                    self.canvas.modifyStyle(r, {fillAttr:self.gainTwoColor})
                if cn == 3:
                    self.canvas.modifyStyle(r, {fillAttr:self.gainColor})
                if cn==2:
                    self.canvas.modifyStyle(r, {fillAttr:self.normalColor})
                if cn==1:
                    self.canvas.modifyStyle(r, {fillAttr:self.lossColor})
                if cn==0:
                    self.canvas.modifyStyle(r, {fillAttr:self.totalLossColor})
                if np.isnan(cn):
                    self.canvas.modifyStyle(r, {fillAttr:self.missingColor})

                g.append(r)

                withinX+=size

            currentX += chromosomePixelWidth+self.chromosomeMargin

        ##### plot the label:
        text = self.canvas.getText(label, currentX, y+0.5*self.heightPerState)
        text.set('text-anchor','begin')
        text.set('dominant-baseline','middle')
        text.set('font-family','Helvetica')
        text.set('font-size', str(self.bigFontSize))


        g.append(text)

        self.canvas.setWidth( max(self.canvas.width, currentX+self.chromosomeMargin+20+ self.offset[0] ) )
        self.canvas.setHeight( max(self.canvas.height, y+self.heightPerState*1.5+ self.offset[1] ))

        return {'x':x}
