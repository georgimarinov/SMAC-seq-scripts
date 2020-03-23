##################################
#                                #
# Last modified 2018/11/26       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import numpy as np
import os
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sets import Set
import matplotlib, copy
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.collections import PatchCollection
import random

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def run():

    if len(sys.argv) < 8:
        print 'usage: python %s methylation_reads_all.tsv peak_list chrFieldID leftFieldID rightFieldID strandFieldID tabix_path outfile_prefix [-resize factor] [-subset N] [-minCov fraction] [-minPassingBases fraction] [-minReads N] [-unstranded] [-minAbsLogLike float] [-scatterPlot colorscheme minScore maxScore color|none] [-window bp] [-binarize threshold]' % sys.argv[0]
        print '\tUse the [-subset] option if you want only N of the fragments; the script will pick the N fragments best covering each region, and will discard regions with fewer than N covering fragments'
        print '\tThe [-heatmap] option will generate png heatmaps instead of text file matrices'
        print '\tThe [-minCov] option will remove all fragments that cover the region at less than the specified fraction'
        sys.exit(1)

    reads_file = sys.argv[1]
    peaks = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    strandFieldID = int(sys.argv[6])
    tabix = sys.argv[7]
    outprefix = sys.argv[8]

    print 'finished inputting genomic sequence'

    minAbsLogLike = 0
    if '-minAbsLogLike' in sys.argv:
        minAbsLogLike = float(sys.argv[sys.argv.index('-minAbsLogLike') + 1])
        print 'will omit bases with absolute loglike values less than', minAbsLogLike

    doBinary = False
    if '-binarize' in sys.argv:
        doBinary = True
        BinaryThreshold = float(sys.argv[sys.argv.index('-binarize') + 1])
        print 'will binarize scores at the ', BinaryThreshold, 'threshold'

    Nsub = 1
    doSubset = False
    if '-subset' in sys.argv:
        doSubset = True
        Nsub = int(sys.argv[sys.argv.index('-subset')+1])
        print 'will only output', Nsub, ' random complete fragments'

    window = 1
    if '-window' in sys.argv:
        window = int(sys.argv[sys.argv.index('-window')+1])
        print 'will average scores across windows of size', window, 'bp'

    MR = 0
    if '-minReads' in sys.argv:
        MR = int(sys.argv[sys.argv.index('-minReads')+1])
        print 'will only output regions with at least', MR, 'reads covering the region'

    MFC = 0
    if '-minCov' in sys.argv:
        MFC = float(sys.argv[sys.argv.index('-minCov')+1])
        print 'will only output fragments with at least', MFC, 'fractional coverage of the region'

    MFCB = 0
    if '-minPassingBases' in sys.argv:
        MFCB = float(sys.argv[sys.argv.index('-minPassingBases')+1])
        print 'will only output fragments with at least', MFCB, 'fractional coverage of the region'

    doScatterPlot = False
    if '-scatterPlot' in sys.argv:
        doScatterPlot = True
        print 'will output scatter plot'
        SPcs = sys.argv[sys.argv.index('-scatterPlot')+1]
        SPmin = float(sys.argv[sys.argv.index('-scatterPlot')+2])
        SPmax = float(sys.argv[sys.argv.index('-scatterPlot')+3])
        SPedge = sys.argv[sys.argv.index('-scatterPlot')+4]
        resize = 1
        if '-resize' in sys.argv:
            resize = float(sys.argv[sys.argv.index('-resize') + 1])
            print 'will resize scatter plots by a factor of', resize

    doNS = False
    if '-unstranded' in sys.argv:
        doNS = True
        print 'will not treat regions as stranded'

    RegionDict = {}

    if peaks.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + peaks
    elif peaks.endswith('.gz'):
        cmd = 'gunzip -c ' + peaks
    elif peaks.endswith('.zip'):
        cmd = 'unzip -p ' + peaks
    else:
        cmd = 'cat ' + peaks
    RN = 0
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        if line == '':
            break
        if line.startswith('#'):
            continue
        linefields = line.strip().split('\t')
        RN += 1
        chr = linefields[chrFieldID]
        start = max(0,int(linefields[leftFieldID]))
        end = int(linefields[rightFieldID])
        if doNS:
            strand = '+'
        else:
            strand = linefields[strandFieldID]
        if RegionDict.has_key(chr):
            pass
        else:
            RegionDict[chr]={}
        if RegionDict[chr].has_key(start):
            pass
        else:
            RegionDict[chr][start] = []
        RegionDict[chr][start].append((end,strand))

    print 'finished inputting regions'

    chrList = RegionDict.keys()
    Matrix = []
    chrList.sort()
    for chr in chrList:
        posList = RegionDict[chr].keys()
        posList.sort()
        for RL in posList:
            for (RR,RS) in RegionDict[chr][RL]:
                start = RL
                end = RR
                cmd = tabix + ' ' + reads_file + ' ' + chr + ':' + str(RL) + '-' + str(RR)
                p = os.popen(cmd, "r")
                reads = []
                line = 'line'
                while line != '':
                    line = p.readline().strip()
                    if line == '':
                        break
                    fields = line.strip().split('\t')
                    chr = fields[0]
                    read_left = int(fields[1])
                    read_right = int(fields[2])
                    strand = fields[3]
                    read = fields[4]
                    cgs = fields[6]
                    loglike = fields[7]
                    if (min(read_right,RR) - max(read_left,RL))/(RR - RL + 0.0) >= MFC:
                        reads.append((chr,read_left,read_right,strand,read,cgs,loglike))
                print chr, RL, RR, len(Matrix)
                maxScoresLength = 0
                for (chr,readleft,readright,readstrand,read,cgs,loglike) in reads:
                    scores = []
                    CGscoreDict = {}
                    CGs = cgs.split(',')
                    LLs = loglike.split(',')
                    for i in range(len(CGs)):
                        CG = int(CGs[i])
                        LL = float(LLs[i])
                        CGscoreDict[CG] = LL
                    MS = 0.0
                    for i in range(RL,RR):
                        if CGscoreDict.has_key(i):
                            MS += 1
                    maxScoresLength = max(MS,maxScoresLength)
                for (chr,readleft,readright,readstrand,read,cgs,loglike) in reads:
                    scores = []
                    CGscoreDict = {}
                    CGs = cgs.split(',')
                    LLs = loglike.split(',')
                    for i in range(len(CGs)):
                        CG = int(CGs[i])
                        LL = float(LLs[i])
                        if math.fabs(LL) < minAbsLogLike:
                            pass
                        else:
                            CGscoreDict[CG] = LL
                    MS = 0.0
                    for i in range(RL,RR):
                        if CGscoreDict.has_key(i):
                            scores.append(CGscoreDict[i])
                            MS += 1
                        else:
                            scores.append(0)
                    scoresset = list(Set(scores))
                    if len(scoresset) == 1 and scoresset[0] == 0:
                        continue
                    if MS/maxScoresLength < MFCB:
                        continue
                    if RS == '-':
                        scores.reverse()
                    Matrix.append((MS,scores))

    if doSubset:
        Matrix = random.sample(Matrix,Nsub)

    Matrix.sort()
    Matrix.reverse()
    NewMatrix = []
    for (MS,scores) in Matrix:
        newscores = []
        for i in range(0,len(scores),window):
            s = 0.0
            ss = []
            for j in range(i,min(i+window,len(scores))):
                if scores[j] == 0:
                    pass
                else:
                    ss.append(scores[j])
            if len(ss) > 0:
                s += np.mean(ss)
            if doBinary:
                if s != 0:
                    if s >= BinaryThreshold:
                        newscores.append(1.0)
                    else:
                        newscores.append(0.01)
                else:
                    newscores.append(0)
            else:    
                newscores.append(s)
        NewMatrix.append(newscores)

    print len(NewMatrix), 'reads retained'

    X = []
    for scores in NewMatrix:
        X.append(scores)
    Z = linkage(X, method='ward', metric='euclidean', optimal_ordering=True)
    clusters = fcluster(Z, 0, criterion='distance')
    CDict = {}
    for i in range(len(clusters)):
        C = clusters[i]
        if CDict.has_key(C):
            pass
        else:
            CDict[C] = []
        CDict[C].append(i)
    Cs = CDict.keys()
    Cs.sort()

    if doScatterPlot:
        Xaxis = []
        Yaxis = []
        ColorScores = []

        i=0
        for C in Cs:
            for k in CDict[C]:
                i+=1
                scores = X[k]
                j = RL
                for s in scores:
                    j += window
                    Xaxis.append(j)
                    Yaxis.append(i)
                ColorScores += scores

        newXaxis = []
        newYaxis = []
        newColorScores = []
        for i in range(len(ColorScores)):
            if math.fabs(ColorScores[i]) > 0.0:
                newColorScores.append(ColorScores[i])
                newXaxis.append(Xaxis[i])
                newYaxis.append(Yaxis[i])

        aspectRatio = ((RR - RL)/(window+0.0))/len(NewMatrix)

        rect = 0.10,0.10,0.8,0.8
        fig = figure(figsize=(20*resize, 20*resize/aspectRatio))
        ax = fig.add_subplot(1,1,1,aspect='equal')
        ax = fig.add_axes(rect)
        lowerlimitX = min(Xaxis) - 2*window
        upperlimitX = max(Xaxis) + 2*window
        lowerlimitY = min(Yaxis) - 1
        upperlimitY = max(Yaxis) + 1
        if SPedge == 'none':
            ax.scatter(newXaxis, newYaxis, marker='o', c=newColorScores, vmin=SPmin, vmax=SPmax, cmap=SPcs)
        else:
            ax.scatter(newXaxis, newYaxis, marker='o', edgecolor=SPedge, c=newColorScores, vmin=SPmin, vmax=SPmax, cmap=SPcs)
        ax.set_xlim(lowerlimitX,upperlimitX)
        ax.set_ylim(lowerlimitY,upperlimitY)
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()

        savefig(outprefix + '.scatter.png')
        savefig(outprefix + '.scatter.eps', format='eps')

run()

