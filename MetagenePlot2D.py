##################################
#                                #
# Last modified 2018/08/20       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import numpy as np
import os
from sets import Set
import gzip

def run():

    if len(sys.argv) < 9:
        print 'usage: python %s inputfilename chrFieldID leftFieldID rightFieldID strandField window_size step_size NMI/JSD_file_list outputfilename -unstranded' % sys.argv[0]
        print '\tnote: the step size should be the same used to calculate the NMI/HSD matrix'
        print '\tnote: assumed NMI file name format:'
        print '\tNMI.100cov_50bp_1000bp.chrII_270000_320000, i.e. the script will split by dots and then by spaces to get the chromosome name'
        print '\tby default the script will print out the average matrix across genes longer than 2*window size, showing a window of size |window| around the TSS and the TTS'
        sys.exit(1)
    
    regionfilename = sys.argv[1]
    chrFieldID = int(sys.argv[2])
    leftFieldID = int(sys.argv[3])
    rightFieldID = int(sys.argv[4])
    strandFieldID = int(sys.argv[5])
    window = int(sys.argv[6])
    StepSize = int(sys.argv[7])
    NMIfiles = sys.argv[8]
    outfilename = sys.argv[9]

    noStrand = False
    if '-unstranded' in sys.argv:
        print 'will treat all regions as + strand'
        noStrand = True

    RegionDict = {}
    listoflines = open(regionfilename)
    k=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.replace('\x00','').strip().split('\t')
        chr=fields[chrFieldID]
        left=int(fields[leftFieldID])
        right=int(fields[rightFieldID])
        if noStrand:
            strand='+'
        else:
            strand=fields[strandFieldID]
        if RegionDict.has_key(chr):
            pass
        else:
            RegionDict[chr]=[]
        RegionDict[chr].append((left,right,strand))

    print 'Finished importing regions'

    MatrixDict = {}

    filelinelist = open(NMIfiles)
    for fileline in filelinelist:
        file = fileline.strip().split('\t')[0]
#        chr = file.split('/')[-1].split('.')[2].split('_')[0]
        chr = file.split('/')[-1].split('.')[-2].split('_')[0]
        if RegionDict.has_key(chr):
            pass
        else:
            continue
        if MatrixDict.has_key(chr):
            pass
        else:
            MatrixDict[chr] = {}
        print chr, file
        if file.endswith('.gz'):
            linelist = gzip.open(file)
        else:
            linelist = open(file)
        for line in linelist:
            fields = line.strip().split('\t')
            if line.startswith('#'):
                posDict = {}
                for i in range(1,len(fields)):
                    posDict[i] = int(fields[i])
                continue
            try:
                pos = int(fields[0])
            except:
                print 'skipping line'
                print fields
                continue
            if MatrixDict[chr].has_key(pos):
                pass
            else:
                MatrixDict[chr][pos] = {}
            for i in range(1,len(fields)):
                if fields[i] == 'nan' or fields[i] == 'na' or fields[i] == 'n' or fields[i] == '':
                    continue
                score = float(fields[i]) 
                p2 = posDict[i]
                if MatrixDict[chr][pos].has_key(p2):
                    pass
                else:
                    MatrixDict[chr][pos][p2] = []
                MatrixDict[chr][pos][p2].append(score)

    print 'finished parsing matrices'

    for chr in MatrixDict.keys():
        for pos in MatrixDict[chr].keys():
            for p2 in MatrixDict[chr][pos].keys():
                newscore = sum(MatrixDict[chr][pos][p2])/len(MatrixDict[chr][pos][p2])
                MatrixDict[chr][pos][p2] = newscore

    print 'finished normalizing matrices'

    AverageGeneMatrix = {}
    for i in range(-window/2,2*window,StepSize):
        AverageGeneMatrix[i] = {}
        for j in range(-window/2,2*window,StepSize):
            AverageGeneMatrix[i][j] = []

    for chr in RegionDict.keys():
        if MatrixDict.has_key(chr):
            pass
        else:
            continue
        for (L,R,strand) in RegionDict[chr]:
            if R - L < window:
                continue
            print L,R,strand
            if strand == '+':
                for i in range(L - window/2, L + window/2, StepSize):
                    p1 = i - (i % StepSize)
                    a1 = i - L
                    for j in range(i, L + window/2, StepSize):
                        p2 = j - (j % StepSize)
                        a2 = j - L
#                        AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        except:
                            print 'problem parsing matrix', a1, a2, chr, p1, p2, strand
                            pass
                    for j in range(R - window/2, R+window/2, StepSize):
                        p2 = j - (j % StepSize)
                        a2 = j - R + window + window/2
#                        AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        except:
                            print 'problem parsing matrix', a1, a2, chr, p1, p2, strand
                            pass
                for i in range(R - window/2, R + window/2, StepSize):
                    p1 = i - (i % StepSize)
                    a1 = i - R + window + window/2
                    for j in range(i, R + window/2, StepSize):
                        p2 = j - (j % StepSize)
                        a2 = j - R + window + window/2
#                        AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p1][p2])
                        except:
                            pass
            if strand == '-':
                for i in range(L + window/2, L - window/2, -StepSize):
                    p1 = i - (i % StepSize)
                    a1 = L - i
                    for j in range(i, L - window/2, -StepSize):
                        p2 = j - (j % StepSize)
                        a2 = L - j
#                        AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        except:
                            pass
                    for j in range(R + window/2, R - window/2, -StepSize):
                        p2 = j - (j % StepSize)
                        a2 = R - j + window + window/2
#                       AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        except:
                            pass
                for i in range(R + window/2, R - window/2, -StepSize):
                    p1 = i - (i % StepSize)
                    a1 = R - i + window + window/2
                    for j in range(i, R - window/2, -StepSize):
                        p2 = j - (j % StepSize)
                        a2 = R - j + window + window/2
#                       AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        try:
                            AverageGeneMatrix[a1][a2].append(MatrixDict[chr][p2][p1])
                        except:
                            pass
                    
    print 'finished compiling averaged matrix'

    for p1 in AverageGeneMatrix.keys():
        for p2 in AverageGeneMatrix[p1].keys():
            newlist = [b for b in AverageGeneMatrix[p1][p2] if str(b) != 'nan']
            newlist.sort()
            AverageGeneMatrix[p1][p2] = newlist
#            print p1, p2, newlist
            if len(AverageGeneMatrix[p1][p2]) > 0:
                normscore = sum(AverageGeneMatrix[p1][p2])/(len(AverageGeneMatrix[p1][p2]) + 0.0)
                AverageGeneMatrix[p1][p2] = normscore
            else:
                AverageGeneMatrix[p1][p2] = 'nan'

    outfile = open(outfilename,'w')
    outline = '#pos'
    for i in range(-window/2,2*window,StepSize):
        outline = outline + '\t' + str(i)
    outfile.write(outline + '\n')

    for i in range(-window/2,2*window,StepSize):
        outline = str(i)
        for j in range(-window/2,2*window,StepSize):
            if AverageGeneMatrix[i][j] == 'nan':
                outline = outline + '\t' + AverageGeneMatrix[i][j]
            else:
                outline = outline + '\t' + "{0:.3f}".format(AverageGeneMatrix[i][j])
        outfile.write(outline + '\n')

    outfile.close()
   
run()
