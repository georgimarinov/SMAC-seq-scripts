##################################
#                                #
# Last modified 2018/11/27       # 
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

    if len(sys.argv) < 8:
        print 'usage: python %s inputfilename chrFieldID posfieldID strandFieldID window_size step_size NMI/JSD_file_list outputfilename [-unstranded]' % sys.argv[0]
        print '\tnote: the step size should be the same used to calculate the NMI/HSD matrix'
        print '\tnote: assumed NMI file name format:'
        print '\tNMI.100cov_50bp_1000bp.chrII_270000_320000, i.e. the script will split by dots and then by spaces to get the chromosome name'
        sys.exit(1)
    
    regionfilename = sys.argv[1]
    chrFieldID = int(sys.argv[2])
    posFieldID = int(sys.argv[3])
    strandFieldID = int(sys.argv[4])
    window = int(sys.argv[5])
    StepSize = int(sys.argv[6])
    NMIfiles = sys.argv[7]
    outfilename = sys.argv[8]

    doUnstranded = False
    if '-unstranded' in sys.argv:
        doUnstranded = True

    RegionDict = {}
    listoflines = open(regionfilename)
    k=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        pos = int(fields[posFieldID])
        if doUnstranded:
            strand = '+'
        else:
            strand = fields[strandFieldID]
        if RegionDict.has_key(chr):
            pass
        else:
            RegionDict[chr] = []
        RegionDict[chr].append((pos,strand))

    print 'Finished importing regions'

    MatrixDict = {}

    filelinelist = open(NMIfiles)
    for fileline in filelinelist:
        file = fileline.strip().split('\t')[0]
        chr = file.split('/')[-1].split('.')[-2].split('-')[0]
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
#                print fields[i]
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
                if len(MatrixDict[chr][pos][p2]) > 0:
                    newscore = sum(MatrixDict[chr][pos][p2])/len(MatrixDict[chr][pos][p2])
                    MatrixDict[chr][pos][p2] = newscore
                else:
                    del MatrixDict[chr][pos][p2]

    print 'finished normalizing matrices'

    outfile = open(outfilename,'w')
    outline = '#pos'
    for i in range(-window,window,StepSize):
        outline = outline + '\t' + str(i)
    outfile.write(outline + '\n')

    for chr in RegionDict.keys():
        if MatrixDict.has_key(chr):
            pass
        else:
            continue
        for (pos,strand) in RegionDict[chr]:
            newpos = pos - pos % StepSize
            final_array = []
            IsFull = True
#            print chr, pos, newpos
            if strand == '+':
                for i in range(newpos - window, newpos + window, StepSize):
                    score = 'nan'
                    try:
                        score = MatrixDict[chr][newpos][i]
                    except:
                        try:
                            score = MatrixDict[chr][i][newpos]
                        except:
                            pass
#                            score = 'nan'
#                            IsFull = False
#                            break
                    final_array.append(score)
            if strand == '-':
                for i in range(newpos + window, newpos - window, -StepSize):
                    score = 'nan'
                    try:
                        score = MatrixDict[chr][newpos][i]
                    except:
                        try:
                            score = MatrixDict[chr][i][newpos]
                        except:
                            pass
#                            score = 'nan'
#                            IsFull = False
#                            break
                    final_array.append(score)
#            if not IsFull:
#                continue
            outline = chr + ':' + str(pos) + ':' + strand
            for score in final_array:
                outline = outline + '\t' + str(score)
            outfile.write(outline + '\n')

    outfile.close()
   
run()
