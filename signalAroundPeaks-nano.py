##################################
#                                #
# Last modified 2019/05/22       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import gzip
import numpy as np
import os
from sets import Set

def run():

    if len(sys.argv) < 8:
        print 'usage: python %s inputfilename chrFieldID posField strandField radius window methylation_single_base.tsv outputfilename [-bismark.cov] [-bed] [-minCov N] [-unstranded] [-ERANGE_hts] [-narrowPeak] [-first number]' % sys.argv[0]
        print '\tInput format: <fields .. tabs> chr <tab> position <tab> strandField' 
        print '\tThis script outputs the average signal over all regions within the given radius' 
        print '\tif the the -bed option is used, the middle point of a bed region will be used; specifiy the posField as the left coordinate of the region' 
        print '\tif the the -narrowPeak option is used, the posFielf will be ignored and strand will be assumed to be +' 
        print '\tNote: the script will normalize only against the number of windows that have a CpG ot other signal in the methylation file' 
        print '\tuse the [-bismark.cov] option if you want to use the script for a Bismark bedcov output file'
        sys.exit(1)
    
    regionfilename = sys.argv[1]
    chrFieldID = int(sys.argv[2])
    posFieldID = int(sys.argv[3])
    strandFieldID = int(sys.argv[4])
    radius = int(sys.argv[5])
    window = int(sys.argv[6])
    metfilename = sys.argv[7]
    outfilename = sys.argv[8]

    doBismarkCov = False
    if '-bismark.cov' in sys.argv:
        doBismarkCov = True
        print 'will treat methylaiton file as a Bismark coverage file'

    doFirst=False
    if '-first' in sys.argv:
        firstN=int(sys.argv[sys.argv.index('-first')+1])
        doFirst=True
        print 'will only look at the first', firstN, 'locations'

    doNarrowPeak=False
    if '-narrowPeak' in sys.argv:
        doNarrowPeak=True
        print 'will treat regions as being in narrowPeak format'

    doBed=False
    if '-bed' in sys.argv:
        print 'will treat input as bed file and center around the midpoint of reigons'
        doBed=True

    noStrand = False
    if '-unstranded' in sys.argv:
        print 'will treat all regions as + strand'
        noStrand = True

    doERANGE=False
    if '-ERANGE_hts' in sys.argv:
        doERANGE=True

    MC = 1
    if '-minCov' in sys.argv:
        MC = int(sys.argv[sys.argv.index('-minCov')+1])
        print 'will discard positions with coverage less than', MC

    RegionDict={}
    ScoreDict={}
    if regionfilename.endswith('.gz'):
        listoflines = gzip.open(regionfilename)
    else:
        listoflines = open(regionfilename)
    k=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.replace('\x00','').strip().split('\t')
        if doNarrowPeak:
            pass
        else:
            if len(fields) < max(chrFieldID, posFieldID, strandFieldID, 3):
                continue
        k+=1
        if doFirst and k > firstN:
            continue
        if len(fields)<3:
           continue
        if doBed:
            chr=fields[chrFieldID]
            left=int(fields[posFieldID])
            right=int(fields[posFieldID+1])
            pos=int((right+left)/2.0)
            if noStrand:
                strand='+'
            else:
                strand=fields[strandFieldID]
        elif doERANGE:
            chr=fields[1]
            pos=int(fields[9])
            strand='+'
        elif doNarrowPeak:
            chr=fields[0]
            pos=int(fields[1]) + int(fields[9])
            strand='+'
        else:
            chr=fields[chrFieldID]
            pos=int(fields[posFieldID])
            if noStrand:
                strand='+'
            else:
                strand=fields[strandFieldID]
        RegionDict[(chr,pos,strand)]=[]
        if ScoreDict.has_key(chr):
            pass
        else:
            ScoreDict[chr]={}
        for i in range(pos-radius,pos+radius):
            ScoreDict[chr][i] = ''
        if k % 10000 == 0:
            print k

    print 'Finished importing regions', len(RegionDict.keys()), 'regions in total'
    print 'Importing methylation scores'    
    
    if metfilename.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + metfilename
    elif metfilename.endswith('gz'):
        cmd = 'zcat ' + metfilename
    else:
        cmd = 'cat ' + metfilename
    p = os.popen(cmd, "r")
    line = 'line'
    j=0
    while line != '':
        line = p.readline()
        if line.startswith('#'):
            continue
        if line.startswith('track type'):
            continue
        if line.startswith('chrom\tstart\tend\tmeth\tunmeth\tcov'):
            continue
        if line == '':
            break
        j+=1
        if j % 1000000 == 0:
            out = str(j/1000000) + 'M lines processed'
            print out
        if ' ' in line:
            fields=line.strip().split(' ')
        else:
            fields=line.strip().split('\t')
        chr = fields[0]
        start = int(float(fields[1]))
        stop = int(float(fields[2]))
        if stop == start:
            stop += 1
        if doBismarkCov:
            meth = float(fields[4])
            unmeth = float(fields[5])
            cov = meth + unmeth + 0.0
        else:
            meth = int(float(fields[3]))
            unmeth = int(float(fields[4]))
            cov = float(fields[5])
        if cov < MC:
            continue
        score = meth/cov
        for i in range(start,stop):
            if ScoreDict.has_key(chr): 
                if ScoreDict[chr].has_key(i):
                    ScoreDict[chr][i] = score

    print 'Finished importing methylation scores'    
    print 'Outputting final stats'    

    for (chr,pos,strand) in RegionDict.keys():
        for i in range(pos-radius,pos+radius):
            RegionDict[(chr,pos,strand)].append(ScoreDict[chr][i])
        if strand == 'R' or strand == '-':
            RegionDict[(chr,pos,strand)].reverse()

    keys=RegionDict.keys()
    keys.sort()

    NormDict = {}
    FinalDict = {}
    for i in range(0-radius,0+radius,window):
        FinalDict[i] = 0.0
        NormDict[i] = 0.0

    for (chr,pos,strand) in keys:
        for i in range(-radius,+radius,window):
            FScores = []
            for j in range(i-radius,i-radius+window):
                if RegionDict[(chr,pos,strand)][j] != '':
                    FScores.append(RegionDict[(chr,pos,strand)][j])
            if len(FScores) > 0:
                FinalDict[i] += np.mean(FScores)
                NormDict[i] += 1

    outfile=open(outfilename,'w')
    outline = '#pos\t1-meth\tmeth\tN_sites'
    outfile.write(outline + '\n')

    keys = FinalDict.keys()
    keys.sort()
    for i in keys:
#        print i, FinalDict[i], NormDict[i]
        if NormDict[i] > 0:
            outline = str(i) + '\t' + str(1 - FinalDict[i]/NormDict[i]) + '\t' + str(FinalDict[i]/NormDict[i]) + '\t' + str(NormDict[i])
        else:
            outline = str(i) + '\t' + 'nan' + '\t' + 'nan' + '\t' + '0'
        outfile.write(outline + '\n')

    outfile.close()
   
run()
