##################################
#                                #
# Last modified 2018/12/03       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import numpy as np
from scipy.stats import entropy
from scipy.stats import fisher_exact
from scipy.stats import beta
from scipy.stats import binom
from sklearn.metrics import normalized_mutual_info_score as NMIS
import random
import os
import math
from sets import Set

def run():

    if len(sys.argv) < 10:
        print 'usage: python %s methylation_reads_all.tsv region.bed chrFieldID leftField rightFieldID minCoverage windowsize stepsize tabix_location outfileprefix [-subsample N] [-expectedMaxDist bp] [-label fieldID]' % sys.argv[0]
        print '\Note: the script assumes Tombo 1.3 probabilities, and a tabix indexed reads file'
        print '\Note: the [-subsample] option will sample the reads in all comparisons down to the minCoverage level; the N parameter indicates how many such subsamplings should be averaged for the final value'
        print '\Note: the [-expectedMaxDist] option will change the initial window over which the required minimum number of reads is to be search for; default: 2kb'
        sys.exit(1)

    reads = sys.argv[1]
    peaks = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    minCov = int(sys.argv[6])
    window = int(sys.argv[7])
    step = int(sys.argv[8])
    tabix = sys.argv[9]
    outprefix = sys.argv[10]

    alph = 10
    bet = 10
    PSS = 100

    SS = 1
    doSS = False
    if '-subsample' in sys.argv:
        SS = int(sys.argv[sys.argv.index('-subsample') + 1])
        doSS = True
        print 'will subsample all comparisons down to', minCov, 'reads'
        print 'will take the average outcome of', SS, 'subsamplings'

    doLabel = False
    if '-label' in sys.argv:
        labelFieldID = int(sys.argv[sys.argv.index('-label') + 1])
        doLabel = True

    EMD = 2000
    if '-expectedMaxDist' in sys.argv:
        EMD = int(sys.argv[sys.argv.index('-expectedMaxDist') + 1])
        print 'will use an expected maximum distance of', EMD

    PeakDict = {}
    if peaks.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + peaks
    elif peaks.endswith('.gz') or peaks.endswith('.bgz'):
        cmd = 'zcat ' + peaks
    elif peaks.endswith('.zip'):
        cmd = 'unzip -p ' + peaks
    else:
        cmd = 'cat ' + peaks
    RN = 0
    P = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = P.readline().strip()
        if line == '':
            break
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        RL = int(fields[leftFieldID])
        RR = int(fields[rightFieldID])
        if doLabel:
            label = fields[labelFieldID]
        else:
            label = chr + '-' + str(RL) + '-' + str(RR)
        if PeakDict.has_key(chr):
            pass
        else:
            PeakDict[chr] = []
        PeakDict[chr].append((RL,RR,label))

    print 'finished inputting peaks'

    for chr in PeakDict.keys():
        PeakDict[chr].sort()
        for (left,right,label) in PeakDict[chr]:
            print chr,left,right,label
            Matrix = {}
            outfile = open(outprefix + '.' + label,'w')
            for i in range(left,right,step):
                RC = 0
                cmd = tabix + ' ' + reads + ' ' + chr + ':' + str(i) + '-' + str(min(right,i + EMD))
                p = os.popen(cmd, "r")
                line = 'line'
                while line != '':
                    line = p.readline().strip()
                    if line == '':
                        break
                    fields = line.strip().split('\t')
                    read_left = int(fields[1])
                    read_right = int(fields[2])
                    if read_left <= i and read_right >= min(right,i + EMD):
                        RC += 1
                    else:
                        continue
                print RC, left, right
                if RC >= minCov:
                    RCC = minCov
                    LJ = min(right,i + EMD)
                    while RCC >= minCov and LJ < right:
                        LJ += window
                        RRR = 0
                        cmd = tabix + ' ' + reads + ' ' + chr + ':' + str(i) + '-' + str(LJ)
                        p = os.popen(cmd, "r")
                        line = 'line'
                        while line != '':
                            line = p.readline().strip()
                            if line == '':
                                break
                            fields = line.strip().split('\t')
                            read_left = int(fields[1])
                            read_right = int(fields[2])
                            if read_left <= i and read_right >= LJ:
                                RRR += 1
                            else:
                                continue
                        RCC = RRR
                    rightLimit = LJ - window
                else: 
                    RCC = RC
                    LJ = min(right,i + EMD)
                    while RCC <= minCov and LJ > i:
                        LJ = LJ - window
                        RRR = 0
                        cmd = tabix + ' ' + reads + ' ' + chr + ':' + str(i) + '-' + str(LJ)
                        p = os.popen(cmd, "r")
                        line = 'line'
                        while line != '':
                            line = p.readline().strip()
                            if line == '':
                                break
                            fields = line.strip().split('\t')
                            read_left = int(fields[1])
                            read_right = int(fields[2])
                            if read_left <= i and read_right >= LJ:
                                RRR += 1
                            else:
                                continue
                        RCC = RRR
#                        print RCC, i, LJ
                    rightLimit = LJ
                rightLimit = min(right,rightLimit)
                print chr, i, rightLimit
                RegionReads = []
                cmd = tabix + ' ' + reads + ' ' + chr + ':' + str(i) + '-' + str(rightLimit)
                p = os.popen(cmd, "r")
                line = 'line'
                while line != '':
                    line = p.readline().strip()
                    if line == '':
                        break
                    fields = line.strip().split('\t')
                    chr = fields[0]
                    read_left = int(fields[1])
                    read_right = int(fields[2])
                    if read_left <= i and read_right >= rightLimit:
                        pass
                    else:
                        continue
                    strand = fields[3]
                    read = fields[4]
                    try:
                        cgs = fields[6].split(',')
                        loglike = fields[7].split(',')
                        RegionReads.append((cgs,loglike))
                    except:
                        continue
                AccDict = {}
                for j in range(i,rightLimit,window):
                    AccDict[j] = []
                    Matrix[j] = {}
                    for k in range(j,rightLimit,window):
                        Matrix[j][k] = 0
                print len(RegionReads)
                if len(RegionReads) < minCov:
                    continue
                for S in range(SS):
                    if doSS:
                        RegionReadsSubSampled = random.sample(RegionReads,minCov)
                    else:
                        RegionReadsSubSampled = RegionReads
                    for (cgs,loglike) in RegionReadsSubSampled:
                        t = zip(cgs,loglike)
                        RD = dict((int(x), float(y)) for x, y in t)
                        for j in range(i,rightLimit,window):
                            (A,B) = (alph,bet)
                            for pos in range(j, j + window):
                                if RD.has_key(pos):
                                    p = RD[pos]
                                    Z = int(PSS*p)
                                    A = A + Z
                                    B = B + PSS - Z
                            if beta.mean(A,B) > 0.5:
                                final_p = 1
                            else:
                                final_p = 0
                            AccDict[j].append(final_p)
                    for j in range(i,rightLimit,window):
                        print j
                        for k in range(j,rightLimit,window):
                            P = np.array(AccDict[j])
                            Q = np.array(AccDict[k])
                            NMI = NMIS(P,Q)
                            matches = np.sum(P==Q)
                            if matches < 0.5*len(P):
                                NMI = -NMI
                            Matrix[j][k] += NMI/SS
            outline = '#'
            for i in range(left,right,window):
                outline = outline + '\t' + str(i)
            outfile.write(outline + '\n')
            for i in range(left,right,window):
                outline = str(i)
                for j in range(left,right,window):
                    if Matrix.has_key(i):
                        if Matrix[i].has_key(j):
                            outline = outline + '\t' + "{0:.2f}".format(Matrix[i][j])
                        else:
                            outline = outline + '\t' + 'nan'
                    else:
                        outline = outline + '\t' + 'nan'
                outfile.write(outline + '\n')
            outfile.close()
            
run()

