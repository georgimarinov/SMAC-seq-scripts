##################################
#                                #
# Last modified 2018/07/05       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import beta
from scipy.stats import binom
import random
import os
import math
from sets import Set

def run():

    if len(sys.argv) < 9:
        print 'usage: python %s methylation_reads_all.tsv peaks chrFieldID leftFiled RightFieldID minCoverage maxDist tabix_location outfile [-subsample N]' % sys.argv[0]
        print '\Note: the script assumes Tombo 1.3 probabilities, a tabix indexed reads file, and uses a beta distribution prior of (10,10) by default'
        print '\Note: the subsample option will sample the reads in all comparisons down to the minCoverage level; the N parameter indicates how many such subsamplings should be averaged for the final value'
        sys.exit(1)

    reads = sys.argv[1]
    peaks = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    minCov = int(sys.argv[6])
    maxDist = int(sys.argv[7])
    tabix = sys.argv[8]
    outfilename = sys.argv[9]

    doSS = False
    if '-subsample' in sys.argv:
        SS = int(sys.argv[sys.argv.index('-subsample') + 1])
        doSS = True
        print 'will subsample all comparisons down to', minCov, 'reads'
        print 'will take the average outcome of', SS, 'subsamplings'

    alph = 10
    bet = 10
    PSS = 100

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
        if PeakDict.has_key(chr):
            pass
        else:
            PeakDict[chr] = []
        PeakDict[chr].append((RL,RR))

    print 'finished inputting peaks'

    outfile = open(outfilename,'w')
    outline = '#chr\tpeak1_left\tpeak1_right\tpeak1_open\tpeak1_closed\tpeak1_fraction\tpeak2_left\tpeak2_right\tpeak2_open\tpeak2_closed\tpeak2_fraction\tp_val'
    outfile.write(outline + '\n')

    for chr in PeakDict.keys():
        PeakDict[chr].sort()
        for i in range(len(PeakDict[chr])-1):
            (RL1,RR1) = PeakDict[chr][i]
            for j in range(i+1,len(PeakDict[chr])):
                (RL2,RR2) = PeakDict[chr][j]
                print chr,RL1,RR1,RL2,RR2,
                if RR2 - RL1 > maxDist:
                    break
                cmd = tabix + ' ' + reads + ' ' + chr + ':' + str(RL1) + '-' + str(RR2)
                p = os.popen(cmd, "r")
                RegionReads = []
                line = 'line'
                while line != '':
                    line = p.readline().strip()
                    if line == '':
                        break
                    fields = line.strip().split('\t')
                    chr = fields[0]
                    read_left = int(fields[1])
                    read_right = int(fields[2])
                    if read_left <= RL1 and read_right >= RR2:
                        pass
                    else:
                        continue
                    strand = fields[3]
                    read = fields[4]
                    cgs = fields[6].split(',')
                    loglike = fields[7].split(',')
                    RegionReads.append((cgs,loglike))
                print len(RegionReads)
                if len(RegionReads) < minCov:
                    continue
                if doSS:
                    p_av = 0.0
                    CR1_av = 0.0
                    OR1_av = 0.0
                    CR2_av = 0.0
                    OR2_av = 0.0
                    for S in range(SS):
                        RegionReadsSubSampled = random.sample(RegionReads,minCov)
                        OpenOrClosed = []
                        for (cgs,loglike) in RegionReadsSubSampled:
                            t = zip(cgs,loglike)
                            RD = dict((int(x), float(y)) for x, y in t)
                            (A,B) = (alph,bet)
                            for pos in range(RL1,RR1):
                                if RD.has_key(pos):
                                    p = RD[pos]
                                    Z = int(PSS*p)
                                    A = A + Z
                                    B = B + PSS - Z
                            if beta.mean(A,B) > 0.5:
                                final_p1 = 1
                            else:
                                final_p1 = 0
                            (A,B) = (alph,bet)
                            for pos in range(RL2,RR2):
                                if RD.has_key(pos):
                                    p = RD[pos]
                                    Z = int(PSS*p)
                                    A = A + Z
                                    B = B + PSS - Z
                            if beta.mean(A,B) > 0.5:
                                final_p2 = 1
                            else:
                                final_p2 = 0
                            OpenOrClosed.append((final_p1,final_p2))
                        C00 = OpenOrClosed.count((0,0))
                        C01 = OpenOrClosed.count((0,1))
                        C10 = OpenOrClosed.count((1,0))
                        C11 = OpenOrClosed.count((1,1))
                        CR1 = C01 + C00
                        OR1 = C10 + C11
                        CR2 = C10 + C00
                        OR2 = C01 + C11
                        oddsratio, pvalue = fisher_exact([[C00, C01], [C10, C11]])
                        logp = -math.log10(pvalue)	
                        p_av += logp
                        CR1_av += CR1
                        OR1_av += OR1
                        CR2_av += CR2
                        OR2_av += OR2
                    (pvalue, CR1, OR1, CR2, OR2) = (p_av/SS,CR1_av/SS,OR1_av/SS,CR2_av/SS,OR2_av/SS)
                    outline = chr + '\t' + str(RL1) + '\t' + str(RR1)
                    outline = outline + '\t' + str(OR1) + '\t' + str(CR1) + '\t' + str(OR1/(OR1 + CR1 + 0.0))
                    outline = outline + '\t' + str(RL2) + '\t' + str(RR2)
                    outline = outline + '\t' + str(OR2) + '\t' + str(CR2) + '\t' + str(OR2/(OR2 + CR2 + 0.0)) + '\t' + str(pvalue)
                    outfile.write(outline + '\n')
                else:
                    OpenOrClosed = []
                    for (cgs,loglike) in RegionReads:
                        t = zip(cgs,loglike)
                        RD = dict((int(x), float(y)) for x, y in t)
                        (A,B) = (alph,bet)
                        for pos in range(RL1,RR1):
                            if RD.has_key(pos):
                                p = RD[pos]
                                Z = int(PSS*p)
                                A = A + Z
                                B = B + PSS - Z
                        if beta.mean(A,B) > 0.5:
                            final_p1 = 1
                        else:
                            final_p1 = 0
                        (A,B) = (alph,bet)
                        for pos in range(RL2,RR2):
                            if RD.has_key(pos):
                                p = RD[pos]
                                Z = int(PSS*p)
                                A = A + Z
                                B = B + PSS - Z
                        if beta.mean(A,B) > 0.5:
                            final_p2 = 1
                        else:
                            final_p2 = 0
                        OpenOrClosed.append((final_p1,final_p2))
                    C00 = OpenOrClosed.count((0,0))
                    C01 = OpenOrClosed.count((0,1))
                    C10 = OpenOrClosed.count((1,0))
                    C11 = OpenOrClosed.count((1,1))
                    OR1 = C01 + C00
                    CR1 = C10 + C11
                    OR2 = C10 + C00
                    CR2 = C01 + C11
                    oddsratio, pvalue = fisher_exact([[C00, C01], [C10, C11]])
                    outline = chr + '\t' + str(RL1) + '\t' + str(RR1)
                    outline = outline + '\t' + str(OR1) + '\t' + str(CR1) + '\t' + str(OR1/(OR1 + CR1 + 0.0))
                    outline = outline + '\t' + str(RL2) + '\t' + str(RR2)
                    if pvalue == 0:
                        pvalue = 1e-300
                    outline = outline + '\t' + str(OR2) + '\t' + str(CR2) + '\t' + str(OR2/(OR2 + CR2 + 0.0)) + '\t' + str(-math.log10(pvalue))
                    outfile.write(outline + '\n')

    outfile.close()
            
run()

