##################################
#                                #
# Last modified 2018/11/30       # 
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
from scipy.stats import norm
import random
import os
import math
from sets import Set
from sklearn.metrics import normalized_mutual_info_score as NMIS

def run():

    if len(sys.argv) < 10:
        print 'usage: python %s methylation_reads_all.tsv peaks chrFieldID leftFiled RightFieldID minCoverage maxDist N_samplings tabix_location outfile [-subsample N] [-quantiles N]' % sys.argv[0]
        print '\Note: the script assumes Tombo 1.3 probabilities, a tabix indexed reads file, and uses a beta distribution prior of (10,10) by default'
        print '\Note: the subsample option will sample the reads in all comparisons down to the minCoverage level; the N parameter indicates how many such subsamplings should be averaged for the final value'
        print '\Note: the subsample option IS REQUIRED AT THE MOMENT'
        sys.exit(1)

    reads = sys.argv[1]
    peaks = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    minCov = int(sys.argv[6])
    maxDist = int(sys.argv[7])
    Nsamp = int(sys.argv[8])
    tabix = sys.argv[9]
    outfilename = sys.argv[10]

    QU = 5
    if '-quantiles' in sys.argv:
        QU = int(sys.argv[sys.argv.index('-quantiles') + 1])
        print 'will split reads into', QU, 'quantiles/bins instead of the default 5'

    doSS = False
    if '-subsample' in sys.argv:
        SS = int(sys.argv[sys.argv.index('-subsample') + 1])
        doSS = True
        print 'will subsample all comparisons down to', minCov, 'reads'
        print 'will take the average outcome of', SS, 'subsamplings'
        if minCov % QU != 0:
            print 'minCov value must be divisble by the number of quantiles, which is', QU, ', exiting'
            sys.exit(1)

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
    outline = '#chr\tpeak1_left\tpeak1_right\tpeak1_open\tpeak1_closed\tpeak1_fraction\tpeak2_left\tpeak2_right\tpeak2_open\tpeak2_closed\tpeak2_fraction\tFisher_test_p_val\tEmpirical_p-val\tMax_upper_empirical_p-val\tMax_lower_empirical_p-val\tNMI'
    outfile.write(outline + '\n')

    for chr in PeakDict.keys():
        PeakDict[chr].sort()
        for i in range(len(PeakDict[chr])-1):
            (RL1,RR1) = PeakDict[chr][i]
            for j in range(i+1,len(PeakDict[chr])):
                (RL2,RR2) = PeakDict[chr][j]
                print 'testing:', chr, RL1, RR1, RL2, RR2
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
                    LLs = []
                    NLLs = sum(float(L) > 0.5 for L in loglike)/(len(loglike) + 0.0)
                    RegionReads.append((NLLs,cgs,loglike))
                print 'found:', len(RegionReads), 'reads', 'needed:', minCov, 'reads'
                if len(RegionReads) < minCov:
                    continue
                if doSS:
                    emp_p_av = 0.0
                    max_possible_upper_emp_p_av = 0.0
                    max_possible_lower_emp_p_av = 0.0
                    emp_p_av = 0.0
                    NMI_av = 0.0
                    p_av = 0.0
                    CR1_av = 0.0
                    OR1_av = 0.0
                    CR2_av = 0.0
                    OR2_av = 0.0
                    for S in range(SS):
                        RegionReadsSubSampled = random.sample(RegionReads,minCov)
                        OpenOrClosed = []
                        for (NLLs,cgs,loglike) in RegionReadsSubSampled:
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
                            OpenOrClosed.append((NLLs,final_p1,final_p2))

#                        (Z1,Z2,Z3) = zip(*OpenOrClosed)
#                        print 'before sorting:', Z1
#                        OpenOrClosed.sort()
#                        (Z1,Z2,Z3) = zip(*OpenOrClosed)
#                        print 'after sorting:', Z1

                        # get empirical sampling distribution of coaccessibility values:
                        OpenOrClosed.sort()
                        ObservedMatches = []
                        STEP = len(OpenOrClosed)/QU
                        OMlist = []
                        for i in range(QU):
                            tempOMlist = []
                            for j in range(i*STEP,(i+1)*STEP):
                                tempOMlist.append(OpenOrClosed[j])
                            [Z,P,Q] = zip(*tempOMlist)
                            P = list(P)
                            Q = list(Q)
                            OMlist.append((P,Q))
                        for NS in range(Nsamp):
                            ObsMatchesInChunks = 0
                            for (P,Q) in OMlist:
                                Ps = np.array(P)
                                Qs = np.array(Q)
                                random.shuffle(Ps)
                                random.shuffle(Qs)
                                ObsMatchesInChunks += np.sum(Ps==Qs)
                            ObservedMatches.append(ObsMatchesInChunks)

                        # calculate normalized mutual information on coaccessibility:

                        [Z,P,Q] = zip(*OpenOrClosed)
                        OpenOrClosed = zip(P,Q)
                        P = list(P)
                        Q = list(Q)
                        P = np.array(P)
                        Q = np.array(Q)
                        NMI = NMIS(P,Q)

                        # calculate empirical p-values (assuming normal distribution of sampling coaccessibilities):
                        matches = np.sum(P==Q)
                        if matches < 0.5*len(P):
                            NMI = -NMI
                        ObservedMatches = np.array(ObservedMatches)
                        OMm = np.mean(ObservedMatches)
                        OMstd = np.std(ObservedMatches)
                        if OMstd == 0:
                            OMstd = 0.1
                        if matches < OMm:
                            emp_p_val = norm.cdf(matches,OMm,OMstd)
                            if emp_p_val == 0:
                                emp_p_val = -300
                            else:
                                emp_p_val = math.log10(emp_p_val)
                        else:
                            emp_p_val = 1 - norm.cdf(matches,OMm,OMstd)
                            if emp_p_val == 0:
                                emp_p_val = 300
                            else:
                                emp_p_val = -math.log10(emp_p_val)

                        # calculate maximum possible p-values:
                        max_possible_upper_emp_p = 1 - norm.cdf(len(P),OMm,OMstd)
                        if max_possible_upper_emp_p == 0:
                            max_possible_upper_emp_p = 300
                        else:
                            max_possible_upper_emp_p = -math.log10(max_possible_upper_emp_p)
                        max_possible_lower_emp_p = norm.cdf(0,OMm,OMstd)
                        if max_possible_lower_emp_p == 0:
                            max_possible_lower_emp_p = -300
                        else:
                            max_possible_lower_emp_p = math.log10(max_possible_lower_emp_p)

#                        print matches, len(P), OMm, OMstd, np.mean(P), np.mean(Q)
#                        print norm.cdf(len(P),OMm,OMstd), norm.cdf(0,OMm,OMstd)
#                        print norm.cdf(matches,OMm,OMstd), 1 - norm.cdf(matches,OMm,OMstd)

                        # calculate Fisher test p-value:
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
                        emp_p_av += emp_p_val
                        max_possible_upper_emp_p_av += max_possible_upper_emp_p 
                        max_possible_lower_emp_p_av += max_possible_lower_emp_p 
                        NMI_av += NMI
                        p_av += logp
                        CR1_av += CR1
                        OR1_av += OR1
                        CR2_av += CR2
                        OR2_av += OR2
                    (emp_p, max_possible_upper_emp_p, max_possible_lower_emp_p, NMI, pvalue, CR1, OR1, CR2, OR2) = (emp_p_av/SS, max_possible_upper_emp_p_av/SS, max_possible_lower_emp_p_av/SS,
                                                                                                                          NMI_av/SS, p_av/SS, CR1_av/SS, OR1_av/SS, CR2_av/SS, OR2_av/SS)

                    if str(emp_p) == 'nan':
                        print emp_p, max_possible_upper_emp_p, max_possible_lower_emp_p, NMI, pvalue, CR1, OR1, CR2, OR2

                    outline = chr + '\t' + str(RL1) + '\t' + str(RR1)
                    outline = outline + '\t' + str(OR1) + '\t' + str(CR1) + '\t' + str(OR1/(OR1 + CR1 + 0.0))
                    outline = outline + '\t' + str(RL2) + '\t' + str(RR2)
                    outline = outline + '\t' + str(OR2) + '\t' + str(CR2) + '\t' + str(OR2/(OR2 + CR2 + 0.0)) + '\t' + str(pvalue)
                    outline = outline + '\t' + str(emp_p)
                    outline = outline + '\t' + str(max_possible_upper_emp_p)
                    outline = outline + '\t' + str(max_possible_lower_emp_p)
                    outline = outline + '\t' + str(NMI)
                    outfile.write(outline + '\n')
                    print outline
#                else:
#                    OpenOrClosed = []
#                    for (cgs,loglike) in RegionReads:
#                        t = zip(cgs,loglike)
#                        RD = dict((int(x), float(y)) for x, y in t)
#                        (A,B) = (alph,bet)
#                        for pos in range(RL1,RR1):
#                            if RD.has_key(pos):
#                                p = RD[pos]
#                                Z = int(PSS*p)
#                                A = A + Z
#                                B = B + PSS - Z
#                        if beta.mean(A,B) > 0.5:
#                            final_p1 = 1
#                        else:
#                            final_p1 = 0
#                        (A,B) = (alph,bet)
#                        for pos in range(RL2,RR2):
#                            if RD.has_key(pos):
#                                p = RD[pos]
#                                Z = int(PSS*p)
#                                A = A + Z
#                                B = B + PSS - Z
#                        if beta.mean(A,B) > 0.5:
#                            final_p2 = 1
#                        else:
#                            final_p2 = 0
#                        OpenOrClosed.append((final_p1,final_p2))
#                    C00 = OpenOrClosed.count((0,0))
#                    C01 = OpenOrClosed.count((0,1))
#                    C10 = OpenOrClosed.count((1,0))
#                    C11 = OpenOrClosed.count((1,1))
#                    OR1 = C01 + C00
#                    CR1 = C10 + C11
#                    OR2 = C10 + C00
#                    CR2 = C01 + C11
#                    oddsratio, pvalue = fisher_exact([[C00, C01], [C10, C11]])
#                    outline = chr + '\t' + str(RL1) + '\t' + str(RR1)
#                    outline = outline + '\t' + str(OR1) + '\t' + str(CR1) + '\t' + str(OR1/(OR1 + CR1 + 0.0))
#                    outline = outline + '\t' + str(RL2) + '\t' + str(RR2)
#                    if pvalue == 0:
#                        pvalue = 1e-300
#                    outline = outline + '\t' + str(OR2) + '\t' + str(CR2) + '\t' + str(OR2/(OR2 + CR2 + 0.0)) + '\t' + str(-math.log10(pvalue))
#                    outfile.write(outline + '\n')

    outfile.close()
            
run()

