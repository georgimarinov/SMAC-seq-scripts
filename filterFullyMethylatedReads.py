##################################
#                                #
# Last modified 2018/07/25       # 
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
import re
import os
import math
from sets import Set

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def HasStretchOfFullyMethylatedBases(RD,W,F,step):

    positions = RD.keys()
    positions.sort()
    left = min(positions)
    right = max(positions)
    IsFullyMethylated = False
    for i in range (left,right-W,step):
        M = 0.0
        U = 0.0
        for j in range(i,i+W):
            if RD.has_key(j):
                p = RD[j]
                if p > 0.5:
                    M += 1
                else:
                    U += 1
                if M > (1-F)*W and U > (1-F)*W:
                    break
#        print i, i+W, M, U, M/(M+U), U/(M+U), F, M/(M+U) >= F
#        if M/(M+U) >= F:
        if U/(M+U) >= F:
            IsFullyMethylated = True
            break

    return IsFullyMethylated
        
def run():

    if len(sys.argv) < 3:
        print 'usage: python %s methylation_reads_all.tsv WindowSize minFraction [-keepShort] [-missingBasesFilter genome.fa basecontexts(comma-separated) minFraction [-doMBFSet]] ' % sys.argv[0]
        print '\Note: the script assumes Tombo 1.3 probabilities'
        print '\Note: the script will remove all reads for which a window of the indicated size is found with a fraction of methylated or unmethylated bases higher than the indicated'
        print '\Note: it will also remove all reads shorter than the WindowSize parameter unless the [-keepShort] option is used'
        print '\Note: it will also print to stdout by default'
        sys.exit(1)

    reads = sys.argv[1]
    W = int(sys.argv[2])
    F = float(sys.argv[3])

    step = int(((1-F)*W)/2)

    doNotKS = True
    if '-keepShort' in sys.argv:
        doNotKS = False

    doMBF = False
    if '-missingBasesFilter' in sys.argv:
        doMBF = True
        fasta = sys.argv[sys.argv.index('-missingBasesFilter') + 1]
        BaseContexts = sys.argv[sys.argv.index('-missingBasesFilter') + 2].split(',')
        MBFminFrac = float(sys.argv[sys.argv.index('-missingBasesFilter') + 3])
        GenomeDict={}
        sequence=''
        inputdatafile = open(fasta)
        for line in inputdatafile:
            if line[0]=='>':
                if sequence != '':
                    GenomeDict[chr] = ''.join(sequence).upper()
                chr = line.strip().split('>')[1]
                sequence=[]
                Keep=False
                continue
            else:
                sequence.append(line.strip())
        GenomeDict[chr] = ''.join(sequence).upper()
        doMBFSet = False
        if '-doMBFSet' in sys.argv:
            doMBFSet = True

    if reads.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + reads
    elif reads.endswith('.gz') or reads.endswith('.bgz'):
        cmd = 'zcat ' + reads
    elif reads.endswith('.zip'):
        cmd = 'unzip -p ' + reads
    else:
        cmd = 'cat ' + reads
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
        if len(fields) <= 6:
            continue
        cgs = fields[6].split(',')
        loglike = fields[7].split(',')
        t = zip(cgs,loglike)
        try:
            RD = dict((int(x), float(y)) for x, y in t)
        except:
#            print 'skipping read'
#            print line.strip()
            continue
        if doMBF:
            chr = fields[0]
            strand = fields[3]
            CoveredBases = len(RD.keys())
            left = min(RD.keys())
            right = max(RD.keys())
            BasesInRead = 0
            if doMBFSet:
                BasePositions = []
                if strand == '+':
                    for BC in BaseContexts:
                        BasePositions += [m.start() for m in re.finditer('(?=' + BC + ')', GenomeDict[chr][left:right+1])]
                        BasesInRead += GenomeDict[chr][left:right+1].count(BC)
                if strand == '-':
                    for BC in BaseContexts:
                        BasePositions += [m.start() for m in re.finditer('(?=' + getReverseComplement(BC) + ')', GenomeDict[chr][left:right+1])]
                        BasesInRead += GenomeDict[chr][left:right+1].count(getReverseComplement(BC))
                BasePositions = list(BasePositions)
                BasePositions = np.array(BasePositions)
                BasePositions = BasePositions + left
                BasePositions = list(BasePositions)
                BasePositions = Set(BasePositions)
                RSet = Set(RD.keys())
                RBS = BasePositions.intersection(RSet)
                if len(RBS)/(len(BasePositions) + 0.0) < MBFminFrac:
                    continue
                else:
                    pass
#                print chr, left, right, strand, CoveredBases, BasesInRead, CoveredBases/(BasesInRead + 0.0), len(BasePositions), len(RSet), len(RBS), len(RBS)/(len(BasePositions) + 0.0)
            else:
                if strand == '+':
                    for BC in BaseContexts:
                        BasesInRead += GenomeDict[chr][left:right+1].count(BC)
                if strand == '-':
                    for BC in BaseContexts:
                        BasesInRead += GenomeDict[chr][left:right+1].count(getReverseComplement(BC))
                if CoveredBases/(BasesInRead + 0.0) < MBFminFrac:
                    continue
                else:
                    pass
        left = int(fields[1])
        right = int(fields[2])
        if doNotKS:
            if right - left < W:
                continue
        if HasStretchOfFullyMethylatedBases(RD,W,F,step):
            continue
        else:
            print line
            
run()

