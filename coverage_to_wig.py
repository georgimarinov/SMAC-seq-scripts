##################################
#                                #
# Last modified 2019/05/04       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os
import math
from sets import Set

def run():

    if len(sys.argv) < 8:
        print 'usage: python %s coverage window step chrField MfieldID UfieldID chrom.sizes outprefix [-minCov N_reads]  [-CpGonly genome.fa]' % sys.argv[0]
        print '\Note: assumed format:'
        print '\t\t#chr\tstart\tend\tmeth\tunmeth\tcov'
        sys.exit(1)

    reads = sys.argv[1]
    W = int(sys.argv[2])
    step = int(sys.argv[3])
    chrFieldID = int(sys.argv[4])
    MFieldID = int(sys.argv[5])
    UFieldID = int(sys.argv[6])
    outprefix = sys.argv[8]

    chrominfo = sys.argv[7]
    chromInfoDict = {}
    linelist = open(chrominfo)
    for line in linelist:
        fields = line.strip().split('\t')
        chr = fields[0]
        chromInfoDict[chr] = int(fields[1])

    CoverageDict = {}

    doCpGonly = False
    if '-CpGonly' in sys.argv:
        doCpGonly = True
        print 'will only include Cs in CpG context'
        fasta = sys.argv[sys.argv.index('-CpGonly') + 1]
        GenomeDict={}
        sequence=''
        inputdatafile = open(fasta)
        for line in inputdatafile:
            if line[0]=='>':
                if sequence != '':
                    GenomeDict[chr] = ''.join(sequence).upper()
                chr = line.strip().split('>')[1]
                print chr
                sequence=[]
                Keep=False
                continue
            else:
                sequence.append(line.strip())
        GenomeDict[chr] = ''.join(sequence).upper()

        print 'finished inputting genomic sequence'

    MC = 1
    if '-minCov' in sys.argv:
        MC = int(sys.argv[sys.argv.index('-minCov') + 1])
        print 'will exclude regions with less than', MC, 'reads'
    
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
        RN += 1
        if RN % 1000000 == 0:
            print str(RN / 1000000) +  'M lines processed'
        chr = fields[chrFieldID]
        left = int(fields[chrFieldID + 1])
        right = int(fields[chrFieldID + 2])
        if left == right:
            right += 1
        if doCpGonly:
            if GenomeDict[chr][left-1:left+1] != 'CG' and GenomeDict[chr][left-2:left] != 'CG':
                continue
        M = int(fields[MFieldID])
        U = int(fields[UFieldID])
        if CoverageDict.has_key(chr):
            pass
        else:
            CoverageDict[chr] = {}
        for i in range(left,right):
            CoverageDict[chr][i] = (M,U)

    print 'finished inputting reads'

    chromosomes = CoverageDict.keys()
    chromosomes.sort()

    outfile1 = open(outprefix + '.meth.wig','w')
    outfile2 = open(outprefix + '.cov.wig','w')

    K=0
    for chr in chromosomes:
        print chr
        positions = CoverageDict[chr].keys()
        if len(positions) == 0:
            continue
        P1 = min(positions)
        PN = min(max(positions),chromInfoDict[chr]-1)
        for i in range ((P1/W)*W,min(chromInfoDict[chr],(PN/W)*W),step):
            M = 0
            U = 0
            C = []
            for pos in range(i, i + W):
                if CoverageDict[chr].has_key(pos):
                    (m,u) = CoverageDict[chr][pos]
                    M += m
                    U += u
                    C.append(m+u)
            if M+U >= MC and i + int(W/2.0) < chromInfoDict[chr]:
                S = (0.0 + M)/(M + U)
                outline1 = chr + '\t' + str(i + int(W/2.0)) + '\t' + str(min(chromInfoDict[chr],i + int(W/2.0) + step)) + '\t' + str(S)
                outline2 = chr + '\t' + str(i + int(W/2.0)) + '\t' + str(min(chromInfoDict[chr],i + int(W/2.0) + step)) + '\t' + str(sum(C)/(len(C) + 0.0))
                outfile1.write(outline1 + '\n')
                outfile2.write(outline2 + '\n')
            
    outfile1.close()
    outfile2.close()

            
run()

