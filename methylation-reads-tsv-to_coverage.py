#################################
#                                #
# Last modified 2018/07/16       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import numpy as np
from scipy.stats import beta
from scipy.stats import binom
import random
import os
import math
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s methylation_reads_all.tsv loglike_threshold outfile [-stranded +|-] [-minAbsLogLike float] [-minAbsPValue float] [-BayesianIntegration window(bp) step alpha beta pseudosamplesize] [-N6mAweight pseudosamplesize genome.fa] [-saveNewSingleMoleculeFile filename]' % sys.argv[0]
        print '\tNote: the [-BayesianIntegration] option requires the [-minAbsPValue] option'
        print '\tNote: the [-saveNewSingleMoleculeFile] option requires the [-BayesianIntegration] option'
        print '\tNote: the [-N6mAweight] option only works together with the -BayesianIntegration option'
        sys.exit(1)

    reads = sys.argv[1]
    LLthreshold = float(sys.argv[2])
    outfilename = sys.argv[3]

    doStranded = False
    if '-stranded' in sys.argv:
        doStranded = True
        WantedStrand = sys.argv[sys.argv.index('-stranded') + 1]
        print 'will only output coverage from reads on the', WantedStrand, 'strand'

    minAbsLogLike = 0
    if '-minAbsLogLike' in sys.argv:
        minAbsLogLike = float(sys.argv[sys.argv.index('-minAbsLogLike') + 1])
        print 'will ignore bases with absolute loglikelihood values less than', minAbsLogLike

    doP = False
    minAbsPValue = 0
    if '-minAbsPValue' in sys.argv:
        doP = True
        minAbsPValue = float(sys.argv[sys.argv.index('-minAbsPValue') + 1])
        print 'will ignore bases with p-values higher than', minAbsPValue, 'and lower than', 1 - minAbsPValue

    doSaveNewFile = False
    doBI = False
    if '-BayesianIntegration' in sys.argv:
        if not doP:
            print 'data not specified to be in probability space, exiting'
            sys.exit(1)
        doBI = True
        window = int(sys.argv[sys.argv.index('-BayesianIntegration') + 1])
        step = int(sys.argv[sys.argv.index('-BayesianIntegration') + 2])
        alph = float(sys.argv[sys.argv.index('-BayesianIntegration') + 3])
        bet = float(sys.argv[sys.argv.index('-BayesianIntegration') + 4])
        PSS = int(sys.argv[sys.argv.index('-BayesianIntegration') + 5])
        print 'will integrate accessibility probabilities over windows of', window, 'bp in size, step size', step, 'bp, using (', alph, bet, ') as beta priors, and a pseudosample size of', PSS
        if '-saveNewSingleMoleculeFile' in sys.argv:
            doSaveNewFile = True
            NewFile = open(sys.argv[sys.argv.index('-saveNewSingleMoleculeFile') + 1],'w')
            print 'will save integrated basepair accessibilities probabilities into a new file:', sys.argv[sys.argv.index('-saveNewSingleMoleculeFile') + 1]
        doAweight = False
        if '-N6mAweight' in sys.argv:
            doAweight = True
            N6mAweight = int(sys.argv[sys.argv.index('-N6mAweight') + 1])
            genome_fasta = sys.argv[sys.argv.index('-N6mAweight') + 2]
            print 'will used different weight for N6mA', sys.argv[sys.argv.index('-saveNewSingleMoleculeFile') + 1]
            GenomeDict={}
            sequence=''
            inputdatafile = open(genome_fasta)
            for line in inputdatafile:
                if line[0]=='>':
                    if sequence != '':
                        GenomeDict[chr] = ''.join(sequence).upper()
                    chr = line.strip().split('>')[1]
                    sequence=[]
                    continue
                else:
                    sequence.append(line.strip())
            GenomeDict[chr] = ''.join(sequence).upper()

    CoverageDict = {}
    
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
        if line.startswith('chromosome\tstart\tend\tread_name'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 8:
            print 'skipping:', fields
            continue
        RN += 1
        if RN % 10000 == 0:
            print RN, 'lines processed'
        chr = fields[0]
        strand = fields[3]
        if doStranded:
            if strand != WantedStrand:
                 continue
        Ps = fields[6].split(',')
        LLs = fields[7].split(',')
        if CoverageDict.has_key(chr):
            pass
        else:
            CoverageDict[chr] = {}
        if doBI:
            PDict = {}
            for i in range(len(Ps)):
                pos = int(Ps[i])
                p = float(LLs[i])
                PDict[pos] = p
            positions = PDict.keys()
            minPos = min(positions)
            maxPos = max(positions)
            if doSaveNewFile:
                NewPos = []
                NewLLs = []
            for pos in range(minPos + window/2, maxPos - window/2, step):
                (A,B) = (alph,bet)
                for j in range(pos - window/2, pos + window/2):
                    if PDict.has_key(j):
                        p = PDict[j]
                        if doAweight:
                            if strand == '+' and GenomeDict[chr][j] == 'A':
                                Z = int(N6mAweight*p)
                                A = A + Z
                                B = B + N6mAweight - Z
                            elif strand == '-' and GenomeDict[chr][j] == 'T':
                                Z = int(N6mAweight*p)
                                A = A + Z
                                B = B + N6mAweight - Z
                            else:
                                Z = int(PSS*p)
                                A = A + Z
                                B = B + PSS - Z
                        else:
                            Z = int(PSS*p)
                            A = A + Z
                            B = B + PSS - Z
                final_p = beta.mean(A,B)
                newpos = pos - (pos % step)
                if doSaveNewFile:
                    NewPos.append(str(newpos))
                    NewLLs.append("{0:.2f}".format(final_p))
                if final_p > minAbsPValue and final_p < 1-minAbsPValue:
                    continue
                else:
                    pass
                if CoverageDict[chr].has_key(newpos):
                    pass
                else:
                    CoverageDict[chr][newpos] = [0,0]
                if final_p < LLthreshold:
                    CoverageDict[chr][newpos][1] += 1
                else:        
                    CoverageDict[chr][newpos][0] += 1
            if doSaveNewFile:
                outline = fields[0] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5]
                outline = outline + '\t' + ','.join(NewPos) + '\t' + ','.join(NewLLs)
                NewFile.write(outline + '\n')
        else:
            try:
                for i in range(len(Ps)):
                    pos = int(Ps[i])
                    ll = float(LLs[i])
                    if math.fabs(ll) >= minAbsLogLike:
                        pass
                    else:
                        continue
                    if doP:
                        if ll > minAbsPValue and ll < 1-minAbsPValue:
                            continue
                        else:
                            pass
                    if CoverageDict[chr].has_key(pos):
                        pass
                    else:
                        CoverageDict[chr][pos] = [0,0]
                    if ll < LLthreshold:
                        CoverageDict[chr][pos][1] += 1
                    else:        
                        CoverageDict[chr][pos][0] += 1
            except:
                print 'skipping read'
                print fields
                continue

    print 'finished inputting reads'

    if doSaveNewFile:
        NewFile.close()

    chromosomes = CoverageDict.keys()
    chromosomes.sort()

    outfile = open(outfilename,'w')

    outline = '#chr\tstart\tend\tmeth\tunmeth\tcov'
    outfile.write(outline + '\n')

    K=0
    for chr in chromosomes:
        print chr
        positions = CoverageDict[chr].keys()
        positions.sort()
        for pos in positions:
            outline = chr + '\t' + str(pos) + '\t' + str(pos + 1)
            outline = outline + '\t' + str(CoverageDict[chr][pos][1])
            outline = outline + '\t' + str(CoverageDict[chr][pos][0])
            outline = outline + '\t' + str(CoverageDict[chr][pos][0] + CoverageDict[chr][pos][1])
            outfile.write(outline + '\n')
            
    outfile.close()

            
run()

