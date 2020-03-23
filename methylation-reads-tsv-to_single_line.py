##################################
#                                #
# Last modified 2020/03/15       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import numpy as np
import random
import os
from sets import Set

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s methylation_reads_all.tsv outfile [-flipSign]' % sys.argv[0]
        sys.exit(1)

    reads = sys.argv[1]
    outfilename = sys.argv[2]

    doFS = False
    if '-flipSign' in sys.argv:
        doFS = True
        print 'will flip signs'

    ReadDict = {}
    
    if reads.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + reads
    elif reads.endswith('.gz'):
        cmd = 'gunzip -c ' + reads
    elif reads.endswith('.zip'):
        cmd = 'unzip -p ' + reads
    else:
        cmd = 'cat ' + reads
    RN = 0
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        if line == '':
            break
        if line.startswith('chromosome\tstart\tend\tread_name'):
            continue
        fields = line.strip().split('\t')
        RN += 1
        if RN % 100000 == 0:
            print str(RN/1000000.) +  'M lines processed'
        chr = fields[0]
        left = int(fields[1])
        right = int(fields[2])
        read = fields[3]
        loglike = fields[4]
        if doFS:
            loglike = str((-1)*float(loglike))
        if ReadDict.has_key(chr):
            pass
        else:
            ReadDict[chr] = {}
        if ReadDict[chr].has_key(read):
            pass
        else:
            ReadDict[chr][read] = {}
            ReadDict[chr][read]['ps'] = []
            ReadDict[chr][read]['lls'] = []
        for i in range(left,right+1):
            ReadDict[chr][read]['ps'].append(i)
            ReadDict[chr][read]['lls'].append(loglike)

    print 'finished inputting reads'

    chromosomes = ReadDict.keys()
    chromosomes.sort()

    outfile = open(outfilename,'w')

    K=0
    for chr in chromosomes:
        print chr
        for readID in ReadDict[chr].keys():
            K+=1
            if K % 100000 == 0:
                print K
            left = min(ReadDict[chr][readID]['ps'])
            right = max(ReadDict[chr][readID]['ps'])
            outline = chr + '\t' + str(left) + '\t' + str(right) + '\t' + '.' + '\t' + readID + '\t' + 'nan' + '\t'
            for p in ReadDict[chr][readID]['ps']:
                outline = outline + str(p) + ','
            outline = outline[0:-1] + '\t'
            for L in ReadDict[chr][readID]['lls']:
                outline = outline + str(L) + ','
            outline = outline[0:-1]
            outfile.write(outline + '\n')
            
    outfile.close()

            
run()

