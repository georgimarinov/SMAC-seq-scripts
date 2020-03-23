##################################
#                                #
# Last modified 2018/03/12       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import numpy as np

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s bismark.cov.gz chrom.sizes radius' % sys.argv[0]
        print '\tInput track can ge compressed'
        print '\tthe script will print to stdout by default'
        sys.exit(1)

    input = sys.argv[1]
    chromsizes = sys.argv[2]
    ChromSizeDict = {}
    linelist=open(chromsizes)
    for line in linelist:
        fields = line.strip().split('\t')
        chr = fields[0]
        end = int(fields[1])
        ChromSizeDict[chr] = end
    radius = int(sys.argv[3])

    if input.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + input
    elif input.endswith('.gz'):
        cmd = 'gunzip -c ' + input
    elif input.endswith('.zip'):
        cmd = 'unzip -p ' + input
    else:
        cmd = 'cat ' + input
    p = os.popen(cmd, "r")
    line = 'line'
    currentChr = ''
    currentWindow = 0
    currentWindowValues = []
    while line != '':
        line = p.readline().strip()
        if line.startswith('#'):
            continue
        if line.startswith('chrom\tstart\tend\tmeth\tunmeth\tcov'):
            continue
        if line == '':
            break
        fields = line.strip().split('\t')
        chr = fields[0]
        pos = int(float(fields[1]))
        meth = int(fields[3])
        unmeth = int(fields[4])
        cov = float(fields[5])
        M = meth/cov
        window = (pos/radius)*radius
        if currentChr == '':
            currentChr = chr
            currentWindow = window
            currentWindowValues = []
            currentWindowValues.append(M)
        elif currentChr != chr or window != currentWindow:
            outline = currentChr + '\t' + str(currentWindow) + '\t' + str(min(currentWindow + radius,ChromSizeDict[currentChr])) + '\t' + str(np.mean(currentWindowValues))
            print outline
            currentChr = chr
            currentWindow = window
            currentWindowValues = []
            currentWindowValues.append(M)
        else:
            currentWindowValues.append(M)

    outline = currentChr + '\t' + str(currentWindow) + '\t' + str(min(currentWindow + radius,ChromSizeDict[currentChr])) + '\t' + str(np.mean(currentWindowValues))
    print outline

run()
