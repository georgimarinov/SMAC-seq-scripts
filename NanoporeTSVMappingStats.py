##################################
#                                #
# Last modified 2018/11/23       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os
import numpy as np
# from scipy import stats

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s methylation_reads_all.tsv(.bgz/.gz/.bz2/.zip) outfile_namee' % sys.argv[0]
        sys.exit(1)

    reads_file = sys.argv[1]
    outfilename = sys.argv[2]

    RLs = []
  
    KK = 0
    if reads_file.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + reads_file
    elif reads_file.endswith('.gz') or reads_file.endswith('.bgz'):
        cmd = 'gunzip -c ' + reads_file
    elif reads_file.endswith('.zip'):
        cmd = 'unzip -p ' + reads_file
    else:
        cmd = 'cat ' + reads_file
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        if line == '':
            break
        if line.startswith('#'):
            continue
        KK += 1
        if KK % 100000 == 0:
            print KK, 'reads processed'
        linefields = line.strip().split('\t')
        start = int(linefields[1])
        end = int(linefields[2])
        RL = end - start 
        RLs.append(RL)

    outfile = open(outfilename, 'w')

    R = np.array(RLs)
    Rmedian = np.median(R)

    outline = 'Total reads:\t' + str(len(RLs))
    outfile.write(outline + '\n')
    outline = 'Total bases:\t' + str(sum(RLs))
    outfile.write(outline + '\n')
    outline = 'Mean read length:\t' + str((sum(RLs) + 0.0)/len(RLs))
    outfile.write(outline + '\n')
    outline = 'Median read length:\t' + str(Rmedian)
    outfile.write(outline + '\n')
#    outline = 'Mode read length:\t' + str(Rmode)
#    outfile.write(outline + '\n')

    outfile.close()

run()

