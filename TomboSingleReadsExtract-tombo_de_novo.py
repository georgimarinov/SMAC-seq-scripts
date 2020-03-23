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
from sets import Set
import h5py
import numpy as np    
import os
import re

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s tombo.per_read_stats genome.fa outfile_prefix [-m5C-only] [-m6A-only] [-CG-only] [-CG-CG-only] [-GC-only] [-m6A-CG-only] [-m6A-GC-only] [-m6A-GC-CG-only] [-doT] [-T-only] [-generic bases(comma-separated)] [-excludeContext string(,string2,string3,...,stringN) radius] [-excludeChr chr1[,chr2,...,chrN]] [-chrPrefix string]' % sys.argv[0]
        print '\tnote: by default, As in all contexts and Cs in all contexts will be printed out'
        print '\tnote: the [-m6A-CG-only option] will print out As in all contexts and Cs in CpG context'
        print '\tnote: the [-m6A-GC-only option] will print out As in all contexts and Cs in GpC context'
        print '\tnote: the [-m6A-CG-GC-only option] will print out As in all contexts and Cs in GpC or CpG context'
        print '\tnote: the [-m5C-only] option will print out Cs in all contexts'
        print '\tnote: the [-CG-only] and [-GC-only] [-CG-GC-only options only apply if they [-m5C-only] has been specified'
        print '\tnote: the [-doT] option is off by defaultt; the [-T-only] option only applies if the [-doT] option has been enabled'
        print '\tnote: the [-generic] option can not be used in combination with any of the other options'
        sys.exit(1)

    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}

    do5C = True
    do6A = True
    CGonly = False
    GConly = False
    CGGConly = False

    doChrPrefix = False
    if '-chrPrefix' in sys.argv:
        doChrPrefix = True
        chrPrefix = sys.argv[sys.argv.index('-chrPrefix') + 1]

    ExcludedChrs = {}
    if '-excludeChr' in sys.argv:
        for chr in sys.argv[sys.argv.index('-excludeChr') + 1].split(','):
            ExcludedChrs[chr] = 1

    if '-m5C-only' in sys.argv:
        do6A = False
        GConly = False
        CGonly = False
        print 'will only output m5C positions'
        if '-CG-only' in sys.argv:
            CGonly = True
            print 'will only output m5C positions in CpG context'
        if '-GC-only' in sys.argv:
            GConly = True
            print 'will only output m5C positions in GpC context'
        if GConly and CGonly:
            print 'incompatible options, exiting'
            sys.exit(1)
        if '-CG-GC-only' in sys.argv:
            CGGConly = True
            print 'will only output m5C positions in GpC or CpG context'

    dom6AGCCGonly = False
    if '-m6A-GC-CG-only' in sys.argv:
        dom6AGCCGonly = True
        print 'will only output m6A, GpC and CpG positions'

    dom6AGConly = False
    if '-m6A-GC-only' in sys.argv:
        dom6AGConly = True
        print 'will only output m6A and GpC positions'

    dom6ACGonly = False
    if '-m6A-CG-only' in sys.argv:
        dom6ACGonly = True
        print 'will only output m6A and CpG positions'

    if '-m6A-only' in sys.argv:
        do5C = False
        print 'will only output m6A positions'

    doT = False
    if '-doT' in sys.argv:
        doT = True
        print 'will output T positions too'
        if '-T-only' in sys.argv:
            do5C = False
            do6A = False
            print 'will only output T positions'
        if CGonly or GConly or CGGConly or dom6AGCCGonly or dom6AGConly or '-m6A-only' in sys.argv:
            print 'incompatible options detected', exiting
            sys.exit(1)

    GENERIC = False
    if '-generic' in sys.argv:
        GENERIC = True
        if CGonly or GConly or CGGConly or dom6AGCCGonly or dom6AGConly or '-m6A-only' in sys.argv:
            print 'incompatible options detected', exiting
            sys.exit(1)
        genericDict = {}
        genericDict['+'] = {}
        genericDict['-'] = {}
        bases = sys.argv[sys.argv.index('-generic') + 1].split(',')
        for B in bases:
            genericDict['+'][B] = True
            genericDict['-'][DNA[B]] = True
        print genericDict
        
    tombo = sys.argv[1]
    fasta = sys.argv[2]
    outprefix = sys.argv[3]

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

    for chr in ExcludedChrs.keys():
        if GenomeDict.has_key(chr):
           del GenomeDict[chr]

    print 'finished inputting genomic sequence'

    doExcludeContext = False
    if '-excludeContext' in sys.argv:
        doExcludeContext = True
        contexts = []
        for cont in sys.argv[sys.argv.index('-excludeContext') + 1].split(','):
            contexts.append(cont.upper())
        radius = int(sys.argv[sys.argv.index('-excludeContext') + 2])
        print 'will filter out the following sequence contexts:', contexts, 'with a radius of', radius, 'bp around each match'

        for chr in GenomeDict.keys():
            print 'masking', chr
            ToBeMasked = {}
            for C in contexts:
                print C, chr, 
                for m in re.finditer(C,GenomeDict[chr]):
                    pos1 = m.start()
                    pos2 = m.end()
                    for i in range(pos1-radius,pos2+radius):
                        ToBeMasked[i] = 1
            print len(ToBeMasked), len(GenomeDict[chr]), len(ToBeMasked)/(len(GenomeDict[chr]) + 0.0)
            NewSeq = []
            for i in range(len(GenomeDict[chr])):
                if ToBeMasked.has_key(i):
                    NewSeq.append('X')
                else:
                    NewSeq.append(GenomeDict[chr][i])
            GenomeDict[chr] = ''.join(NewSeq)

        print 'finished masking genomic sequence'

    ReadDict = {}

    fMfile = h5py.File(tombo, 'r')

    groupC = fMfile['Statistic_Blocks']
    for block in groupC.itervalues():
#        print '..'
        reads =  block['read_ids']
        RIDtoReadNameDict = {}
        for rid in reads.attrs.iteritems():
            RIDtoReadNameDict[rid[1]] = rid[0]
        for BB in block.attrs.iteritems():
            if BB[0] == 'chrm':
                chr = BB[1].split('|')[0]
                if doChrPrefix:
                    chr = chr[len(chrPrefix):]
            if BB[0] == 'strand':
                strand = BB[1]
            if BB[0] == 'start':
                start = BB[1]
#        print chr, start
        if GenomeDict.has_key(chr):
            pass
        else:
            print 'chromosome not found in genome.fa file, skipping', chr
            continue
        for (pos,ll,ridNum) in block['block_stats']:
            rid = RIDtoReadNameDict[ridNum]
            if ReadDict.has_key((chr,strand,rid)):
                pass
            else:
                ReadDict[(chr,strand,rid)] = {}
#            seq = GenomeDict[chr][pos-1:pos+2]
            seq1 = GenomeDict[chr][pos:pos+2]
            seq2 = GenomeDict[chr][pos-1:pos+1]
            if strand == '+':
                if GENERIC:
                    try:
                        if genericDict['+'][seq1[0]]:
                            if ReadDict[(chr,strand,rid)].has_key(pos):
                                print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                            ReadDict[(chr,strand,rid)][pos] = ll
                        continue
                    except:
                        continue
                else:
                    pass
                if GConly:
                    if seq2 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGonly:
                    if seq1 == 'CG':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGGConly:
                    if seq1 == 'CG' or seq2 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGCCGonly:
                    if seq1 == 'CG' or seq2 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGConly:
                    if seq2 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6ACGonly:
                    if seq1 == 'CG':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                else:
                    if do5C and (seq1[0] == 'C'):
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif do6A and seq1[0] == 'A':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif doT and seq1[0] == 'T':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, 'T+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
            if strand == '-':
                if GENERIC:
                    try:
                       if genericDict['-'][seq1[0]]:
                            if ReadDict[(chr,strand,rid)].has_key(pos):
                                print 'positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1)
                            ReadDict[(chr,strand,rid)][pos] = ll
                       continue
                    except:
                        continue
                else:
                    pass
                if GConly:
                    if seq1 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGonly:
                    if seq2 == 'CG':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGGConly:
                    if seq2 == 'CG' or seq1 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGCCGonly:
                    if seq2 == 'CG' or seq1 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGConly:
                    if seq1 == 'GC':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6ACGonly:
                    if seq2 == 'CG':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                else:
                    if do5C and (seq1[0] == 'G'):
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '5C-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif do6A and seq1[0] == 'T':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif doT and seq1[0] == 'A':
                        if ReadDict[(chr,strand,rid)].has_key(pos):
                            print 'positions already seen', (chr,strand,rid), pos, 'T-, exiting', sys.exit(1)
                        ReadDict[(chr,strand,rid)][pos] = ll

    print 'finished processing file'

    reads = ReadDict.keys()
    reads.sort()

    outfile = open(outprefix + '.reads.tsv', 'w')

    for (chr,strand,rid) in reads:
        positions = ReadDict[(chr,strand,rid)].keys()
        positions.sort()
        if len(positions) == 0:
            print 'skipping:', (chr,strand,rid), ReadDict[(chr,strand,rid)]
            continue
        outline = chr + '\t' + str(positions[0]) + '\t' + str(positions[-1]) + '\t' + strand + '\t' + rid + '\t' + '.'
        Ps = ''
        LLs = ''
        for pos in positions:
            Ps = Ps + str(pos) + ','
            LLs = LLs + "{0:.2f}".format(ReadDict[(chr,strand,rid)][pos]) + ','
        outline = outline + '\t' + Ps[0:-1]
        outline = outline + '\t' + LLs[0:-1]
        outfile.write(outline + '\n')

    outfile.close()

    
run()
