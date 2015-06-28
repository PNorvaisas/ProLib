#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import sys
import textwrap as tw
import itertools as IT
 


def readseq(ifile):
	data=open(ifile,'r')
	idata=data.read().split('\n')
	idata=[ line for line in idata if line!='' and not '|' in line]
	datastring=''.join(idata)
	data.close()
	return datastring

def longest_substring(s1, s2):
    t = [[0]*(1+len(s2)) for i in range(1+len(s1))]
    l, xl = 0, 0
    for x in range(1,1+len(s1)):
        for y in range(1,1+len(s2)):
            if s1[x-1] == s2[y-1]:
                t[x][y] = t[x-1][y-1] + 1
                if t[x][y]>l:
                    l = t[x][y]
                    xl  = x
            else:
                t[x][y] = 0
    return s1[xl-l: xl]


ifile1=sys.argv[1]
ifile2=sys.argv[2]
seq1=readseq(ifile1)
seq2=readseq(ifile2)

matrix = matlist.blosum62
gap_open = -10
gap_extend = -0.5

alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
 
top_aln = alns[0]
aln_1, aln_2, score, begin, end = top_aln
longest=longest_substring(aln_1,aln_2)


aligned=''
for a,b in IT.izip(tw.wrap(aln_1),tw.wrap(aln_2)):
	aligned=aligned+a+'\n'+b+'\n\n'
print aln_1,'\n',aln_2
output=open(sys.argv[1]+'_'+sys.argv[2]+'_alignment.txt','w')
out="{}\nRaw input 1: {}\n{}\nRaw input 2: {}\nScore: {}\nBegin: {}\nEnd: {}\n\n{}\n\nLongest match:\n{}".format(ifile1,'\n'.join(tw.wrap(seq1)),ifile2,'\n'.join(tw.wrap(seq2)),score,begin,end,aligned,'\n'.join(tw.wrap(longest)))
output.write(out)
output.close()

out_match=open(sys.argv[1]+'_'+sys.argv[2]+'_match.txt','w')
out_match.write(longest)
out_match.close()
