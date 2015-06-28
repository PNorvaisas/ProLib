#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')


for mod in ['pip','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','Bio.SeqIO']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)		
		#pass # module doesn't exist, deal with it.

import Bio.SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import textwrap as tw
import itertools as IT
from os.path import basename

try:
	import MDAnalysis as MD
except ImportError, e:
	print "Module MDAnalysis not found"


aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
aa_list=aas.split(',')
aa_list=[aa.split() for aa in aa_list]
aa={ a.strip().upper() : b.strip() for a,b in aa_list }




def filename(ifile):
	base=os.path.basename(ifile)
	iname=os.path.splitext(base)[0]
	itype=os.path.splitext(base)[1].replace('.','')
	ipat=ifile.replace(base,'')
	return ipat, iname, itype

def dircheck(somedir):
	while True:   	
		if os.path.exists(somedir):
			qstn = "Directory %(somedir)s already exists! Delete, quit, continue or provide a new name (d/q/c/<type name>): " % vars()
			answ = raw_input(qstn)
			if answ == "d":
				shutil.rmtree(somedir)
				os.makedirs(somedir)
				break
			elif answ == "q":
				sys.exit("Have a nice day!")
			elif answ == "c":
				break
			else:
				somedir=answ
				continue
		else:
			os.makedirs(somedir)
			break
	return somedir

def align(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -10
	gap_extend = -0.5
	#alns = pairwise2.align.localxx(seq1, seq2) 
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend) 
	top_aln = alns[0]
	return top_aln

def alignlocal(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -20
	gap_extend = -20
	#alns = pairwise2.align.localxx(seq1, seq2) 
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend) 
	top_aln = alns[0]
	return top_aln

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
    return s1[xl-l: xl], xl-l, xl

def chunkstring(string, length):
    return [string[0+i:length+i] for i in range(0, len(string), length)]

def crop(ifile,mstring,aminoacids):
	ipath, iname, itype=filename(ifile)
	print iname
	u=MD.Universe(ifile)
	chain=u.selectAtoms('segid A')
	hetatm=u.selectAtoms('not protein')
	p=chain.selectAtoms('protein')
	ress=''.join([ aminoacids[res] for res in p.resnames() if res in aminoacids.keys()])
	print "Length of PDB chain: {}".format(len(ress))
	results=align(mstring,ress)
	#record = p.sequence(format='string')
	longest1=longest_substring(results[1], results[0])
	res2=alignlocal(ress,longest1[0])
	longest, start, stop=longest_substring(res2[1], res2[0])
	print longest, start, stop
	zerores=chain.residues.resids()[0]
	start=start+zerores
	stop=stop+zerores-1
	
	segment=chain.selectAtoms("resid {}:{}".format(start,stop))
	segg=''.join([ aminoacids[res] for res in segment.resnames() if res in aminoacids.keys()])
	
	aligned=''
	#print chunkstring(results[1],50)
	#print results[1]
	for a,b in IT.izip(chunkstring(results[1],50),chunkstring(results[0],50)):
		aligned=aligned+'PDB '+a+'\nRef '+b+'\n\n'
	print aligned
	print "Length of longest matching chain: {}\n".format(len(segment.resnames()))
	print "Longest matching chain:\n{}".format('\n'.join(tw.wrap(segg,50)))
	if len(hetatm)>0:
		chop=MD.Merge(segment,hetatm)
	else:
		chop=segment
	chop.atoms.write(iname+'_'+mstring[0:4]+str(len(mstring))+'_cropped.pdb')
	logf=open(iname+'_'+mstring[0:4]+str(len(mstring))+'_log.txt','w')
	log="Input PDB: {}\nReference amino acid chain:\n{}\n\nAlignment:\n{}\nLongest match: {}\nLongest matching chain:\n{}\n".format(ifile,'\n'.join(tw.wrap(mstring,50)),aligned,len(segment.resnames()),'\n'.join(tw.wrap(segg,50)))
	logf.write(log)
	logf.close()
	summary=mstring[0:4]+str(len(mstring))+'_summary.csv'
	if os.path.exists(summary):
		sumff=open(summary,"a")
		sumf=csv.writer(sumff, dialect='excel')
	else:
		sumff=open(summary,"a")
		sumf=csv.writer(sumff, dialect='excel')
		sumf.writerow(['Ref',len(mstring),mstring,mstring])
	sumf.writerow([iname,len(segment.resnames()),results[1],segg])
	sumff.close()


def readseq(ifile):
	data=open(ifile,'r')
	idata=data.read().split('\n')
	idata=[ line for line in idata if line!='' and not '|' in line]
	datastring=''.join(idata)
	data.close()
	return datastring


ifile=sys.argv[1]

matchingstr=readseq(sys.argv[2])
crop(ifile,matchingstr,aa)



