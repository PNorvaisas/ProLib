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

try:
	import MDAnalysis as MD
except ImportError, e:
	print "Module MDAnalysis not found"


aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
aa_list=aas.split(',')
aa_list=[aa.split() for aa in aa_list]
aa={ a.strip().upper() : b.strip() for a,b in aa_list }




def filename(ifile):
	if ifile.split('.')[0]=='':
		ipat=''
		iname=''
		itype=ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep="\\"
		elif "/" in ifile.split('.')[0]:
			sep="/"
		else:
			ipat=''
			iname=ifile.split('.')[0]
			itype=ifile.split('.')[1]
			return ipat, iname, itype
		allpath=ifile.split('.')[0]	
		iname=allpath.split(sep)[-1]
		ipath=allpath.split(sep)[:-1]
		ipat='/'.join(ipath)
		itype=ifile.split('.')[1]
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

def longest_match(s1,s2):
	s1_frags=s1.strip('_').split('-')
	s2_frags=s2.strip('-').split('-')
	starts=[]
	stops=[]
	if len(s1_frags)==1 and len(s2_frags)==1:
		for i,s in enumerate(s2):
			#print i, s
			if i>0:
				if s!='-' and s2[i-1]=='-':
					starts.append(i)
				if s=='-' and s2[i-1]!='-':
					stops.append(i-1)
			if i==0 and s!='-':
				starts.append(i)
		seq=s1[starts[0]:stops[0]]	
	else:
		print 'Non-continuous align - skipping!'
		sys.exit(1)
	return seq, starts[0], stops[0]	

def crop(ifile,mstring,aminoacids):
	ipath, iname, itype=filename(ifile)
	u=MD.Universe(ifile)
	segments=list(set([s.name for s in u.segments]))
	print 'File: {}, segments: {}'.format(iname,segments)
	sg=u.selectAtoms('segid {}'.format(segments[0]))
	p=sg.selectAtoms('protein')
	het=sg.selectAtoms('not protein')
	ress=''.join([ aminoacids[res] for res in p.resnames() if res in aminoacids.keys()])
	print "Length of PDB chain: {}".format(len(ress))
	print "Length of ref chain: {}".format(len(mstring))
	results=align(ress,mstring)

	aligned=''
	for a,b in IT.izip(tw.wrap(results[0],50),tw.wrap(results[1],50)):
		aligned=aligned+'PDB '+a+'\nRef '+b+'\n\n'
	print aligned

	#record = p.sequence(format='string')
	longest, start, stop=longest_match(results[0], results[1])
	#print longest,start, stop
	shift=p.residues[0].id
	start=start+shift
	stop=stop+shift
	segment=sg.selectAtoms("resid {}:{}".format(start,stop))
	segg=''.join([ aminoacids[res] for res in segment.resnames() if res in aminoacids.keys()])
	

	print "Length of longest matching chain: {}\n".format(len(segment.resnames()))
	print "Longest matching chain:\n{}\n\n".format('\n'.join(tw.wrap(segg,50)))
	#sys.exit(1)
	new=MD.Merge(segment, het)
	new.atoms.write(iname+'_'+mstring[0:4]+str(len(mstring))+'_cropped.pdb')
	logf=open(iname+'_'+mstring[0:4]+str(len(mstring))+'_log.txt','w')
	log="Input PDB: {}\nReference amino acid chain:\n{}\n\nAlignment:\n{}\nLongest match: {}\nLongest matching chain:\n{}".format(ifile,'\n'.join(tw.wrap(mstring,50)),aligned,len(segment.resnames()),'\n'.join(tw.wrap(segg,50)))
	logf.write(log)
	logf.close()

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



