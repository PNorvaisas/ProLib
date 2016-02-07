#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
#Get sequences of PDB structures

import platform
platform.python_version()



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




ifile=sys.argv[1]
summ=sys.argv[2]

ipath, iname, itype=filename(ifile)

u=MD.Universe(ifile)
chain=u.selectAtoms('segid A')
p=chain.selectAtoms('protein')
#print p.resnames

ress=''.join([ aa[res] for res in p.resnames() if res in aa.keys()])
sumff=open(summ,"a")
sumf=csv.writer(sumff, dialect='excel')
sumf.writerow([iname,len(ress),ress])
sumff.close()




