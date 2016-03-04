#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')

#Script to crop PDB structures according to given sequence


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

def writecsv(data,ofile,delim='\t'):
	f=open(ofile,'wb')
	ofile=csv.writer(f, delimiter=delim) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in data:
		#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()



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
	#Gap open and extend scores determine penalty for creating a gap in alignment
	gap_open = -10
	gap_extend = -0.5
	#alns = pairwise2.align.localxx(seq1, seq2)
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
	top_aln = alns[0]
	return top_aln

def alignlocal(seq1,seq2):
	matrix = matlist.blosum62
	#Gap open and extend scores determine penalty for creating a gap in alignment
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

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return lcs_set



def longmatch(seq,ref):
    matches={}
    lseq=[]
    i=0
    lstart=0
    lend=0
    mlen=0
    lk=['0:0']
    for rc in ref:
        if seq[0]==rc:
            lseq.append(rc)
            lend=i
        else:
            if len(lseq)>0:
                key='{}:{}'.format(lstart,lend)
                matches[key]=''.join(lseq)
                if len(lseq)>mlen:
                    lk=[key]
                    mlen=len(lseq)
                elif len(lseq)>mlen:
                    lk.append(key)
                #print '{} - {}'.format(key,len(lseq))
                lseq=[]
            lstart=i+1
        seq=seq[1:]
        i=i+1

    if len(lseq)>0:
        key='{}:{}'.format(lstart,lend)
        matches[key]=''.join(lseq)
        lk=[key]
    return matches,lk


def chunkstring(string, length):
    return [string[0+i:length+i] for i in range(0, len(string), length)]


def crop(ifile,mstring,odir,l,ph,keepmols=False,onlyfull=True):
	aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
	aa_list=aas.split(',')
	aa_list=[ aa.split() for aa in aa_list]
	aminoacids={ a.strip().upper() : b.strip() for a,b in aa_list }
	outdata=[]
	info=[]
	ipath, iname, itype=filename(ifile)
	#print iname
	u=MD.Universe(ifile)
	segids=list(set([segg.id for segg in u.segments]))
	for sname in segids:
		chain=u.selectAtoms('segid {}'.format(sname))
		#print "Segment: {}".format(sname)
		p=chain.selectAtoms('protein')
		if len(p.resnames())>0:
			ress=''.join([ aminoacids[res] if res in aminoacids.keys() else 'X'  for res in p.resnames()])
			hetatm=chain.selectAtoms('not protein')
			# First make alignment with the reference sequence
			results=align(mstring,ress)
			matches1,lk1=longmatch(results[1], results[0])
			longest1=matches1[lk1[0]]
			if onlyfull and len(longest1)!=len(mstring):
				print '{} - PDB chain shorter than reference - skipping!'.format(iname)
				continue

			start1,stop1=(int(nr) for nr in lk1[0].split(':'))
			#print 'Unadjusted: {}...{}'.format(longest1[:10],longest1[-10:])
			#print 'Unadjusted start: {} stop: {} length: {}'.format(start1,stop1,stop1-start1+1)

			res2=alignlocal(longest1,ress)
			matches,lk=longmatch(res2[1], res2[0])
			longest=matches[lk[0]]

			start,stop=(int(nr) for nr in lk[0].split(':'))
			#print 'Adjusted: {}...{}'.format(longest[:10],longest[-10:])
			#print 'Adjusted start: {} stop: {} length: {}'.format(start,stop,stop-start+1)
			#Adjust position of residues
			zerores=p.residues.resids()[0]
			#print 'Start in PDB: {}'.format(zerores)
			#print p.residues.resids()
			resseq=p.residues.resids()[start:stop+1]
			#start=start+zerores
			#stop=stop+zerores+1


			segment=p.selectAtoms("resid {}:{}".format(resseq[0],resseq[-1]))
			#Test
			segg=''.join([ aminoacids[res] for res in segment.resnames() ])#if res in aminoacids.keys()
			#print 'Final: {}...{}'.format(segg[:10],segg[-10:])
			#print 'Final start: {} stop: {} length: {}'.format(start,stop,stop-start)

			aligned=''
			#print chunkstring(results[1],50)
			#print results[1]
			for a,b in IT.izip(chunkstring(results[1],50),chunkstring(results[0],50)):
				aligned=aligned+'PDB '+a+'\nRef '+b+'\n\n'
			# print aligned
			#print "Length of longest matching chain: {}\n".format(len(segment.resnames()))
			#print "Longest matching chain:\n{}".format('\n'.join(tw.wrap(segg,50)))
			#print "Match length: {}".format(len(segg))
			#l for ligands
			if len(l)>0:
				lsel=[rid for rid in list(set(hetatm.resnames())) if rid in l]
			else:
				lsel=[]
			#ph for protein heteroatoms
			if len(ph)>0:
				psel=[rid for rid in list(set(hetatm.resnames())) if rid in ph]
			else:
				psel=[]
			#All heteroatoms
			hsel=lsel+psel

			if len(hsel)>0:
				hcrops=[hetatm.selectAtoms('resname {}'.format(hname)) for hname in hsel]
				#print '{} {}'.format(len(l),len(hcrops))
				if len(hcrops)==1:
					hatoms=hcrops[0]
				elif len(hcrops)>1:
					hcrop=MD.Merge(*hcrops)
					hatoms=hcrop
				else:
					hatoms=''
			else:
				hatoms=''

			if len(hatoms)>0:
				chop=MD.Merge(segment,hatoms)
			elif len(hetatm)>0 and keepmols:
				chop=MD.Merge(segment,hetatm)
			else:
				chop=segment

			chop.atoms.write('{}/{}.pdb'.format(odir,iname)) #,mstring[0:4],str(len(mstring))
			logf=open('{}/{}_{}{}_log.txt'.format(odir,iname,mstring[0:4],str(len(mstring))),'w')
			log="Input PDB: {}\nReference amino acid chain:\n{}\n\nAlignment:\n{}\nLongest match: {}\nLongest matching chain:\n{}\nPDB:\n{}\nReference:\n{}".format(ifile,'\n'.join(tw.wrap(mstring,50)),aligned,len(segment.resnames()),'\n'.join(tw.wrap(segg,50)),results[1],results[0])
			logf.write(log)
			logf.close()
			outdata.append([iname,'Ref',len(segment.resnames()),results[0],segg])
			outdata.append([iname,sname,len(segment.resnames()),results[1],segg])
			info.append([iname,sname,sname,';'.join(lsel),';'.join(psel)])

	return outdata, info



def readseq(ifile):
	data=open(ifile,'r')
	idata=data.read().split('\n')
	idata=[ line for line in idata if line!='' and not '|' in line]
	datastring=''.join(idata)
	data.close()
	return datastring



