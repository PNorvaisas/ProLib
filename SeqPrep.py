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

from collections import Counter
#Needs Biopython
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


from Sequencing import *


aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
aa_list=aas.split(',')
aa_list=[aa.split() for aa in aa_list]
aa={ a.strip().upper() : b.strip() for a,b in aa_list }

#Let's go to the directory where PDB structures are kept
os.chdir('/Users/povilas/Projects/ProLib/CAII_raw')#
#Define file containing reference sequence for protein
ref_file='CAII-P00918_frag.fasta'
ref=readseq(ref_file)
pdbfiles=glob.glob('*.pdb')


data=[]
matches=[]
header=['PDB ID','Segment','Segments','Salts_and_water','Heteroatoms','Length','Longest match','Longest match sequence']
salts=['HOH','ZN','SO4','NI','CO','CU','MN','HG','HGB','CMH','GOL','MBO','CL','BEZ','DMS','BCN','DMS','MES','BME']
data.append(header)
for ifile in pdbfiles:
    print "{} - {}/{}".format(ifile,pdbfiles.index(ifile),len(pdbfiles))
    ipath, iname, itype=filename(ifile)
    #Create MDAnalysis atom universe object
    u=MD.Universe(ifile)
	#Iterate through segments (A,B...) to find ones containing protein sequences
    for seg in u.segments:
        p=seg.selectAtoms('protein')
        np=seg.selectAtoms('not protein')
        #Check whether segment contains protein groups
        if len(p.resnames())>0:
	        #print "Segment: {}, residues: {}".format(seg.name,len(p.resnames()))
			ress=''.join([ aa[res] for res in p.resnames() if res in aa.keys()])
	        nres=sorted(list(set(np.resnames())))
            nresc=[r for r in nres if r not in salts]
			#Find longest matching sequence in the PDB
			lmatches=lcs(ref,ress)
	        lmatch=[m for m in lmatches]
            longest=len(lmatch[0])
            matches.append(lmatch[0])
            #print 'Longest: {}'.format(longest)
			data.append([iname,seg.name,len(u.segments),';'.join(nres),';'.join(nresc),len(ress),longest,lmatch[0]])


writecsv(data,'Summary_sequences.csv',delim=',')

#Get summary on the matching sequences
stat=Counter(matches)
stat2={k:i for i,k in stat.iteritems()}

#Generate table for this summary
seqmatches=[]
seqmatches.append(['Occurences','Length','Sequence'])
print 'Occurences : Length - Sequence'
for i in sorted(stat2.keys(),reverse=True):
	print '{} : {} - {}'.format(i,len(stat2[i]),stat2[i])
	seqmatches.append([i,len(stat2[i]),stat2[i]])
writecsv(seqmatches,'Summary_matches.csv',delim=',')

#Select the sequence which occurs in most structures. This works great with CAII as it is, however ther might be cases where making sequence
#couple of amino acids shorter will increase the number of usable structures.
bk=sorted(stat2.keys(),reverse=True)[0]
bm=stat2[bk]
crop_file='{}-{}_{}.fasta'.format(bk,len(bm),bm[0:4])
bestmatch=open(crop_file,'w')
bestmatch.write(bm)
bestmatch.close()



#You can use other sequence file for cropping
#crop_file='CAII-P00918_frag.fasta'
#cropseq=readseq(crop_file)
cropseq=bm


#Name the summary file according the first 4 letters of the seuence and its length (useful to remember what was happening)
sfile='Summary_{}{}.csv'.format(cropseq[0:4],str(len(cropseq)))

results=[]
rhead=['PDB ID','Segment','Match length','Alignment','Cropped']
results.append(rhead)
results.append(['Ref','Ref',len(cropseq),cropseq,cropseq])

odir=dircheck('Cropped')
#Crop loop for each PDB structure.
# crop function has default options:
#keepmols=False - don't keep heteroatoms
#onlyfull=True - make crop only for completely matching structures
for ifile in pdbfiles:
    res=crop(ifile,cropseq,odir)
    results.extend(res)

writecsv(results,sfile,delim=',')
