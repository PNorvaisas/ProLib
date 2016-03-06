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
from voltool2 import *




# pdbf=open('/Users/Povilas/Projects/ProLib/AminoAcids.csv','r')
# rdr=csv.reader(pdbf, delimiter=',')
# aadata=[row for row in rdr]
# pdbf.close()
#
# aastr=','.join(['{} {}'.format(row[0],row[1]) for row in aadata[1:]])



aas="ALA A,ARG R,ASN N,ASP D,CYS C,GLN Q,GLU E,GLY G,HIS H,ILE I,LEU L,LYS K,MET M,PHE F,PRO P,SER S,THR T,TRP W,TYR Y,VAL V,ABA X,ASH D,CIR R,CME C,CMT C,CSD C,CSO C,CSW C,CSX C,CYM C,CYX C,DDE H,GLH G,HID H,HIE H,HIP H,HSD H,HSE H,HSP H,IAS D,KCX K,LYN K,MHO M,MLY K,MSE M,OCS C,PFF F,PTR Y,SEP S,TPO T"
aa_list=aas.split(',')
aa_list=[aa.split() for aa in aa_list]
aa={ a.strip().upper() : b.strip() for a,b in aa_list }

#Let's go to the directory where PDB structures are kept
os.chdir('/Users/povilas/Projects/ProLib/Hsp90aN')#
#Define file containing reference sequence for protein
ref_file='CAII-P00918_frag.fasta'
ref_file='Hsp90aN_seq.fasta'
ref_file='Crop_159.fasta'
ref=readseq(ref_file)
pdbfiles=glob.glob('*.pdb')


#u=MD.Universe('1fsn.pdb')

splitdir=dircheck('Split_to_fragments')

problems=[]

data=[]
matches=[]
header=['PDB ID','Segment','Segments','Salts_and_water','Heteroatoms','Length','Longest match','Longest match sequence']
salts=['HOH','ZN','SO4','NI','CO','CU','MN','HG','HGB','CMH','GOL','MBO','CL','BEZ','DMS','BCN','DMS','MES','BME']
data.append(header)
for ifile in pdbfiles:
	ipath, iname, itype=filename(ifile)
	#Create MDAnalysis atom universe object
	u=MD.Universe(ifile)
	#Iterate through segments (A,B...) to find ones containing protein sequences
	segids=list(set([segg.id for segg in u.segments]))
	print "{}, S: {} - {}/{}".format(ifile,len(segids),pdbfiles.index(ifile),len(pdbfiles))
	for segid in segids:
		seg=u.selectAtoms('segid {}'.format(segid))
		p=seg.selectAtoms('protein')
		np=seg.selectAtoms('not protein')
		#Check whether segment contains protein groups
		if len(p.resnames())>0:
			#print "Segment: {}, residues: {}".format(seg.name,len(p.resnames()))
			ress=''.join([ aa[res] for res in p.resnames() if res in aa.keys()])
			probs=[ res for res in p.resnames() if res not in aa.keys()]
			if len(probs)>0:
				problems.append(probs)
				print probs
			nres=sorted(list(set(np.resnames())))
			nresc=[r for r in nres if r not in salts]
			#Find longest matching sequence in the PDB
			lmatches=lcs(ref,ress)
			lmatch=[m for m in lmatches]
			longest=len(lmatch[0])
			matches.append(lmatch[0])
			#print 'Longest: {}'.format(longest)
			seg.atoms.write('{}/{}_{}.pdb'.format(splitdir,iname,segid))
			data.append([iname,segid,len(segids),';'.join(nres),';'.join(nresc),len(ress),longest,lmatch[0]])

writecsv(data,'Summary_sequences.csv',delim=',')


#os.chdir(splitdir)
#Now it's a good time to align and protonate fragments!
#Run AandP.py in fragments folder




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
crop_file='{}-{}_{}-best.fasta'.format(bk,len(bm),bm[0:4])
bestmatch=open(crop_file,'w')
bestmatch.write(bm)
bestmatch.close()


os.chdir('/Users/povilas/Projects/ProLib/CDK/CDK2')

#You can use other sequence file for cropping
#crop_file='CAII-P00918_frag.fasta'
crop_file='Crop_159.fasta'
cropseq=readseq(crop_file)
#cropseq=bm

pdbfiles=glob.glob('*.pdb')
pdbfiles=['1ol1.pdb']

#Name the summary file according the first 4 letters of the seuence and its length (useful to remember what was happening)
sfile='Summary_{}{}.csv'.format(cropseq[0:4],str(len(cropseq)))

results=[]
rhead=['PDB ID','Segment','Match length','Alignment','Cropped']
results.append(rhead)
results.append(['Ref','Ref',len(cropseq),cropseq,cropseq])

odir=dircheck('Cropped',dansw='c')
#Crop loop for each PDB structure.
# crop function has default options:
#keepmols=False - don't keep heteroatoms
#onlyfull=True - make crop only for completely matching structures
for ifile in pdbfiles:
    res=crop(ifile,cropseq,odir,[],[])
    results.extend(res)

writecsv(results,sfile,delim=',')
