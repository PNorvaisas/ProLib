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
from os.path import basename

from multiprocessing import Pool, Manager
import multiprocessing as mp


try:
	import itertools as IT
except ImportError, e:
	print "Module itertools not found"

try:
	import MDAnalysis as MD
except ImportError, e:
	print "Module MDAnalysis not found"

from Sequencing import *
from voltool2 import *



help_message = '''
SeqSplit.py
PDB analysis and split to segments

Flags:
	-h  Display this help message
	-v  Verbose output
Arguments:
	-i <file>    Name of input file (structure or text containing list of files), or .pdb to get all pdb files
	-r <file>    Reference sequence file
	-c <cores>   Number of cores to use
	-o <dir>     Output directory
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def usage():
	print help_message % vars()
	return

optionsset='''

Options:
<--------------------------------------------->
    Input:	%(ifile)s
Reference:	%(rfile)s
      Out:  %(odir)s
    Cores:  %(cores)s
<--------------------------------------------->
	'''



def main(argv=None):
	odir="Split_to_fragments"
	ifile=""
	#Script file
	rfile=""
	batches=5
	cores=2
	ncores=mp.cpu_count()
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:r:o:c:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return
			if option in ("-i", "--input"):
				ifile=value
			if option in ("-r", "--odir"):
				rfile=value
			if option in ("-o", "--odir"):
				odir=value
			if option in ("-c", "--cores"):
				#Set the number of cores to use
				cores=int(value)


		# for argument in args:
		# 	if argument in ("pqr", "--onlypqr"):
		# 		#Only generate pqr files
		# 		pqr = True
		# 		mcvol=False
		# 		doanalyze=False
		# 	if argument in ("mcvol", "--mcvol"):
		# 		#Perform colume calculation on previously generated pqr files
		# 		mcvol = True
		# 		pqr=False
		# 		doanalyze=False
		# 	if argument in ("analyze", "--analyze"):
		# 		#Analyse output of mcvol
		# 		doanalyze = True
		# 		mcvol=False
		# 		pqr=False
		# 		load=False



	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	#Check for the input integrity


	#Sanity check for input
	try:
		if ifile!="":
			ipath, iname, itype = filename(ifile)
			if iname == "" and itype=="" : #
				raise Exception("No input file specified.")
			elif iname=='' and itype!='':
				print "Working directory is expected to contain %(itype)s files!" % vars()
				pqrfl=glob.glob("./*.%(itype)s" % vars())
				pqrfl=sorted(pqrfl)
				if len(pqrfl)==0:
					raise Exception("No *.%(itype)s files found in working directory!" % vars())
				else:
					print "Found {} files in working directory!".format(len(pqrfl))
			elif not os.path.exists(ifile):
				raise Exception("Input file or directory does not exist.")
			elif os.path.isfile(ifile):
				if not itype in ("pdb"):
					print "File %(ifile)s is not supported" % vars()
				else:
					print "File %(ifile)s is expected to be a structure file" % vars()
		if not os.path.isfile(rfile):
			raise Exception("Reference sequence not found!")



	except Exception, e:
		print e
		#print "!!!-------------------------------------------!!!"
		#usage()
		sys.exit(1)




	aas="ALA A,ARG R,ASN N,ASP D,CYS C,GLN Q,GLU E,GLY G,HIS H,ILE I,LEU L,LYS K,MET M,PHE F,PRO P,SER S,THR T,TRP W,TYR Y,VAL V,ABA X,ASH D,CIR R,CME C,CMT C,CSD C,CSO C,CSW C,CSX C,CYM C,CYX C,DDE H,GLH G,HID H,HIE H,HIP H,HSD H,HSE H,HSP H,IAS D,KCX K,LYN K,MHO M,MLY K,MSE M,OCS C,PFF F,PTR Y,SEP S,TPO T"
	aa_list=aas.split(',')
	aa_list=[aa.split() for aa in aa_list]
	aa={ a.strip().upper() : b.strip() for a,b in aa_list }



	if cores>ncores:
		cores=ncores

	print optionsset %vars()

	print "Using {}/{} cores!".format(cores,ncores)
	#Let's go to the directory where PDB structures are kept
	#os.chdir('/Users/povilas/Projects/ProLib/CAII_raw')#
	#Define file containing reference sequence for protein
	#rfile='CAII-P00918_frag.fasta'

	ilist=genlist(ifile)
	ref=readseq(rfile)
	#ilist=glob.glob('*.pdb')



	#splitdir=dircheck('Split_to_fragments')
	splitdir=dircheck(odir)

	data=[]
	matches=[]
	header=['PDB ID','Segment','Segments','Salts_and_water','Heteroatoms','Length','Longest match','Longest match sequence']
	salts=['HOH','ZN','SO4','NI','CO','CU','MN','HG','HGB','CMH','GOL','MBO','CL','BEZ','DMS','BCN','DMS','MES','BME']
	data.append(header)

	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]
	inps=ilist


	blist=slice_list(ilist,cores)
	#list(chunks(ilist,batches))
	#print blist

	for ifiles in blist:
		args.append((ifiles,ref,aa,splitdir,salts,q))
	#print args
	#sys.exit(1)
	result = p.map_async(split_handler, args)
	start=time.time()
	prcprev=0
	comp=[]
	lcomp=0
	while True:
		if result.ready():
			break
		else:
			#print q.qsize()
			#print q
			lcomp=q.qsize()
			prc = float(lcomp)*100/float(len(inps))

			if prc>prcprev:
				#res=q.get()
				#comp.extend(res)
				timepassed=float(time.time()-start)
				#rem=[j for j in inps if j not in comp]
				if prc!=0:
					total=int((timepassed*100/prc))
					remaining=int((timepassed*100/prc)-timepassed)
				else:
					total=int('inf')
					remaining=int('inf')
				print "Done {0:3d}%, {1:>5}/{2:<5} remaining: {3:<5} total: {4:<5}".format(int(prc),str(lcomp),str(len(inps)),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
				#print rem
				# if len(rem)<5:
				# 	print 'Working on: {}'.format(rem)
				prcprev=prc
		time.sleep(2)
	print 'Collecting results....'
	results=result.get()
	#print results
	# #print len(results.keys())
	# itr=0
	data
	itr=0
	for res in results:
		itr=itr+1
		prc=itr*100/len(results)
		print "{0:3d}%".format(int(prc))
		data.extend(res[1])
		matches.extend(res[2])



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
	crop_file='{}{}-best.fasta'.format(bm[0:4],len(bm))
	bestmatch=open(crop_file,'w')
	bestmatch.write(bm)
	bestmatch.close()


def split_handler((ifiles,ref,aa,splitdir,salts,q)):
	data=[]
	matches=[]
	for ifile in ifiles:
		ipath, iname, itype=filename(ifile)
		#Create MDAnalysis atom universe object
		u=MD.Universe(ifile)
		#Iterate through segments (A,B...) to find ones containing protein sequences
		segids=list(set([segg.id for segg in u.segments]))
		#print "{}, S: {} - {}/{}".format(ifile,len(segids),ilist.index(ifile),len(ilist))
		for segid in segids:
			seg=u.selectAtoms('segid {}'.format(segid))
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
				seg.atoms.write('{}/{}_{}.pdb'.format(splitdir,iname,segid))
				data.append([iname,segid,len(segids),';'.join(nres),';'.join(nresc),len(ress),longest,lmatch[0]])
		q.put(ifile)
	return ifiles, data, matches

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def slice_list(input, size):
    input_size = len(input)
    slice_size = input_size / size
    remain = input_size % size
    result = []
    iterator = iter(input)
    for i in range(size):
        result.append([])
        for j in range(slice_size):
            result[i].append(iterator.next())
        if remain:
            result[i].append(iterator.next())
            remain -= 1
    return result


if __name__ == "__main__":
	sys.exit(main())