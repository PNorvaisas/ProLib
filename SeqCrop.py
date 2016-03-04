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
	-i <file>    Name of input file (structure or text), or .pdb/.pqr to get all pdb/pqr files
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
Reference:  %(rfile)s
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
				raise Exception("No input file specified!")
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




	aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
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

	odir=dircheck(odir)
	ilist=genlist(ifile)
	ref=readseq(rfile)
	if '../' in odir:
		od='../'
	else:
		od=''
	sfile='{}Summary_{}{}.csv'.format(od,ref[0:4],str(len(ref)))

	data=[]
	header=['PDB ID','Segment','Match length','Alignment','Cropped']
	data.append(header)
	data.append(['Ref','Ref',len(ref),ref,ref])


	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]
	inps=ilist


	blist=slice_list(ilist,cores)

	for ifiles in blist:
		args.append((ifiles,ref,odir,q))

	result = p.map_async(crop_handler, args)
	start=time.time()
	prcprev=0

	while True:
		if result.ready():
			break
		else:
			#print q.qsize()
			#print q
			lcomp=q.qsize()
			prc = float(lcomp)*100/float(len(inps))

			if prc>prcprev:
				timepassed=float(time.time()-start)

				if prc!=0:
					total=int((timepassed*100/prc))
					remaining=int((timepassed*100/prc)-timepassed)
				else:
					total=int('inf')
					remaining=int('inf')
				print "Done {0:3d}%, {1:>5}/{2:<5} remaining: {3:<5} total: {4:<5}".format(int(prc),str(lcomp),str(len(inps)),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
				prcprev=prc
		time.sleep(2)
	print 'Collecting results....'
	results=result.get()
	#print results
	# #print len(results.keys())
	# itr=0
	itr=0
	for res in results:
		#print res[1]
		itr=itr+1
		prc=itr*100/len(results)
		print "{0:3d}%".format(int(prc))
		data.extend(res[1])
	writecsv(data,sfile,delim=',')


def crop_handler((ifiles,ref,odir,q)):
	data=[]
	matches=[]
	for ifile in ifiles:
		res=crop(ifile,ref,odir)
		q.put(ifile)
		data.extend(res)
	return ifiles, data

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