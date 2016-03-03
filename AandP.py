#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
"""
APChimera.py

Created by Povilas Norvaisas on 2016-03-02.
Copyright (c) 2016. All rights reserved.

"""




for mod in ['pip','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','time']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)
		#pass # module doesn't exist, deal with it.


from collections import OrderedDict
from collections import defaultdict
import numpy as np

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

help_message = '''
Use McVol in smart fashion.

Flags:
	-h  Display this help message
	-v  Verbose output
Arguments:
	-i <file>    Name of input file (structure or text), or .pdb/.pqr to get all pdb/pqr files
	-s <file>    Python script to execute in Chimera
	-b <batch>   Number of files for one instance of chimera
	-c <cores>   Number of cores to use

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
   Script:	%(sfile)s
  Batches:  %(batches)s
    Cores:  %(cores)s
<--------------------------------------------->
	'''



def main(argv=None):
	ifile=""
	#Script file
	sfile=''
	batches=3
	cores=2
	ncores=mp.cpu_count()

	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:b:s:c:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return
			if option in ("-i", "--input"):
				ifile=value
			if option in ("-b", "--input"):
				batches=int(value)
			if option in ("-s", "--script"):
				sfile=value
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
					print "Found {} files working directory!".format(len(pqrfl))
			elif not os.path.exists(ifile):
				raise Exception("Input file or directory does not exist.")
			elif os.path.isfile(ifile):
				if not itype in ("pqr", "pdb"):
					print "File %(ifile)s is not supported" % vars()
				else:
					print "File %(ifile)s is expected to be a structure file" % vars()
		if sfile!='' and os.path.isfile(sfile):
			spath, sname, stype = filename(sfile)
			print spath
			if spath!='':
				print 'Creating temporary script file in working directory....'
				stemp="{}_temp.{}".format(sname,stype)
				shutil.copyfile(sfile,stemp)
			else:
				stemp=sfile
		else:
			raise Exception("No script file specified.")



	except Exception, e:
		print e
		#print "!!!-------------------------------------------!!!"
		#usage()
		sys.exit(1)

	cores=2
	ncores=mp.cpu_count()
	#----------------------------
	if cores>ncores:
		cores=ncores

	print optionsset %vars()

	print "Using {}/{} cores!".format(cores,ncores)
	# os.chdir('/Users/Povilas/Projects/ProLib/CAII_raw/Aligned-Protonated')
	# sys.path.append('/Applications/Chimera.app/Contents/MacOS')
	# sfile='AP_script.py'
	# ifile='.pdb'
	ilist=genlist(ifile)


	odir=dircheck("Aligned-Protonated")
	#sys.exit(1)
	results, log=chimera_M(ilist,stemp,cores,batches)
	#print outputs
	if sfile!=stemp:
		os.remove(stemp)
		os.remove(stemp.replace('.py','.pyc'))

	lfile=open('Script_execution_log.txt','w')
	print 'Writing log file...'
	for item in log:
		lfile.write("%s\n" % item)
	lfile.close()





#-------------Functions------------


def runcmd(cmd):
	#Wrapper for execution of programs
	failure, output = commands.getstatusoutput(cmd)
	if failure:
		print '''Running failed \n %(cmd)s \n %(output)s'''.encode('utf-8') % vars();
	return failure, output



class NestedDict(dict):
	#Infinite map class which allows to organise data
	def __getitem__(self, key):
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())


def genlist(ifile):
	#Generates list of input files, checks their existance
	ilist=[]
	ipath, iname, itype=filename(ifile)
	if itype in ['pdb','pqr'] and os.path.isfile(ifile):
		ilist.append(ifile)
	elif itype in ['txt',''] and os.path.isfile(ifile):
		ifl=open(ifile,'r')
		idata=ifl.read().split('\n')
		idata=[fl.trim() for fl in idata if fl!='']
		for fld in idata:
			ilist.extend(genlist(fld))
	elif os.path.isdir(ifile):
		ilist.extend(glob.glob(ifile+"/*.pqr"))
	elif iname=='' and itype in ['pdb','pqr']:
		if itype in ['pdb','pqr']:
			ffiles=glob.glob('*.%(itype)s' % vars())
			#print ffiles
			ilist.extend(ffiles)
		elif itype=='txt':
			for tfile in glob.glob('*.%(itype)s' % vars()):
				ilist.extend(genlist(tfile))
		else:
			print "Bad file type %(inp)s!" % vars()
	return sorted(ilist)

def chimerahandler((ifiles,sfile,q)):
	#Multi-threaded execution McVol handler

	comm="chimera --nogui {} {}".format(' '.join(ifiles),sfile)
	fail, out=runcmd(comm)
	#time.sleep(2)
	#print ifiles
	#print fail, out
	if len(ifiles)>1:
		del ifiles[0]
	q.put(' '.join(ifiles))#'{}'.format(' '.join(ifiles))
	#, fail, out
	return ifiles, fail, out

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]



def chimera_M(ilist,sfile,cores,batches):
	#Multi-threaded Chimera execution
	#outputs=NestedDict()

	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]
	#All input
	inps=ilist

	ref=ilist[0]

	data=ilist[1:]

	blist=[[ref]]+[[ref]+bat for bat in list(chunks(data,batches))]


	for ifiles in blist:
		args.append((ifiles,sfile,q))
	#print args
	#sys.exit(1)
	result = p.map_async(chimerahandler, args)
	start=time.time()
	prcprev=0
	comp=[]
	while True:
		if result.ready():
			break
		else:
			#print q.qsize()
			if q.qsize()>0:
				res=q.get()
				#print 'Result: {}'.format(res)
				if ' ' in res:
					comp.extend(res.split(' '))
				else:
					comp.append(res)
				prc = float(len(comp))*100/float(len(inps))
				timepassed=float(time.time()-start)
				if prc>prcprev:
					#print completed
					rem=[j for j in inps if j not in comp]
					if prc!=0:
						total=int((timepassed*100/prc))
						remaining=int((timepassed*100/prc)-timepassed)
					else:
						total=int('inf')
						remaining=int('inf')
					print "Done {0:3d}%, {1:>5}/{2:<5} remaining: {3:<5} total: {4:<5}".format(int(prc),str(len(comp)),str(len(inps)),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
					#print comp
					#print rem
					if len(rem)<5:
						print 'Working on: {}'.format(rem)
					prcprev=prc
		time.sleep(2)
	print 'Collecting results....'
	results=result.get()
	#print results
	# #print len(results.keys())
	# itr=0
	log=[]
	itr=0
	for res in results:
		itr=itr+1
		prc=itr*100/len(results)
		print "{0:3d}%".format(int(prc))
		if not res[1]:
			#
			log.append(res[2])
		else:
			for ifl in res[0].split():
				print '{} has failed!'.format(ifl)



	return results, log



def filename(ifile):
	#Filename parser that separates path, name and extension
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
	#Checks for the presence of given directory and allows user to make selection what to do with it
	#Delete and overwrite, overwrite contents, change output directory
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





def install(package):
	#Installs packages if they are not present
	pip.main(['install', package])



#----------------------------------

if __name__ == "__main__":
	sys.exit(main())

