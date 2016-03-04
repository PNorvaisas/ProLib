#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
"""
voltool.py

Created by Povilas Norvaisas on 2014-04-11.
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
	-s <file>    Name of McVol settings file (if not given default settings will be used)
	-p <file>    PDB info file to use (contains information about the structure of PDB files and ligand choices)
				 Will be generated in the directory for each run, default name "PDB_info.txt"
	-c <cores>   Number of cores to use

Options:
	pqr	     Only convert selected files to pqr 
	mcvol	     Calculate volumes
	analyze	     Analyze McVol output and generate tables
	 

Settings:
Add McVol to your PATH!
Add obabel to your PATH!


!!!--------Requires following packages--------!!!
		     McVol
		   OpenBabel
		      pip
	     MDAnalysis Python module
!!!-------------------------------------------!!!
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
 Settings:	%(sfile)s
   PQR it:	%(pqr)s 
  Volumes:	%(mcvol)s
  Analyze:	%(doanalyze)s
    Cores:  %(cores)s
<--------------------------------------------->
	'''



def main(argv=None):
	#This is a wrapper containing main program workflow
	#At first I provide default variable values
	dset={'nmc':200, 'surfPT':2500, 'probe':1.3, 'membZmin':-100, 'membZmax':100, 'startgridspacing':1.0, 'membgridspacing':2.0, 'cavgridspacing':0.1, 'minVol':1, 'waterVol':18, 'DummyRad':1.7, 'blab':0, 'MembDim':4, 'CoreZMin':0, 'CoreZMax':0, 'CoreDim':0, 'CleftDim':3, 'CleftRel':70, 'CleftMethod':2, 'SurfaceCluster':1.5}
	#Settings template file
	settmp=''
	#Generate settings file from default values provided in dset map object
	for ds in dset.keys():
		#print ds, dset[ds]
		settmp=settmp+'{} {}\n'.format(ds,dset[ds])
	#Input file
	ifile=""
	#Settings file
	sfile=''
	rname=''
	cores=2
	#reference=False
	pqr=True
	doanalyze=True
	mcvol=True
	pdbinfofile="PDB_info.txt"
	pdbdata={}
	load=False

	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:s:p:c:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return
			if option in ("-i", "--input"):
				#Input file can be one pdb/pqr file, file containing list of pdb/pqr files or wildcard-type input .pdb/.pqr
				ifile=value
			if option in ("-s", "--settings"):
				#Enabled input of custom setting file - good for testing different parameters
				sfile=value
			if option in ("-p", "--pdbinfofile"):
				#User can manually provide PDB_info file
				pdbinfofile=value
			if option in ("-c", "--cores"):
				#Set the number of cores to use
				cores=int(value)


		for argument in args:
			if argument in ("pqr", "--onlypqr"):
				#Only generate pqr files
				pqr = True
				mcvol=False
				doanalyze=False
			if argument in ("mcvol", "--mcvol"):
				#Perform colume calculation on previously generated pqr files
				mcvol = True
				pqr=False
				doanalyze=False
			if argument in ("analyze", "--analyze"):
				#Analyse output of mcvol
				doanalyze = True
				mcvol=False
				pqr=False
				load=False



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
				print "Working directory is expected to contain files of %(itype)s!" % vars()
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

			elif os.path.isdir(ifile):
				print "%(ifile)s is expected to be a directory with pqr files!" % vars()
				pqrfl=glob.glob(ifile+"/*.pqr")
				if len(pqrfl)==0:
					raise Exception("No .pqr files found in directory %(ifile)s!")
				else:
					print "Found {} files in {} directory!".format(len(pqrfl),ifile)
		elif not load:
			raise Exception("No input file specified.")



	except Exception, e:
		print e
		#print "!!!-------------------------------------------!!!"
		#usage()
		sys.exit(1)


	#----------------------------


	print optionsset %vars()


	if sfile!='' and os.path.isfile(ifile):
		spath, sname, stype = filename(sfile)
		setf=open(sfile,'r')
		settmp=setf.read()
		setf.close()



	ilist=genlist(ifile)
	print "Input: {}\n".format(', '.join(ilist))





	#PDB_info file contains PDB file dissection information - which chains to use and which ligands to take
	#PDB_info file is automatically generated if user makes manual selection of chain and ligand for each PDB file
	#If PDB file was provided by user - it's not beeing generated in the current run.
	if os.path.isfile(pdbinfofile):
		ask='Do you want to use PDB info file {}? (Y/N) '.format(pdbinfofile)
		answ=raw_input(ask)
		if answ=='Y':
			pdbdata=readpdbdata(pdbinfofile)
		else:
			print 'PDB info file {} not used!'.format(pdbinfofile)
	else:
		answ='N'

	#--------------------------------
	# Workflow is divided into several main tasks:
	# 1. pqr - collect info on pdb file locations and convert them to pqr
	# 2. mcvol - perform volume calculation
	# 3. doanalyze - analyze output
	#To make script start at the specific point in the process corresponding flags can be provided with the input
	#By default program goes through all steps.
	#--------------------------------------
	if pqr:
		#Collect information on
		links, pdbdata=prepare(ilist,pdbdata,pqr)
		if answ=='N' and not os.path.isfile(pdbinfofile):
			writepdbdata(pdbdata)

	if mcvol:
		if not pqr:
			links, pdbdata=prepare(ilist,pdbdata,pqr)
		outputs=volumeit_M(links,settmp,cores)
		#sys.exit(1)
		saveout(outputs)

	if doanalyze:
		if not pqr and not mcvol:
			links, pdbdata=prepare(ilist,pdbdata,pqr)
			outputs=collect(links)

		#Analyze McVol outputs
		odata=analyze(outputs)
		#Join all cavities/clefts to one pqr file
		print 'Joining cavities or clefts from multiple files to one file'
		odata=functions_M(joincav_handler,odata,cores)

		#Find centers for all cavities/clefts
		print 'Finding centers for all cavities/clefts'
		odata = functions_M(findcenters_handler,odata,cores)#findcenters(odata)
		#Pickle center data for later use
		f = open('Odata_all.pckl', 'w')
		pickle.dump(odata, f)
		f.close()
		#Generate summary tables
		tables=summarize(odata,pdbdata)
		#Write tables
		writesheets(tables)

		print 'Results written!'






#-------------Functions------------
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

def summarize(odata,pdbdata):
	#Makes tables out of the map data structure
	tables=NestedDict()
	summary=[]
	results=[]

	sheader=['File','Fragment', 'Protein segment','Ligand segment','Ligand name','Total surface, A^2','Protein volume, A^3','Solvent accessible volume, A^3', 'VdW volume, A^3','Cavity volume, A^3','Cleft volume, A^3','Waters placed']
	rheader=['File','Fragment','Cavity number','Type','Volume','Waters placed','x','y','z','Link']
	arkeys=['Number','Type','Volume','Waters','Center','Link']

	results.append(rheader)
	summary.append(sheader)
	for k in sorted(odata.keys()):
		frag=odata[k]
		#results=[]

		for f in sorted(frag.keys()):
			general=frag[f]['General']
			pdb=pdbdata[k]
			line=[k,f]

			for c in ['Protein segment','Ligand segment','Ligand name','Surface','Volume','SAV','VdW','Cavity volume', 'Cleft volume','Total waters']:
				if c in general.keys():
					line.append(general[c])
				elif c in pdb.keys():
					line.append(pdb[c])
				else:
					line.append('')
			summary.append(line)
			cavities=frag[f]['Volumes']

			if len(cavities.keys())!=0:
				for cav in sorted(cavities.keys()):
					cavity=cavities[cav]
					rline=[k,f]
					rline=getcavdat(cavity,f,rline,arkeys)
					results.append(rline)

	tables['Volumes']=results
	tables['Summary']=summary
	return tables
#
#
# def summalign(odata,results):
# 	aligntable=[]
# 	header=['Reference','R Frag','Target','T Frag','Alignment','R index','R Type','R Volume','R Waters placed','Rx','Ry','Rz','T index', 'T Type','T Volume','T Waters placed','Tx','Ty','Tz']
# 	aligntable.append(header)
# 	vkeys=['Type','Volume','Waters','Center']
# 	for k in results.keys():
# 		ref,tgt=k.split('_')
# 		for f in results[k].keys():
# 			ref_f,tgt_f=f.split('_')
# 			rdat=odata[ref][ref_f]['Volumes']
# 			tdat=odata[tgt][tgt_f]['Volumes']
# 			for al in results[k][f]:
# 				#print k,f,al,al[0],al[1],al[2]
# 				alin=al[0]
# 				rin=al[1]
# 				tin=al[2]
# 				rline=[ref,ref_f,tgt,tgt_f,alin]
# 				if rin!='':
# 					rline.append(rin)
# 					rline=getcavdat2(rdat,rin,rline,vkeys)
# 				else:
# 					rline.extend(['','','','','','',''])
# 				if tin!='':
# 					rline.append(tin)
# 					rline=getcavdat2(tdat,tin,rline,vkeys)
# 				else:
# 					rline.extend(['','','','','','',''])
#
# 				aligntable.append(rline)
#
# 	return aligntable
#
#
def getcavdat(dataset, fragment, line, dkeys):
	#Generates links to McVol results
	for c in dkeys:
		if c in dataset.keys():
			if c=='Link':
				link=dataset[c]
				link="{}/{}".format(fragment,link)
				line.append(link)
			elif c=='Center':
				cnt=dataset[c]
				line.extend([dataset[c][0],dataset[c][1],dataset[c][2]])
			else:
				line.append(dataset[c])
		else:
			line.append('')
	return line
#
# def getcavdat2(dataset,c,line, dkeys):
# 	#Generates links to McVol results
# 	if c in dataset.keys():
# 		vol=dataset[c]
# 		for d in dkeys:
# 			if d=='Center':
# 				cnt=vol[d]
# 				line.extend([cnt[0],cnt[1],cnt[2]])
# 			else:
# 				line.append(vol[d])
# 	else:
# 		line.extend(['','','','','',''])
#
# 	return line
#

def joincav_handler((odata_frag,q)):
	#frags=odata[k]
	#odata=NestedDict()
	for k in odata_frag.keys():
		frags=odata_frag[k]
		for f in frags.keys():
			cavdata=frags[f]['Volumes']#odata[k][f]['Volumes']
			cavname='{}/{}_cavities.pqr'.format(k,k+'-'+f)
			clfname='{}/{}_clefts.pqr'.format(k,k+'-'+f)
			cavities=open(cavname,'w')
			clefts=open(clfname,'w')
			ica=1
			icl=1

			for cav in sorted(cavdata.keys()):
				tp=cavdata[cav]['Type']
				link=cavdata[cav]['Link']
				link="{}/{}/{}".format(k,f,link)
				inp=open(link,'r')
				idata=inp.read()
				inp.close()
				for row in idata.split('\n'):
					data=[it.strip() for it in row.split()]
					data=[numerize(it) for it in data]
					if len(data)==0 or 'END' in data:
						continue
					else:
						if tp=='cavity':
							data[1]=int(ica)
						elif tp=='cleft':
							data[1]=int(icl)
						#print data
						datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}'.format(*data)
						if tp=='cavity':
							cavities.write(datastr+'\n')
							ica=ica+1
						elif tp=='cleft':
							clefts.write(datastr+'\n')
							icl=icl+1
			cavities.close()
			clefts.close()
			if ica==1:
				os.remove(cavname)
			else:
				odata_frag[k][f]['General']['Cavity link']=cavname
			if icl==1:
				os.remove(clfname)
			else:
				odata_frag[k][f]['General']['Cleft link']=clfname
		q.put(k)
		#print type(odata_frag)
	return odata_frag

def mapmerge(mainmap,map1):
	"""
	Deep merge of maps
	"""
	for key, value in map1.items():
		#print key
		if key in mainmap.keys() and isinstance(value,dict):
			print type(value)
			mainmap[key]=mapmerge(mainmap[key],value)
		else:
			mainmap[key] = value

	return mainmap






def joincav(odata):
	#Joins cavities or clefts from multiple files to one file
	print 'Joining cavities and clefts to single files!'
	oks=odata.keys()
	for k in odata.keys():
		prc=float(oks.index(k)+1)*100/float(len(oks))
		print 'Working on: {0} {1:3}/{2:3} {3:3d}%'.format(k,oks.index(k),len(oks),int(prc))
		frags=odata[k]
		for f in frags.keys():
			cavdata=frags[f]['Volumes']#odata[k][f]['Volumes']
			cavname='{}/{}_cavities.pqr'.format(k,k+'-'+f)
			clfname='{}/{}_clefts.pqr'.format(k,k+'-'+f)
			cavities=open(cavname,'w')
			clefts=open(clfname,'w')
			ica=1
			icl=1

			for cav in sorted(cavdata.keys()):
				tp=cavdata[cav]['Type']
				link=cavdata[cav]['Link']
				link="{}/{}/{}".format(k,f,link)
				inp=open(link,'r')
				idata=inp.read()
				inp.close()
				for row in idata.split('\n'):
					data=[it.strip() for it in row.split()]
					data=[numerize(it) for it in data]
					if len(data)==0 or 'END' in data:
						continue
					else:
						if tp=='cavity':
							data[1]=int(ica)
						elif tp=='cleft':
							data[1]=int(icl)
						#print data
						datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}'.format(*data)
						if tp=='cavity':
							cavities.write(datastr+'\n')
							ica=ica+1
						elif tp=='cleft':
							clefts.write(datastr+'\n')
							icl=icl+1
			cavities.close()
			clefts.close()
			if ica==1:
				os.remove(cavname)
			else:
				odata[k][f]['General']['Cavity link']=cavname
			if icl==1:
				os.remove(clfname)
			else:
				odata[k][f]['General']['Cleft link']=clfname

	return odata




def functions_M(func,odata,cores):
	#print "\nFinding center coordinates for cavities and clefts!\n"
	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]
	inps=odata.keys()

	slist=slice_list(inps,cores)
	#print slist
	for slice in slist:
		odata_frag={k:odata[k] for k in slice}
		args.append((odata_frag,q))
	#print args
	#sys.exit(1)
	result = p.map_async(func, args)
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
	itr=0
	for res in results:
		#print res
		itr=itr+1
		prc=itr*100/len(results)
		print "{0:3d}%".format(int(prc))
		odata=mapmerge(odata,res)
	return odata


def findcenters_handler((odata_frag,q)):
	for k in odata_frag.keys():
		frag=odata_frag[k]
		for f in frag.keys():
			vols=frag[f]['Volumes']
			for v in vols.keys():
				linkt=vols[v]['Link']
				#print rlinkt
				link='{}/{}/'.format(k,f)+linkt
				#print rlink
				coords=centerpqr(link)
				vols[v]['Center']=coords
		q.put(k)
	return odata_frag





def findcenters(odata):
	#Finds center coordinates for cavities and clefts
	print "\nFinding center coordinates for cavities and clefts!\n"
	for k in odata.keys():
		frag=odata[k]
		for f in frag.keys():
			vols=frag[f]['Volumes']
			for v in vols.keys():
				linkt=vols[v]['Link']
				#print rlinkt				
				link='{}/{}/'.format(k,f)+linkt
				#print rlink
				coords=centerpqr(link)
				vols[v]['Center']=coords
	return odata






def centerpqr(pqrfile):
	#Finds the center of body described in PQR file. Approximates it to the nearest point.
	table=tableout(pqrfile)
	ttable=zip(*table)
	#print ttable
	xc=ttable[5]
	yc=ttable[6]
	zc=ttable[7]
	#print 'Min x:', min(xc)
	#print 'Min y:', min(yc)
	#print 'Min z:', min(zc)	
	xcent=sum(xc)/float(len(xc))
	ycent=sum(yc)/float(len(yc))
	zcent=sum(zc)/float(len(zc))
	#xdif=[abs(x-xcent) for x in xc]
	#ydif=[abs(y-ycent) for y in yc]
	#zdif=[abs(z-zcent) for z in zc]
	distances=[math.sqrt(math.pow(x-xcent,2)+math.pow(y-ycent,2)+math.pow(z-zcent,2)) for x,y,z in IT.izip(xc,yc,zc)]
	mind=distances.index(min(distances))
	#print mind
	#print xc[mind]
	coords=(xc[mind],yc[mind],zc[mind])
	#print "Center: ({},{},{})".format(xcent,ycent,zcent) 
	#print "Chosen:", coords

	return coords

def volsave(pqr,loc,ind,rad):
	pqrf=open(loc,'w')
	count=0
	for x in pqr.keys():
		ys=pqr[x]
		for y in ys.keys():
			for z in ys[y]:
				count=count+1
				data=['ATOM',count,'DUM','DUM',ind,x,y,z,ind,rad]
				datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}\n'.format(*data)
				pqrf.write(datastr)

	pqrf.close()
	return count


def inpqr(link,coords,dist=1):
	#Checks whether the given coordinates exist in PQR file
	data=tableout(link)
	ttable=zip(*data)
	xc=ttable[5]
	yc=ttable[6]
	zc=ttable[7]
	distances=[math.sqrt(math.pow(x-coords[0],2)+math.pow(y-coords[1],2)+math.pow(z-coords[2],2)) for x,y,z in IT.izip(xc,yc,zc)]
	#print 'Distance', round_to(min(distances),0.01)
	#Maximal distance between center and any point is 1A
	if min(distances) <= dist:
		mind=distances.index(min(distances))
		return True
	else:
		return False



def round_to(n, precission):
	#Round a number to desired precision
	correction = 0.5 if n >= 0 else -0.5
	return int(n/precission+correction)*precission

def tableout(inp):
	#Read file to a list
	ifl=open(inp, 'r')
	idata=ifl.read()
	ifl.close()
	table=[]
	for row in idata.split('\n'):
		data=[it.strip() for it in row.split()]
		data=[numerize(it) for it in data]
		if len(data)!=0:
			table.append(data)

	return table




def volumeit(links,settings):
	#McVol execution - single threaded
	outputs=NestedDict()

	for k in links.keys():
		for f in links[k]:
			if f!='' and os.path.isfile(f):
				fpath, fname, ftype=filename(f)
				if '_' in fname:
					part=fname.split('_')[1]
				else:
					part='vol'
				wdir=os.path.join(k,part)
				if os.path.isdir(wdir):
					shutil.rmtree(wdir)
				os.makedirs(wdir)
				sfl=fpath+'/'+fname+'.setup'
				sfile=open(sfl,'w')
				sfile.write(settings)
				sfile.close()
				cwwd=os.getcwd()
				print "Calculating: %(fname)s" % vars()
				os.chdir(wdir)
				comm="McVol ../%(fname)s" % vars()
				fail, out=runcmd(comm)
				print fail, out
				os.chdir(cwwd)
				os.remove(sfl)
				if not fail:
					outputs[k][part]=out
				else:
					print out
	return outputs

def mcvolhandler((k,f,settings,q)):
	#Multi-threaded execution McVol handler
	if f!='' and os.path.isfile(f):
		fpath, fname, ftype=filename(f)
		if '_' in fname:
			part=fname.split('_')[-1]
		else:
			part='vol'
		wdir=os.path.join(k,part)
		if os.path.isdir(wdir):
			shutil.rmtree(wdir)
		os.makedirs(wdir)
		sfl=fpath+'/'+fname+'.setup'
		sfile=open(sfl,'w')
		sfile.write(settings)
		sfile.close()
		cwwd=os.getcwd()
		print "Calculating: %(fname)s" % vars()
		os.chdir(wdir)
		comm="McVol ../%(fname)s" % vars()
		fail, out=runcmd(comm)
		#print fail, out
		os.chdir(cwwd)
		os.remove(sfl)

	q.put('{}_{}'.format(k,part))
	return k, part, fail, out



def volumeit_M(links,settings,cores):
	#Multi-threaded McVol execution
	outputs=NestedDict()

	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]
	inps=[]

	for k in links.keys():
		for f in links[k]:
			args.append((k,f,settings,q))
			inps.append('{}_{}'.format(k,f))
	#print args
	result = p.map_async(mcvolhandler, args)
	start=time.time()
	prcprev=0
	comp=[]
	while True:
		if result.ready():
			break
		else:
			if q.qsize()>0:
				res=q.get()
				if isinstance(res,list):
					comp.extend(res)
				else:
					comp.append(res)
				prc = float(len(comp)+1)*100/float(len(args))
				timepassed=float(time.time()-start)
				if prc>prcprev:
					#print completed
					rem=[j for j in inps if j not in comp ]
					if prc!=0:
						total=int((timepassed*100/prc))
						remaining=int((timepassed*100/prc)-timepassed)
					else:
						total=int('inf')
						remaining=int('inf')
					print "Done {0:3d}%, {1:>5}/{2:<5} remaining: {3:<5} total: {4:<5}".format(int(prc),str(len(comp)),str(len(rem)),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
					if len(rem)<5:
						print 'Working on: {}'.format(rem)
					prcprev=prc
		time.sleep(2)
	print 'Collecting results....'
	results=result.get()
	#print len(results.keys())
	itr=0
	for res in results:
		itr=itr+1
		prc=itr*100/len(results)
		print "{0:3d}%".format(int(prc))
		if not res[2]:
			outputs[res[0]][res[1]]=res[3]
		else:
			print '{}, {} has failed!'.format(res[0],res[1])


	return outputs

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

def saveout(outputs):
	#Saves McVol shell output
	print "Saving outputs..."
	for k in outputs.keys():
		newmap=outputs[k]
		for p in newmap.keys():
			ofile=open('%(k)s/%(k)s_%(p)s.txt' % vars(), "w")
			ofile.write(outputs[k][p])
			ofile.close()

def writesheets(sheets):
	#Writes organized data from sheets object to file.
	#odir=dircheck('Split')
	for i in sheets.keys():
		if i in ['Summary','Volumes','Alignment']:
			oname=i+".csv"
		else:
			print 'Unknown key for table: {}'.format(i)

		ofile=csv.writer(open(oname,"wb"), dialect='excel') #,delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
		for row in sheets[i]:
			row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
			ofile.writerow(row)
		#ofile.close()


def prepare(ilist, pdbdata,pqr):
	#Generates list of actual files from the input string
	if len(pdbdata.keys())==0:
		pdbdata=NestedDict()
	links={}
	for f in ilist:
		if os.path.isfile(f):
			fpath, fname, ftype=filename(f)
			if ftype=='pdb':
				if not pqr and os.path.exists(fname):
					links[fname]=[it for it in glob.glob('%(fname)s/%(fname)s*.pqr' % vars()) if not 'cavities' in it and not 'clefts' in it]
				else:
					odir=dircheck(fname)
					pdb_files, pdbdata=splitpdb(f,pdbdata)
					pqrfiles=[]
					for p in pdb_files:
						try:
							pqrfile=pqrit(p)
						except:
							print 'Cannot convert {} to PQR! Skipping!'.format(f)
						pqrfiles.append(pqrfile)
					links[fname]=pqrfiles
			elif ftype=='pqr':
				#if "_" in fname:
				#part not finished!!!
				if fname.count('_')>0:
					fname='_'.join(fname.split('_')[:-1])
				if fname in links.keys():
					#print 'Filename', fname
					#print links
					lst=links[fname]
					lst.append(f)
					#print lst
					links[fname]=lst
				else:
					links[fname]=[f]
		else:
			print "No such file %(f)s!" % vars()

	return links, pdbdata


def writepdbdata(pdbdata):
	#Saves PDB file parameters chosen by the user
	pdbfile=open('PDB_info.txt','w')
	pdbfile.write(",".join(['File','Protein segment','Ligand segment','Ligand name'])+"\n")
	for k in pdbdata.keys():
		general=pdbdata[k]
		for val, var in {'Protein segment': 'pseg', 'Ligand segment': 'lseg', 'Ligand name': 'lname'}.iteritems():
			#print val, var
			if val in general.keys():
				exec('{}="{}"'.format(var,str(general[val])))
			#assert var == general[val]
			else:
				exec(var+'=""')
			#assert var == ""
		datarow=[k,pseg,lseg,lname]
		pdbfile.write(",".join(datarow)+"\n")
	pdbfile.close()

def readpdbdata(pdbfile):
	#Reads saved user choices for ligands in analyzed PDB files
	pdbdata=NestedDict()
	pdbf=open(pdbfile,'r')
	rdr=csv.reader(pdbf, delimiter=',')
	data=[row for row in rdr]
	pdbf.close()
	header=data[0]
	hlabs=['Protein segment','Ligand segment','Ligand name','Complex']
	hmap={hl : header.index(hl) for hl in hlabs if hl in header }
	#psin=header.index('Protein segment')
	#lsin=header.index('Ligand segment')
	#lin=header.index('Ligand name')
	#cm=header.index('Complex')
	for d in data[1:]:
		print d
		for k in hmap.keys():
			pdbdata[d[0]][k]=d[hmap[k]]
	return pdbdata

def collect(links):
	#Collects data for McVol calculations according to the expected links. Used in case only PQR files were generated from PDB input.
	outputs=NestedDict()
	for k in links.keys():
		for f in links[k]:
			if f!='' and os.path.isfile(f):
				fpath, fname, ftype=filename(f)
				if '_' in fname:
					part=fname.split('_')[-1]
				else:
					part='vol'
			if os.path.exists('%(k)s/%(k)s_%(part)s.txt' % vars()):
				ofile=open('%(k)s/%(k)s_%(part)s.txt' % vars(), "r")
				outputs[k][part]=ofile.read()
				ofile.close()
			else:
				print "No results file %(k)s_%(part)s.txt" % vars()
	return outputs

def analyze(outputs):
	#Collects data from McVol results into map data structure which is later used for data saving.
	odata=NestedDict()
	for k in outputs.keys():
		lower=outputs[k]
		for l in lower.keys():
			data=lower[l]
			datalines=data.split('\n')
			cavvol=0
			clefvol=0
			for dl in datalines:
				words=dl.split()
				if 'Checking' in dl:
					vtype=words[1]
					vnum=int(words[2])
					vvol=float(words[4].replace('A^3',''))
					odata[k][l]['Volumes'][vnum]['Type']=vtype
					odata[k][l]['Volumes'][vnum]['Volume']=vvol
					odata[k][l]['Volumes'][vnum]['Number']=vnum
					if vtype=='cavity':
						cavvol=cavvol+vvol
					elif vtype=='cleft':
						clefvol=clefvol+vvol
				if 'printing' in dl:
					vnum=int(words[2])
					link=words[4]
					odata[k][l]['Volumes'][vnum]['Link']=link
				if 'waters placed' in dl:
					vnum=int(words[5])
					vwtr=int(words[0])
					odata[k][l]['Volumes'][vnum]['Waters']=vwtr
				if 'Protein Volume' in dl:
					vol=float(words[2])
					odata[k][l]['General']['Volume']=vol
				if 'Solvent Accessible Volume' in dl:
					vol=float(words[3])
					odata[k][l]['General']['SAV']=vol
				if 'VdW Volume' in dl:
					vol=float(words[2])
					odata[k][l]['General']['VdW']=vol
				if 'Total Surface' in dl:
					surf=float(words[2])
					odata[k][l]['General']['Surface']=surf
				if 'Water molecules placed' in dl:
					wtrs=int(words[0])
					odata[k][l]['General']['Total waters']=wtrs
			odata[k][l]['General']['Cavity volume']=cavvol
			odata[k][l]['General']['Cleft volume']=clefvol
	#print odata
	return odata

def pqrit(ifile):
	#Converts PDB to PQR by the use of obabel

	ipath, iname, itype=filename(ifile)
	#print ipath,iname,itype
	ofile='%(ipath)s/%(iname)s.pqr' % vars()
	print ofile
	command="obabel -i %(itype)s %(ifile)s -o pqr -O %(ofile)s" % vars()
	failure, output=runcmd(command)
	#print output
	if failure:
		ofile=''
	else:
		fixpqr(ofile)
	#Generated pqr file is not recognised by McVol and needs further modeifications done by fixpqr function
	return ofile

def splitpdb(ifile, odata, ligname=''):
	#PDB parsing engine that let's the user select PDB chain, ligand or ligands and saves their choices.
	info=NestedDict()
	ipath, iname, itype=filename(ifile)
	u=MD.Universe(ifile)
	print "PDB file {}".format(iname+'.'+itype)
	if iname in odata.keys():
		pdbinfo=odata[iname]
		pdbi=True
	else:
		pdbi=False

	#Collect PDB structure information
	for s in u.segments:
		print "Segment: {}".format(s.name)
		p=s.selectAtoms('protein')
		np=s.selectAtoms('not protein')
		pres=list(set([r.name for r in p.residues if r.name not in ['','HOH',None]]))
		npres=list(set([r.name for r in np.residues if r.name not in ['','HOH',None]]))
		npids=list(set([r.id for r in np.residues if r.name not in ['','HOH',None]]))
		if len(npids)>len(npres):
			print "\tMore than one residue with the same name!"
			#print "\t",npids
			segmented=open("Segmented.txt", "a")
			segmented.write(ifile+'\n')
			segmented.close()
			return [], odata
		if len(pres)>0 or len(npres)>0:
			info[s.name]['Protein']['Res']=pres
			info[s.name]['Protein']['Size']=int(len(pres))
			info[s.name]['Protein']['Atoms']=p
			info[s.name]['Hetero']['Res']=npres
			info[s.name]['Hetero']['Size']=int(len(npres))
			info[s.name]['Hetero']['Atoms']=np
			print "\tProtein: {}".format(pres)
			print "\tHeteroatoms: {}\n".format(npres)



	psegs=defaultdict(list)
	npsegs=defaultdict(list)

	for k in info.keys():
		if info[k]['Protein']['Size']!=0:
			psegs[info[k]['Protein']['Size']].append(k)
		if info[k]['Hetero']['Size']!=0:
			npsegs[info[k]['Hetero']['Size']].append(k)


	#Get protein segment
	if pdbi and 'Protein segment' in pdbinfo.keys() and pdbinfo['Protein segment']!='':
		pseg=pdbinfo['Protein segment']
		print "\nFound protein segment:", pseg
	else:

		if len(psegs.keys())==1 and len(psegs[psegs.keys()[0]])==1:
			pseg=psegs[psegs.keys()[0]][0]
			print "\nAuto protein segment:", pseg

		elif len(psegs.keys())==0:
			print "No protein found!"
			return []
		else:
			segs=[]
			for g in psegs.keys():
				segs.extend(psegs[g])
			for s in segs:
				print "Protein segment: {}".format(s)
				print "\tProtein: {}".format(info[s]['Protein']['Res'])
				if info[s]['Hetero']['Size']!=0:
					print "\tHeteroatoms: {}".format(info[s]['Hetero']['Res'])
			ask="\nProtein segment to use/q to quit/n to skip: "
			pseg=""
			while pseg not in segs:
				pseg = raw_input(ask)
				if pseg=="q":
					print "Exiting script!"
					sys.exit(1)
				if pseg=="n":
					print "Skipping!"
					return []

	prot=info[pseg]['Protein']['Atoms']




	#Choose ligand segment
	if pdbi and 'Ligand segment' in pdbinfo.keys():
		lseg=pdbinfo['Ligand segment']
		if lseg!='':
			print "\nFound ligand segment:", lseg
			l=info[lseg]['Hetero']['Atoms']
		else:
			print "\nNo ligand segment!"
			ligname=None
	else:
		if ligname!=None:
			if len(npsegs.keys())==1 and len(npsegs[npsegs.keys()[0]])==1:
				lseg=npsegs[npsegs.keys()[0]][0]
				print "Auto ligand segment:", lseg
				l=info[lseg]['Hetero']['Atoms']
			elif len(npsegs.keys())==0:
				print "No ligands found!"
				ligname=None
			else:
				legs=[]
				for g in npsegs.keys():
					legs.extend(npsegs[g])
				for ls in legs:
					print "Ligand segment: {}".format(ls)
					liga=info[ls]['Hetero']['Atoms']
					for lres in info[ls]['Hetero']['Res']:
						liga_n=len(liga.selectAtoms('resname {}'.format(lres)))
						print "\t{} with {} atoms.".format(lres,liga_n)
				lask="Ligand segment to use/q to quit/c to continue/n to skip: "
				print "Protein segment that was chosen: {}".format(pseg)
				lseg=""
				while lseg not in legs:
					lseg = raw_input(lask)
					if lseg=="q":
						print "Exiting script!"
						sys.exit(1)
					elif lseg=="c":
						print "No ligand will be used!"
						ligname=None
						break
					elif lseg=="n":
						print "Skipping!"
						return []
				l=info[lseg]['Hetero']['Atoms']


	#Choose ligand


	if pdbi and 'Ligand name' in pdbinfo.keys():
		if pdbinfo['Ligand name']!='':
			ligname=pdbinfo['Ligand name']
			if ';' in ligname:
				ligname=ligname.split(";")
			print "\nFound ligand:", ligname
		else:
			ligname=None
			print "\nNo ligand!"
	else:
		if ligname!=None:
			pligs=info[lseg]['Hetero']['Res']
			if ligname=='':
				pligs=[pli for pli in pligs if pli not in ['MG','SO4']]
				for lres in pligs:
					liga_n=len(l.selectAtoms('resname {}'.format(lres)))
					print "\tLigand: {} with {} atoms.".format(lres,liga_n)
				if len(pligs)==1: #Recent changes here!!!
					ligname=pligs[0]
				else:
					qstn = "Enter ligand residue name or two ligands separated by 'and'/q to quit/c to continue/n to skip: "
					ligands=['']
					while len([lig for lig in ligands if not lig in pligs])!=0:
						ligname = raw_input(qstn)
						if ligname=="q":
							print "Exiting script!"
							sys.exit(1)
						elif ligname=="c":
							print "No ligand will be used!"
							ligname=None
							break
						elif ligname=="n":
							print "Skipping!"
							return []
						elif 'and' in ligname:
							ligands=[li.strip() for li in ligname.split('and')]
							ligands=list(set(ligands))
						#print ligands
						else:
							ligands=[ligname]

					if ligname!=None:
						ligands=list(set(ligands))
						if len(ligands)==1:
							ligname=ligands[0]
						else:
							ligname=ligands

			else:
				if not ligname in pligs:
					ligname=None
					print "No such ligands {} found!".format(ligname)

	if pdbi and 'Complex' in pdbinfo.keys():
		ionname=pdbinfo['Complex']
		ion=l.selectAtoms('resname {}'.format(ionname))
		prot=prot+ion


	lig_o=''
	plig_o=''
	prot_o=iname+"/"+iname+"_P.pdb"





	if pseg in ['SYSTEM','']:
		prot.set_segid('P')



	if ligname!=None:
		odata[iname]['Ligand segment']=lseg
		if isinstance(ligname, list):
			lig=l.selectAtoms('protein')
			i=1
			odata[iname]['General']['Ligand name']=";".join(ligname)
			for ligs in ligname:
				lt=l.selectAtoms('resname {}'.format(ligs))
				lt.set_resid(i)
				lig=lig+lt
				i=i+1
		else:
			lig=l.selectAtoms('resname {}'.format(ligname))
			lig.set_resid(1)
			odata[iname]['Ligand name']=ligname
		if lseg in ['SYSTEM','']:
			lig.set_segid('L')

		plig=prot+lig
		lig_o=iname+"/"+iname+"_L.pdb"
		plig_o=iname+"/"+iname+"_PL.pdb"
		plig.write(plig_o)
		lig.write(lig_o)




	prot.write(prot_o)
	out=[prot_o, lig_o, plig_o]
	out=[i for i in out if i!='']

	odata[iname]['Protein segment']=pseg


	return out, odata

def fixpqr(ifile):
	#Fixes pqr format to the one understandable by McVol

	#ifl=open(ifile, 'r')
	#idata=ifl.read()
	#ifl.close()
	table=tableout(ifile)
	ifl=open(ifile, 'w')
	for data in table:
		if len([el for el in data if el in ['CONECT','MASTER','COMPND','AUTHOR','HOH']])==0:
			#not 'CONECT' in data and not 'MASTER' in data and not 'COMPND' in data and not 'AUTHOR' in data and not 'HOH' in data:

			#print data
			if len(data)==0:
				continue
			elif 'END' in data:
				#print 'END'
				ifl.write(data[0]+'\n')
			else:
				#print data
				#data[1]=int(data[1])
				#data[4]=int(data[4])
				#print data
				if len(data[2])==4:
					datastr='{0: <5} {1:5d} {2: <4} {3: <4} {5: 4d} {6: 11.3f} {7: 8.3f} {8: 8.3f} {9:.3f} {10:.3f}'.format(*data)
				else:
					datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {5: 4d} {6: 11.3f} {7: 8.3f} {8: 8.3f} {9:.3f} {10:.3f}'.format(*data)
				#print datastr
				ifl.write(datastr+'\n')
	ifl.close()


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

def numerize(s):
	#Convert values to numbers
	try:
		if s=='NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s)==0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s

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

