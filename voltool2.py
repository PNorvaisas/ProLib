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
	-r <file>    Name of reference file to use (only if input and reference files are in PDB!)
	-s <file>    Name of McVol settings file (if not given default settings will be used)
	-p <file>    PDB info file to use (contains information about the structure of PDB files and ligand choices)
				 Will be generated in the directory for each run, default name "PDB_info.txt"

Options:
	pqr	     Only convert selected files to pqr 
	mcvol	     Calculate volumes
	analyze	     Analyze McVol output and generate tables
	load	     Load previously saved data and continue
	 

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
Reference:	%(rfile)s
   PQR it:	%(pqr)s 
  Volumes:	%(mcvol)s
  Analyze:	%(doanalyze)s
     Load:	%(load)s
<--------------------------------------------->
	'''



def main(argv=None):
	dset={'nmc':200, 'surfPT':2500, 'probe':1.3, 'membZmin':-100, 'membZmax':100, 'startgridspacing':1.0, 'membgridspacing':2.0, 'cavgridspacing':0.1, 'minVol':1, 'waterVol':18, 'DummyRad':1.7, 'blab':0, 'MembDim':4, 'CoreZMin':0, 'CoreZMax':0, 'CoreDim':0, 'CleftDim':3, 'CleftRel':70, 'CleftMethod':2, 'SurfaceCluster':1.5}
	settmp=''
	for ds in dset.keys():
		print ds, dset[ds]
		settmp=settmp+'{} {}\n'.format(ds,dset[ds])
	#dsettings='/home/povilas/Scripts/McVol_settings.txt'
	#setftmp=open(dsettings,'r')
	#settmp=setftmp.read()
	#setftmp.close()
	ifile=""
	sfile=''
	rfile=''
	rname=''
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
			opts, args = getopt.getopt(argv[1:], "hi:r:s:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return	
			if option in ("-i", "--input"):
				ifile=value	
			if option in ("-s", "--settings"):
				sfile=value
			if option in ("-p", "--pdbinfofile"):
				pdbinfofile=value
			if option in ("-r", "--reffile"):
				rfile=value
	
	
		for argument in args:		
			if argument in ("pqr", "--onlypqr"):
				pqr = True
				mcvol=False
				doanalyze=False
			if argument in ("mcvol", "--mcvol"):
				mcvol = True
				pqr=False
				doanalyze=False
			if argument in ("analyze", "--analyze"):
				doanalyze = True
				mcvol=False
				pqr=False
				load=False
			if argument in ("load", "--load"):
				load=True
				doanalyze = False
				mcvol=False
				pqr=False
			


	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	#Check for the input integrity



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
		
	
		if rfile!='':
			if rfile=='all':
				rname='all'
			else:
				rpath, rname, rtype = filename(rfile)
				if os.path.isfile(rfile):			
					if not rtype in ("pqr", "pdb"):
						raise Exception("File %(rfile)s is not supported" % vars())
						rname=''
					else:
						print "File %(rfile)s is expected to be a structure file" % vars()
						#reference=True
		

		
			
	
	except Exception, e:
		print e
		#print "!!!-------------------------------------------!!!"
		#usage()
		sys.exit(1)


	#----------------------------


	print optionsset %vars()
	if not load :
		if sfile!='' and os.path.isfile(ifile):
			spath, sname, stype = filename(sfile)
			setf=open(sfile,'r')
			settmp=setf.read()
			setf.close()


		#ilist=[ifile]

		ilist=genlist(ifile)
		print "Input: {}\n".format(', '.join(ilist))
		if rfile!='' and rfile!='all':
			rlist=genlist(rfile)
			print 'Reference:{}\n'.format(rfile)
			if rfile not in ilist:
				ilist.extend(rlist)
		else:
			rlist=[]	
	
	
	

	
	if os.path.isfile(pdbinfofile):
		ask='Do you want to use PDB info file {}? (Y/N) '.format(pdbinfofile)
		answ=raw_input(ask)
		if answ=='Y':
			pdbdata=readpdbdata(pdbinfofile)
		else:
			print 'PDB info file {} not used!'.format(pdbinfofile)
	else:
		answ='N'


	if pqr:
		links, pdbdata=prepare(ilist,pdbdata,pqr)
		if answ=='N' and not os.path.isfile(pdbinfofile):
			writepdbdata(pdbdata)		

	if mcvol:
		if not pqr:
			links, pdbdata=prepare(ilist,pdbdata,pqr)
		outputs=volumeit_M(links,settmp,16)
		#sys.exit(1)
		saveout(outputs)

	if doanalyze:
		if not pqr and not mcvol:			
			links, pdbdata=prepare(ilist,pdbdata,pqr)
			outputs=collect(links)

		odata=analyze(outputs)

		odata=joincav(odata)

		odata=vmdit(odata,rname)

	

		odata = findcenters(odata)

		f = open('Odata_center.pckl', 'w')
		pickle.dump(odata, f)
		f.close()
		
		#results=align_M(odata,rname,20)

		#f = open('Odata_results.pckl', 'w')
		#pickle.dump([odata,results], f)
		#f.close()

		print 'It works'
		#sys.exit(1)
		#odata=finddiff(odata,rname)

		#tables=summarize(odata,pdbdata)
		#writesheets(tables)
	if load or doanalyze:
		if os.path.isfile('Odata_center.pckl'):
			f = open('Odata_center.pckl', 'rb')
			odata=pickle.load(f)
			f.close()
			tables=summarize(odata,pdbdata)
		#results=align_M(odata,rname)
			if os.path.isfile('Odata_results.pckl'):			
				f = open('Odata_results.pckl', 'rb')
				odata,results=pickle.load(f)
				f.close()
				restab=summalign(odata,results)
				tables['Alignment']=restab
		
		
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


def summalign(odata,results):
	aligntable=[]
	header=['Reference','R Frag','Target','T Frag','Alignment','R index','R Type','R Volume','R Waters placed','Rx','Ry','Rz','T index', 'T Type','T Volume','T Waters placed','Tx','Ty','Tz']
	aligntable.append(header)	
	vkeys=['Type','Volume','Waters','Center']	
	for k in results.keys():
		ref,tgt=k.split('_')
		for f in results[k].keys():
			ref_f,tgt_f=f.split('_')
			rdat=odata[ref][ref_f]['Volumes']
			tdat=odata[tgt][tgt_f]['Volumes']
			for al in results[k][f]:
				#print k,f,al,al[0],al[1],al[2]
				alin=al[0]
				rin=al[1]
				tin=al[2]
				rline=[ref,ref_f,tgt,tgt_f,alin]
				if rin!='':
					rline.append(rin)
					rline=getcavdat2(rdat,rin,rline,vkeys)
				else:
					rline.extend(['','','','','','',''])
				if tin!='':
					rline.append(tin)
					rline=getcavdat2(tdat,tin,rline,vkeys)
				else:
					rline.extend(['','','','','','',''])
				
				aligntable.append(rline)
				
	return aligntable


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

def getcavdat2(dataset,c,line, dkeys):
	#Generates links to McVol results					
	if c in dataset.keys():
		vol=dataset[c]
		for d in dkeys:
			if d=='Center':
				cnt=vol[d]
				line.extend([cnt[0],cnt[1],cnt[2]])
			else: 
				line.append(vol[d])
	else:
		line.extend(['','','','','',''])

	return line



def joincav(odata):
	#Joins cavities or clefts from multiple files to one file
	for k in odata.keys():
		frags=odata[k]			
		for f in frags.keys():
			cavdata=odata[k][f]['Volumes']
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


def align_ref(odata,rname):
	#Finds whether cavities or clefts found with McVol are overlapping. Checks whther the center point of one is within a determined distance of the points of the other.
	#odata = findcenters(odata)
	refs=['P']
	tgts=['PL']			

	if rname!='' and rname in odata.keys():
		reference=True
		print "Setting up reference {}.\n".format(rname)
		rfrag=odata[rname]
		refvols=rfrag['P']['Volumes']
		refs=['R','P','R']
		tgts=['P','PL','PL']
		
	else:
		print "No reference used!"
		

	for k in odata.keys():
		#if k==rname:
			#continue
		print '<<<Aligning {}>>>'.format(k)
		frag=odata[k]
		
		num=1
		for ref,tgt in IT.izip(refs,tgts):
			if tgt in frag.keys():
				print "{} to {}".format(ref,tgt)
				rt={}
				if ref=='R':
					rvols=odata[rname]['P']['Volumes']
				else:
					#This will cause problems
					rvols=frag[ref]['Volumes']
				tvols=frag[tgt]['Volumes']
				
			
				for rvol in sorted(rvols.keys()):				
					coords=rvols[rvol]['Center']				
					for tvol in sorted(tvols.keys()):
						tlinkt=tvols[tvol]['Link']
						tlink='{}/{}/'.format(k,tgt)+tlinkt					
						if inpqr(tlink,coords):
							print '\t{}{}={}{}'.format(ref,rvol,tgt,tvol)
							rt[rvol]=tvol				
			
					
				rmiss=[c for c in rvols.keys() if not c in rt.keys()]
				rin=[c for c in rvols.keys() if c in rt.keys()]
				tmiss=[c for c in tvols.keys() if not c in rt.values()]
				tin=[c for c in tvols.keys() if c in rt.values()]
				print 'Overlap: {}:{}/{}:{}, {}diff: {}, {}diff: {}\n------------------------------------------\n'.format(ref,len(rin),tgt,len(tin),ref,len(rmiss),tgt,len(tmiss))	
				#if not reference:
				num=1

				for ri in sorted(rin):
					odata[k][ref][ref+'_'+tgt][num]=rvols[ri]
					#Backwards compatibility
					odata[k][tgt][tgt+'_'+ref][num]=tvols[rt[ri]]
					num=num+1
				for rm in sorted(rmiss):
					odata[k][ref][ref+'_'+tgt][num]=rvols[rm]
					num=num+1
				for tm in sorted(tmiss):
					odata[k][tgt][tgt+'_'+ref][num]=tvols[tm]
					num=num+1
										

	return odata



def finddiff(odata,rname):
	refs=['P']
	tgts=['PL']
	if rname!='' and rname in odata.keys():
		reference=True
		print "Setting up reference {}\n".format(rname)
		refs=['R','P']
		tgts=['PL','PL']
		
	else:
		print "No reference used!"

	for k in odata.keys():
		print '<<<Finding volume differences for {}>>>'.format(k)
		for ref,tgt in IT.izip(refs,tgts):
			if ref in odata[k].keys() and ref+'_'+tgt in odata[k][ref].keys():
				if ref=='R':
					rpath=rname+'/P/'
					tpath='{}/{}/'.format(k,tgt)
				elif tgt=='R':
					rpath='{}/{}/'.format(k,ref)
					tpath=rname+'/P/'
				else:
					rpath='{}/{}/'.format(k,ref)
					tpath='{}/{}/'.format(k,tgt)
				diffdirs=['{}/{}-{}'.format(k,ref,tgt),'{}/{}-{}'.format(k,tgt,ref)]
				for d in diffdirs:
					if not os.path.exists(d):
						os.makedirs(d)
				if ref in ['R','P']:
					common=True
				else:
					common=False
				reftgt=odata[k][ref][ref+'_'+tgt]
				tgtref=odata[k][tgt][tgt+'_'+ref]	
				print 'Calculating {} & {}'.format(ref,tgt)
				for a in range(1,max(reftgt.keys()+tgtref.keys())):
					if a in tgtref.keys() and a in reftgt.keys():
						pqr2=tpath+str(tgtref[a]['Link'])
						pqr2_no=tgtref[a]['Number']
						pqr1=rpath+str(reftgt[a]['Link'])
						pqr1_no=reftgt[a]['Number']
						pqr1v=volumize(pqr1)
						pqr2v=volumize(pqr2)
						d12, d21, a12=voldiff(pqr1v,pqr2v)
						#print d12, d21, a12
						d12_c=volcalc(d12)
						d21_c=volcalc(d21)
						a12_c=volcalc(a12)	
			
						#d12_c=volsave(d12,'{}/{}-{}/{}-{}.pqr'.format(k,ref,tgt,pqr1_no,pqr2_no),pqr1_no,0.2)
						#d21_c=volsave(d21,'{}/{}-{}/{}-{}.pqr'.format(k,tgt,ref,pqr2_no,pqr1_no),pqr2_no,0.2)
						#a12_c=volsave(a12,'{}/{}-{}/{}&{}.pqr'.format(k,ref,tgt,pqr1_no,pqr2_no),pqr1_no,0.2)							
						#dcount, ccount=voldiff(pqr1,pqr2,diffdir,common)
						print '{}-{}={}A {}-{}={}A {}&{}={}A'.format(pqr1_no,pqr2_no,d12_c*0.001,pqr2_no,pqr1_no,d21_c*0.001, pqr1_no,pqr2_no,a12_c*0.001)
					
	


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

def volumize(pqr):
	volume=NestedDict()
	data=tableout(pqr)
	ttable=zip(*data)
	#print ttable
	xc=ttable[5]
	yc=ttable[6]
	zc=ttable[7]
	for x,y,z in IT.izip(xc,yc,zc):
		if len(volume[x][y])==0:
			volume[x][y]=[z]
		else:
			volume[x][y]=volume[x][y]+[z]
			
	return volume

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

def volcalc(pqr):
	pqrf=open(loc,'w')
	count=0
	for x in pqr.keys():
		ys=pqr[x]
		for y in ys.keys():
			for z in ys[y]:
				count=count+1
	return count

def voldiff(pqr1,pqr2):
	d12=copy.deepcopy(pqr1)
	d21=copy.deepcopy(pqr2)
	a12=NestedDict()
	for x1 in pqr1.keys():
		y1s=pqr1[x1]
		for y1 in y1s.keys():
			z1s=y1s[y1]
			for z1 in z1s:
				#print z1
				#print d21[x1][y1]
				if z1 in pqr2[x1][y1]:
					d12[x1][y1]=[z for z in d12[x1][y1] if z!=z1]
					d21[x1][y1]=[z for z in d21[x1][y1] if z!=z1]
					if len(d12[x1][y1])==0:
						d12[x1].pop(y1)
					if len(d12[x1].keys())==0:
						d12.pop(x1)
					if len(d21[x1][y1])==0:
						d21[x1].pop(y1)
					if len(d21[x1].keys())==0:
						d21.pop(x1)
					if len(a12[x1][y1])==0:
						a12[x1][y1]=[z1]
					else:
						a12[x1][y1]=a12[x1][y1]+[z1]

	

	return d12, d21, a12


			


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

def vmdit(odata,rname):
	#Generates VMD script chich allows quick visualization of McVol output
	script={'Header':'#!/usr/bin/tclsh\n',
		'Protein':'mol new %(Protein)s\nset P %(molnr)d\n',
		'Pcav':'mol new %(Pcav)s\nset Pcav %(molnr)d\nmol modstyle 0 $Pcav Lines\n',
		'Pclef':'mol new %(Pclef)s\nset Pclef %(molnr)d\nmol modstyle 0 $Pclef Lines\nmol modcolor 0 $Pclef Element\n',
		'Ligand':'mol new %(Ligand)s\nset L %(molnr)d\nmol modstyle 0 $L Licorice\n',
		'PLcav':'mol new %(PLcav)s\nset PLcav %(molnr)d\nmol modstyle 0 $PLcav Points 16.0\n',
		'PLclef':'mol new %(PLclef)s\nset PLclef %(molnr)d\nmol modstyle 0 $PLclef Points 16.0\nmol modcolor 0 $PLclef Element\n'
		}
	
	files={}
	
	for k in odata.keys():		
		frags=odata[k]
		for f in frags.keys():
			general=odata[k][f]['General']
			if f in ['P','vol']:
				files['Protein']='{}_{}.pdb'.format(k,f)		
			if f=='L':
				files['Ligand']='{}_{}.pdb'.format(k,f)
			if f in ['P','PL']:
				if 'Cavity link' in general.keys():
					files[f+'cav']=general['Cavity link'].split('/')[1]
				if 'Cleft link' in general.keys():
					files[f+'clef']=general['Cleft link'].split('/')[1]
		#print files
		#print general
		scrfile=open('{}/{}.tcl'.format(k,k),'w')
		scrfile.write(script['Header'])
		files['molnr']=0
		for e in ['Protein','Pcav','Pclef','Ligand','PLcav','PLclef']:
			if e in files.keys():
				string=script[e]
				scrfile.write(string % files)
 				files['molnr']=files['molnr']+1
		scrfile.close()
		odata[k][f]['General']['Tcl script']='{}.tcl'.format(k)
	return odata




					
	

def volumeit(links,settings):
	#McVol execution
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
	#McVol execution
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
	#Writes organized data to file.
	#odir=dircheck('Split')
	for i in sheets.keys():
		if i in ['Summary','Volumes','Alignment']:
			oname=i+".csv"
		else:
			print 'Unknown key for table: {}'.format(i)
			

		#else:
		#	if '_' in i:
		#		tp=i.split('_')[1]
		#		nm=i.split('_')[0]
		#	else:
		#		nm=i
		#		tp='results'	
		#	oname="{}/{}_{}.csv".format(nm,nm,tp)
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
						pqrfile=pqrit(p)
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
	psin=header.index('Protein segment')
	lsin=header.index('Ligand segment')
	lin=header.index('Ligand name')
	cm=header.index('Complex')
	for d in data[1:]:
		print d
		pdbdata[d[0]]['Protein segment']=d[psin]
		pdbdata[d[0]]['Ligand segment']=d[lsin]
		pdbdata[d[0]]['Ligand name']=d[lin]
		pdbdata[d[0]]['Complex']=d[cm]
	
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
	#Collects data from McVol results data into map data structure which is later used for data saving.
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

	return ofile

def splitpdb(ifile, odata, ligname=''):
	#PDB parsing engine that let's the user select ligand or ligands and saves their choices.
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
	failure, output = commands.getstatusoutput(cmd)
	if failure:
		print '''Running failed \n %(cmd)s \n %(output)s'''.encode('utf-8') % vars();
	return failure, output

class NestedDict(dict):
	def __getitem__(self, key):         
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

def numerize(s):
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



def install(package):
	pip.main(['install', package])




def align_M(odata,rname,cores):
	print 'Aligning cavities and clefts...'


	outputs=NestedDict()

	p=Pool(cores)
	m = Manager()
	q = m.Queue()
	args=[]

	if rname=='all':
		for ref in odata.keys():
			for tgt in odata.keys():
				args.append((odata,ref,tgt,q))
	elif rname!='' and rname in odata.keys():
		for tgt in odata.keys():
			args.append((odata,rname,tgt,q))

	result = p.map_async(align, args)
	start=time.time()
	prcprev=0

	while True:
		if result.ready():
			break
		else:
			prc = float(q.qsize()+1)*100/float(len(args))
			timepassed=float(time.time()-start)
			if prc>prcprev:
				if prc!=0:
					total=int((timepassed*100/prc))
					remaining=int((timepassed*100/prc)-timepassed)
				else:
					total=int((timepassed*100/prc))
					remaining=int('inf')
				print "Done {0:3d}%, remaining: {1:<5} total: {2:<5}".format(int(prc),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
				prcprev=prc
		time.sleep(2)

	print result
	results=result.get()

	for res in results:
		rn=res[0]
		tn=res[1]
		#print rn,tn,res[2]
		outputs.update(res[2])

	return outputs



def align((odata,rname,k,q)):
	#Finds whether cavities or clefts found with McVol are overlapping. Checks whther the center point of one is within a determined distance of the points of the other.

	#refs=['P']
	#tgts=['PL']

	results=NestedDict()

	if rname!='' and rname in odata.keys():
		reference=True
		#print "Setting up reference {}.\n".format(rname)
		rfrag=odata[rname]
		#refvols=rfrag['P']['Volumes']
		refs=['P','P','PL','PL']

	else:
		print "No reference used!"
	if k!='' and k in odata.keys():
		tgts=['P','PL','P','PL']
		frag=odata[k]




	#if k==rname:
		#continue
	#print '<<<Aligning {} to {}>>>'.format(k,rname)


	num=1
	for ref,tgt in IT.izip(refs,tgts):
		if (tgt in frag.keys()) and (ref in rfrag.keys()):
			#print "{}-{} to {}-{}".format(rname,ref,k,tgt)
			rt={}
			rvalues=[]
			rvols=rfrag[ref]['Volumes']

			tvols=frag[tgt]['Volumes']


			for rvol in sorted(rvols.keys()):
				coords=rvols[rvol]['Center']
				for tvol in sorted(tvols.keys()):
					tlinkt=tvols[tvol]['Link']
					tlink='{}/{}/'.format(k,tgt)+tlinkt
					if inpqr(tlink,coords):
						#print '\t{}{}={}{}'.format(ref,rvol,tgt,tvol)
						if rvol in rt.keys():
							rt[rvol]=rt[rvol]+[tvol]
						else:
							rt[rvol]=[tvol]

						#What if couple target volumes?

			for rv in rt.values():
				if isinstance(rv,list):
					rvalues.extend(rv)
				else:
					rvalues.append(rv)



			rmiss=[c for c in rvols.keys() if not c in rt.keys()]
			rin=[c for c in rvols.keys() if c in rt.keys()]
			#Problems here
			tmiss=[c for c in tvols.keys() if not c in rvalues]
			tin=[c for c in tvols.keys() if c in rvalues]
			#print 'Overlap: {}:{}/{}:{}, {}diff: {}, {}diff: {}\n------------------------------------------\n'.format(ref,len(rin),tgt,len(tin),ref,len(rmiss),tgt,len(tmiss))
			#if not reference:
			num=1
			rtres=[]
			for ri in sorted(rin): # His
				if isinstance(rt[ri],list):
					for tgv in rt[ri]:
						rtres.append([num,ri,tgv])
						num=num+1
				else:
					rtres.append([num,ri,rt[ri]])
					num=num+1
			for rm in sorted(rmiss): #Reference missing
				rtres.append([num,rm,''])
				num=num+1
			for tm in sorted(tmiss): #Target missing
				rtres.append([num,'',tm])
				num=num+1

			results[rname+'_'+k][ref+'_'+tgt]=rtres


	q.put('i')

	return rname, k, results

#----------------------------------

if __name__ == "__main__":
	sys.exit(main())

