#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
"""
Volsave.py

Created by Povilas Norvaisas on 2014-04-11.
Copyright (c) 2014. All rights reserved.

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

ifile=sys.argv[1]
tmp=sys.argv[2]


def csvreader(ifile):
	ifl=open(ifile,'r')
	rdr=csv.reader(ifl, delimiter=',')
	data=[ln for ln in rdr]
	ifl.close()
	return data



def volsave(table,rad,mode):
	
	countca=0
	countcl=0
	headers=table[0]
	xi=headers.index('x')
	yi=headers.index('y')
	zi=headers.index('z')
	indi=headers.index('hcluster')
	if mode=='Cavities':
		cav=open('Cavities.pqr','w')
		clf=open('Clefts.pqr','w')
		tpi=headers.index('Type')
	elif mode=='Clusters':
		clus=open('Clusters.pqr','w')
		
	for r in table[1:]:
		
		x=float(r[xi])
		y=float(r[yi])
		z=float(r[zi])
		ind=int(r[indi])


		if mode=='Cavities':
			tp=r[tpi]
			if tp=='cavity':
				countca=countca+1
				data=['ATOM',countca,'DUM','DUM',ind,x,y,z,ind,rad]
				print data
				datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}\n'.format(*data)
				print datastr
				cav.write(datastr)
			if tp=='cleft':
				countcl=countcl+1
				data=['ATOM',countcl,'DUM','DUM',ind,x,y,z,ind,rad]
				
				datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}\n'.format(*data)
				print datastr
				clf.write(datastr)

		elif mode=='Clusters':
			countcl=countcl+1
			data=['ATOM',countcl,'DUM','DUM',ind,x,y,z,ind,rad]
			
			datastr='{0: <5} {1:5d}  {2: <3} {3: <4} {4: 4d} {5: 11.3f} {6: 7.3f} {7: 7.3f} {8:5.2f} {9:5.2f}\n'.format(*data)
			print datastr
			clus.write(datastr)

	if mode=='Cavities':
		cav.close()
		clf.close()

	elif mode=='Clusters':
		clus.close()


if '.csv' in tmp:
	cavs=csvreader(ifile)		
	volsave(cavs,1.7,'Cavities')
	clusts=csvreader(tmp)
	volsave(clusts,1.7,'Clusters')
elif tmp in ['Cavities','Clusters']:
	data=csvreader(ifile)
	volsave(data,1.7,tmp)

