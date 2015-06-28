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

import textwrap as tw
import itertools as IT

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#from pylab import *
import numpy as np

class NestedDict(dict):
	def __getitem__(self, key):         
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

f=open('All_alignments.csv','r')
cr=csv.reader(f)
idata=[ln for ln in cr]
f.close()

r=open('1uyl/1uyl_aligned.csv','r')
rr=csv.reader(r)
ref=[ln for ln in rr]
r.close()

data=NestedDict()
headers=idata[0]

for ln in idata[1:]:
	if isinstance(data[ln[0]][ln[1]],list):
		data[ln[0]][ln[1]]=data[ln[0]][ln[1]]+[ln[2:]]
	else:
		data[ln[0]][ln[1]]=[]

transposed=[]
ref=list(zip(*ref)[3][1:])
transposed.append(ref)
for fl in data.keys():
	results=data[fl]['R_P']
	refs=zip(*results)[1]
	vols=zip(*results)[9]
	resline=[]
	for r in ref:
		if r in refs:
			resline.append(vols[refs.index(r)])
		else:
			resline.append(None)
	transposed.append(resline)
print ref
print headers


#tr=open('R_P_results.csv','w')
#ctr=csv.writer(tr,dialect='excel')
#for ln in transposed:
#	ctr.writerow(ln)
#tr.close()

volumes=zip(*transposed)
#print volumes[1:]
#volumes=[[float(cell) for cell in ln if cell!=''] for ln in volumes[1:]]
#volumes=transposed.asarray()
transposed=np.asarray(transposed)

labels=transposed[0]
print transposed
print transposed.shape

plt.boxplot(transposed[1:],labels=labels)
plt.show()



