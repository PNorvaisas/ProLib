#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')

# for mod in ['pip','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','time']:
# 	try:
# 		exec "import %(mod)s" % vars()
# 	except ImportError, e:
# 		print "Module not found %(mod)s\nTrying to install!" % vars()
# 		install(mod)
# 	#pass # module doesn't exist, deal with it.


from Sequencing import *
import numpy as np

def clean(val):
	trash=['<','=','>']
	for t in trash:
		val=val.replace(t,'')
	if val!='':
		val=float(val)
	return val


afffile=sys.argv[1]

afpath, afname, aftype=filename(afffile)

aff=open(afffile,'r')
rdr=csv.reader(aff, delimiter=',')
data=[row for row in rdr]
aff.close()


fixed=[]
header=data[0]
hlabs=['PDB ID','Chain ID','HET ID','Ki (nM)','Kd (nM)']
hmap={hl : header.index(hl) for hl in hlabs if hl in header }
data=[row for row in data if len(row)!=0 and (row[hmap['Ki (nM)']]!='' or row[hmap['Kd (nM)']]!='')]
vallen=0
for row in data[1:]:
	print row[hmap['PDB ID']],row[hmap['HET ID']],'Ki:', row[hmap['Ki (nM)']],'Kd:', row[hmap['Kd (nM)']]
	kis=[ val.strip().split(' ')[0].strip() for val in row[hmap['Ki (nM)']].split(',')] #.split(' ')[0].strip()
	kds=[ val.strip().split(' ')[0].strip() for val in row[hmap['Kd (nM)']].split(',')] #.split(' ')[0].strip()
	kis=[ val.split('-') if '-' in val else [val] for val in kis]
	kds=[ val.split('-') if '-' in val else [val] for val in kds]
	kisc=[]
	kdsc=[]
	for val in kis:
		for valu in val:
			valu=clean(valu)
			if valu!='':
				kisc.append(valu)
	for val in kds:
		for valu in val:
			valu=clean(valu)
			if valu!='':
				kdsc.append(valu)
	logki=np.power(10,(np.log10(kisc)-9))
	logkid=np.power(10,(1/np.log10(kdsc))-9 )
	print logki, logkid
	values=list(logki)+list(logkid)
	if len(values)>vallen:
		vallen=len(values)
	fixed.append([row[hmap['PDB ID']],row[hmap['HET ID']]]+values)
	#print np.power(10,(np.log10(kisc)-9)), np.power(10,(1/(np.log10(kdsc)-9)))

fixed=[row+['']*(vallen+2-len(row)) for row in fixed]
fixed=[['PDB ID','HET ID']+['Ki_M']*vallen]+fixed

writecsv(fixed,afname+'_fixed.csv',delim=',')