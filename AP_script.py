import os
import sys
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import chimera


def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)


odir='$odir'
ensure_dir(odir)

fails=[]
nomodels=len(chimera.openModels.list())

for mol in chimera.openModels.list():
	print mol.name
	print mol.id

if nomodels>1:
	models=chimera.openModels.list()[1:]
else:
	models=chimera.openModels.list()

for mol in models:
	print "Working on {}".format(mol.name)
	rc("select " + str(mol.id))
	try:
		rc("addh")
	except:
		print "Protonation unsuccessful!"
		fails.append(mol.name)
	if nomodels>1:
		rc("mm #0 #"+str(mol.id)+" computeSS false")
	rc("write "+str(mol.id)+" "+odir+"/"+str(mol.name))

if len(fails)>0:
	if os.path.exists('Failed_protonation.txt'):
		fprot=open('Failed_protonation.txt','a')
	else:
		fprot=open('Failed_protonation.txt','w')
	for item in fails:
		fails.write("%s\n" % item)
	fprot.close()

rc("close all")

rc("stop now")

