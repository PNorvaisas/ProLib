import os
import sys
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import chimera


def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)


odir="Aligned-Protonated"
ensure_dir(odir)
# change to folder with data files
#sos.chdir("/Users/pett/data")

# gather the names of .pdb files in the folder
#file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

#print file_names
#print os.getcwd()
# loop through the files, opening, processing, and closing each in turn
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

#for modin in range(len(chimera.openModels.list())):
#	rc("select " + str(modin))
#	if modin>0:
#		rc("mm #0 #"+str(modin)+" computeSS false")
#	rc("addh")
#	rc("write "+str(modin)+" "+str(modin))
	
#print len(chimera.openModels.list())

rc("close all")
# uncommenting the line below will cause Chimera to exit when the script is done
#rc("stop now")
# note that indentation is significant in Python; the fact that
# the above command is exdented means that it is executed after
# the loop completes, whereas the indented commands that 
# preceded it are executed as part of the loop.