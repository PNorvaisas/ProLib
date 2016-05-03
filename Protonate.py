import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import chimera


def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)


odir="Aligned-Protonated"
ensure_dir(odir)
# change to folder with data files
#os.chdir("/Users/pett/data")

# gather the names of .pdb files in the folder
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

#print file_names
#print os.getcwd()
# loop through the files, opening, processing, and closing each in turn
fprot=[]
selection=file_names#[:10]
for fn in selection:
	print "------------------------------------------------------------"
	print "{} - {}/{}".format(fn,selection.index(fn)+1,len(selection))
	replyobj.status("Processing " + fn) # show what file we're working on
	rc("open " + fn)
	try:
		rc("addh")
	except:
		print "Protonation unsuccessful!"
		fprot.append(fn)
	nomodels=len(chimera.openModels.list())
	if nomodels>1:
		rc("mm #0 #"+str(nomodels-1)+" computeSS false")
	rc("write "+str(nomodels-1)+" "+odir+"/"+fn)

if len(fprot)>0:
	fprot=open('Failed_protonation.txt','w')
	for item in fails:
		fails.write("%s\n" % item)
	fprot.close()

rc("close all")
