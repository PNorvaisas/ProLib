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

import Bio.SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import textwrap as tw
import itertools as IT
from os.path import basename

try:
	import MDAnalysis as MD
except ImportError, e:
	print "Module MDAnalysis not found"



from Networking import writecsv
from Sequencing import *


aas="Ala A, Arg R, Asn N, Asp D, Cys C, Glu E, Gln Q, Gly G, His H, Ile I, Leu L, Lys K, Met M, Phe F, Pro P, Ser S, Thr T, Trp W, Tyr Y, Val V"
aa_list=aas.split(',')
aa_list=[aa.split() for aa in aa_list]
aa={ a.strip().upper() : b.strip() for a,b in aa_list }




#ifile=sys.argv[1]
data=[]
header=['PDB ID','Segments','Length','Residues']
data.append(header)
for ifile in glob.glob('*.pdb'):
    print ifile
    ipath, iname, itype=filename(ifile)
    u=MD.Universe(ifile)
    #chain=u.selectAtoms('segid A')
    for seg in u.segments:
        print "Segment: {}".format(seg.name)
        p=seg.selectAtoms('protein')
        ress=''.join([ aa[res] for res in p.resnames() if res in aa.keys()])

        data.append([iname,seg.name,len(ress),ress])

writecsv(data,'Summary.csv')


#Find smallest fragment present in most sequences. It's probably shorter than average sequene size, but must be long enough to fully overlap with the active center



#At first get protein sequences from PDB files
os.chdir('/home/povilas/Data/Projects/2014-Prot-Press/Turiai/CA/CAII_IT1')

rfile='../CAII-P00918_frag.fasta'
mstring=readseq(rfile)
sfile=mstring[0:4]+str(len(mstring))+'_summary.csv'

results=[]
rhead=['PDB ID','Segment','Match length','Alignment','Cropped']
results.append(rhead)
results.append(['Ref','Ref',len(mstring),mstring,mstring])
ifiles=glob.glob('*.pdb')
ifiles=['4kap.pdb']

for ifile in sorted(ifiles):
    if ifile not in ['1can.pdb','1cra.pdb','1cao.pdb']:
        ipath, iname, itype=filename(ifile)
        res=crop(ifile,mstring,aa)
        results.extend(res)

writecsv(results,sfile)



pdb='HHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASFK'
ref='-WGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLAHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASF-'

als=align(mstring,pdb)

matches,lk=longmatch(als[1],als[0])







