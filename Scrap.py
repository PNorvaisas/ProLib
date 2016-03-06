


pdb='HHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASFK'
ref='-WGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLAHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASF-'

als=align(cropseq,pdb)

matches,lk=longmatch(als[1],als[0])



odir=dircheck('249')

for r in results[1:]:
    if r[2]==249 and r[1]!='Ref':
        shutil.copy('{}_{}_{}-{}_cropped.pdb'.format(r[0],r[1],cropseq[0:4],str(len(cropseq))),os.path.join(odir,'{}_{}.pdb'.format(r[0],r[1])))



data=[]
header=['File','Segment','Has ZN','Only ligand','Has other ligands','Heteroatoms']
data.append(header)
allress=[]
maxress=0


trash=['ZN','GOL','DMS','HG','SO4','CL','BR','NA','IOD','MBO','MES','UNX']
os.chdir('/home/povilas/Data/Projects/2014-Prot-Press/Turiai/CA/CAII_IT1/249_aligned/Prep')
ifiles=glob.glob('*.pdb')
for ifile in ifiles:
    print ifile
    ipath, iname, itype=filename(ifile)
    u=MD.Universe(ifile)
    #chain=u.selectAtoms('segid A')
    for seg in u.segments:
        print "Segment: {}".format(seg.name)
        p=seg.selectAtoms('protein')
        hetatm=seg.selectAtoms('not protein')
        ress=hetatm.resnames()
        rez=[ i for i in list(ress) if i not in trash]
        haslig=str(len(rez)>0)
        if len(rez)>maxress:
            maxress=len(rez)
        onl=rez[0] if len(rez)==1 else ''
        data.append([iname,seg.name,str('ZN' in ress),onl,haslig,';'.join(['']*(11-len(rez))+sorted(rez))])
        allress.extend(ress)
        #data.append([iname,seg.name,len(ress),ress])

writecsv(data,'Ligands_all.csv')

for d in data[1:]:
    if d[2]=='True':
        shutil.copy('{}.pdb'.format(d[0]),'{}/{}.pdb'.format('Selected',d[0]))

#Counter generates table with number of occurences for each instance of ligand

occs=Counter(allress)
ocsum=[]
for k in occs.keys():
    ocsum.append([k,occs[k]])
writecsv(ocsum,'Occurences_all.csv')

all=list(set(allress))





