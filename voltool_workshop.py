#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
"""
voltool.py

Created by Povilas Norvaisas on 2015-12-6.

"""

from voltool import *



os.chdir('/home/povilas/Data/Projects/2014-Prot-Press/Turiai/CA/CAII_screen')
selected=glob.glob('*.pdb')

os.chdir('/home/povilas/Data/Projects/2014-Prot-Press/Turiai/CA/CAII_IT1/249_aligned/Prep/New_alignment')

clean=glob.glob('*.pdb')

selection=[]
fl=open('Selection.txt','a')
for f in clean:
    if f in selected:
        shutil.copy(f,'Selected/'+f)
        fl.write(f+'\n')
fl.close()
