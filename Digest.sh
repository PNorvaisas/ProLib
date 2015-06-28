#!/bin/bash

for pdb in *.pdb; do
  Protease.py $pdb
  cd ${pdb%.*}
  voltool.py -i .pdb
  cat Summary.csv >> ../Summary_all.csv
  cd ..

done
