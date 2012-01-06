#!/bin/bash

PKAL=/home/whitead/pdb2pqr/pdb2pqr.py
GERS=$HOME/Documents/ProteinSurfaces/lib/gerstein-surf
DSSP=$HOME/Documents/ProteinSurfaces/lib/dssp

#Local Directory for Wenjun Huang's Computer:
#PKAL=/home/wenjunh/pdb2pqr/pdb2pqr.py
#GERS=/home/wenjunh/Documents/ProteinSurfaces/calc-surface
#DSSP=/home/wenjunh/DSSP/dssp-2

for i in `ls *.gz`; 
  do gunzip $i; 
done;

for i in `ls | sed -n 's/\(.*\)\.pdb/\1/p'`; do 
   mkdir $i; 
   mv $i.pdb $i; 
done;

if [[ !(-e chargeList) ]]; then
  echo "pdb charge surface" > chargeList
fi;

for i in `ls`; do 
  if [[ !(-e $i/$i.pka) && $i != "chargeList" ]]; then
    $DSSP -i $i/$i.pdb -o $i/$i.dssp
    if [[ !(-e $i/$i.surf) && !(-e $i/$i.wwater.surf)]]; then
      cat $i/$i.pdb | grep -v " HOH " > $i/$i.temp.pdb
      $GERS -i $i/$i.temp.pdb -o $i/$i.surf;
      $GERS -i $i/$i.pdb -o $i/$i.wwater.surf; 
    fi;
    SF=`$HOME/Documents/ProteinSurfaces/psurf/analyzeSA.py $i/$i.surf 67 $i/sadata;`
    /usr/bin/python $PKAL --noopt --ff=PARSE -v --with-ph=7.0 $i/$i.pdb $i/$i.pka;
    CH=`sed -n 's/REMARK   6 Total charge on this protein: \([-0-9\.]*\) e/\1/p' $i/$i.pka`;
    echo "$i $CH $SF" >> chargeList;
  else
    echo "Skipping pka calculation for $i";
    $DSSP -i $i/$i.pdb -o $i/$i.dssp
    if [[ !(-e $i/$i.surf) && !(-e $i/$i.wwater.surf)]]; then
      cat $i/$i.pdb | grep -v " HOH " > $i/$i.temp.pdb
      $GERS -i $i/$i.temp.pdb -o $i/$i.surf;
      $GERS -i $i/$i.pdb -o $i/$i.wwater.surf;
    fi;
    SF=`$HOME/Documents/ProteinSurfaces/psurf/analyzeSA.py $i/$i.surf 67 $i/sadata;`
    CH=`sed -n 's/REMARK   6 Total charge on this protein: \([-0-9\.]*\) e/\1/p' $i/$i.pka`;
    echo "$i $CH $SF" >> chargeList;
  fi;
done;
