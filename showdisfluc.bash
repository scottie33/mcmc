#!/bin/bash


if [ $# -lt 7 ]; then
	echo "cli: cmd psffile dir index rmin rmax Nrbin bonded? [totNF]"
	exit -1
fi

tstring="resid 1"

dir=$2
index=$3
rmin=$4
rmax=$5
Nrbin=$6

bonded=$7


bmin=0.27
bmax=1.73

dismin=-0.01
dismax=0.03

totNF=10000000
if [ $# -ge 8 ]; then
	let totNF=$8
fi

if [ $bonded -eq 1 ]; then
	cp bonddistance.tcl tbonddistance.tcl
else
	cp nonbonddistance.tcl tbonddistance.tcl
fi


echo "set psffilename \"$1\"" >> tbonddistance.tcl
echo "puts \$psffilename" >> tbonddistance.tcl
echo "mol new \$psffilename type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> tbonddistance.tcl
tempnumer=0
for eachpdb in `find ${dir} -name "*${index}.pdb"`; do	
	echo "mol addfile $eachpdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> tbonddistance.tcl
	let tempnumer=$tempnumer+1
	if [ $tempnumer -ge $totNF ]; then 
		break
	fi
done
echo "bonddistance \"$tstring\" dis.dat disdis.dat ${rmin} ${rmax} ${Nrbin}" >> tbonddistance.tcl
echo "quit" >> tbonddistance.tcl

vmd -dispdev text -e tbonddistance.tcl

if [ $bonded -eq 0 ]; then
	mv factors.dat factors_${index}.dat
fi

echo "imax=`tail -n 1 dis.dat | awk '{print $1}'`" > indexrange.gpl
echo "imin=`head -n 1 dis.dat | awk '{print $1}'`" >> indexrange.gpl
echo "bmin=$bmin" >> indexrange.gpl
echo "bmax=$bmax" >> indexrange.gpl
echo "rmin=${rmin}" >> indexrange.gpl
echo "rmax=${rmax}" >> indexrange.gpl
echo "dismin=$dismin" >> indexrange.gpl
echo "dismax=$dismax" >> indexrange.gpl


gnuplot < draw_bondfluc.gpl

mv bondeddis.dat bondeddis_${index}.dat 
mv dis.dat dis_${index}.dat
mv disdis.dat disdis_${index}.dat

mv bondfluctuation.eps bondfluc_${index}.eps
mv bonddistribution.eps bonddis_${index}.eps


