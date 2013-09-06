#!/bin/bash


if [ $# -lt 4 ]; then
	echo "cli: cmd psffile rmin rmax Nrbin [totNF]"
	exit -1
fi


totNF=10000000
if [ $# -ge 5 ]; then
	let totNF=$5
fi

cp bonddistance.tcl tbonddistance.tcl

tstring="resid 1"
bmin=0.27
bmax=1.73

dismin=-0.05
dismax=0.45

echo "set psffilename \"$1\"" >> tbonddistance.tcl
echo "puts \$psffilename" >> tbonddistance.tcl
echo "mol new \$psffilename type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> tbonddistance.tcl
tempnumer=0
for eachpdb in `ls confT*.pdb`; do	
	echo "mol addfile $eachpdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> tbonddistance.tcl
	let tempnumer=$tempnumer+1
	if [ $tempnumer -ge $totNF ]; then 
		break
	fi
done
echo "bonddistance \"$tstring\" dis.dat disdis.dat $2 $3 $4" >> tbonddistance.tcl
echo "quit" >> tbonddistance.tcl

vmd -dispdev text -e tbonddistance.tcl

echo "imax=`tail -n 1 dis.dat | awk '{print $1}'`" > indexrange.gpl
echo "imin=`head -n 1 dis.dat | awk '{print $1}'`" >> indexrange.gpl
echo "bmin=$bmin" >> indexrange.gpl
echo "bmax=$bmax" >> indexrange.gpl
echo "rmin=$2" >> indexrange.gpl
echo "rmax=$3" >> indexrange.gpl
echo "dismin=$dismin" >> indexrange.gpl
echo "dismax=$dismax" >> indexrange.gpl


gnuplot < draw_bondfluc.gpl

