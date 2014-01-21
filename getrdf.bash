#!/bin/bash


if [ $# -lt 7 ]; then
	echo "cli: cmd psffile dir index deltar radiimax PBCorNot com? [comrlower] [totNF]"
	exit -1
fi

#proc getrdf {seltext1 seltex2 deltar radiimax PBCorNot} {
#set gr0 [measure gofr $sel1 $sel2 delta $deltar rmax $radiimax usepbc 0 selupdate 1 first 0 last -1 step 1]

tstring1="resid 1" ;#"index 27"
tstring2="resid 1" ;#"not index 27"

dir=$2
index=$3
deltar=$4
radiimax=$5
PBCorNot=$6
COMorNot=$7

comrlower=1
if [ $# -ge 8 ]; then
	comrlower=$8
fi

bmin=0.27
bmax=1.73

totNF=10000000
if [ $# -ge 9 ]; then
	let totNF=$9
fi

cp rdf.tcl rdftemp.tcl


echo "set psffilename \"$1\"" >> rdftemp.tcl
echo "puts \$psffilename" >> rdftemp.tcl
echo "mol new \$psffilename type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> rdftemp.tcl
tempnumer=0
for eachpdb in `find ${dir} -name "*${index}.pdb"`; do	
	echo "mol addfile $eachpdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> rdftemp.tcl
	let tempnumer=$tempnumer+1
	if [ $tempnumer -ge $totNF ]; then 
		break
	fi
done
if [ $COMorNot -eq 1 ];then
	echo "getrdfcom \"$tstring1\" \"$tstring2\" $deltar $radiimax $PBCorNot $index $comrlower" >> rdftemp.tcl
else
	echo "getrdf \"$tstring1\" \"$tstring2\" $deltar $radiimax $PBCorNot $index" >> rdftemp.tcl
fi
echo "quit" >> rdftemp.tcl

vmd -dispdev text -e rdftemp.tcl

echo "bmin=$bmin" > indexrange.gpl
echo "bmax=$bmax" >> indexrange.gpl
echo "rmin=0.0" >> indexrange.gpl
echo "rmax=${radiimax}" >> indexrange.gpl
echo "fname='grdf_${index}.dat'" >> indexrange.gpl


gnuplot < draw_rdf.gpl

mv rdf.eps grdf_${index}.eps



