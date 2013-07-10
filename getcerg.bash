#!/bin/bash

if [ $# -lt 2 ];then
	echo " cmd templistfilename chainid"
	exit -1
fi

cid=$2

tempname="RG2"
if [ $cid -lt 10 ]; then
	tempname="RG200${cid}.dat"
elif [ $cid -lt 100 ]; then
	tempname="RG20${cid}.dat"
else 
	tempname="RG2${cid}.dat"
fi
echo " now $tempname"
#cat $tempname | awk '{print $1,$2+$3+$4}' > tempRG2.dat
cp $tempname tempRG2.dat

rm tempREC.dat

touch tempREC.dat

cat $1 | while read line; do
	#getOE filename colid valTemp
	python getOE.py tempRG2.dat 2 $line
	cat OET.dat >> tempREC.dat
done

gnuplot < draw_cerg2.gpl
gnuplot < draw_cerg2_unlog.gpl

exit 0
