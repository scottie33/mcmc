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

mv tempREC.dat cerg2.dat

./derivative cerg2.dat decerg2.dat

tempT=`cat decerg2.dat | sort -k2,2n | tail -n 1 | awk '{print $1}'`
tempC=`cat decerg2.dat | sort -k2,2n | tail -n 1 | awk '{print $2}'`

echo "set arrow 1 from $tempT, second 0.0 to $tempT, second $tempC nohead lt 1 lw 0.1 lc rgb 'gray'" > tempcerg.gpl
echo "set label '{/Times-Italic T_{/Symbol Q}=$tempT}' at $tempT+0.5, second $tempC-0.5 left" >> tempcerg.gpl

gnuplot < draw_cerg2.gpl

exit 0

