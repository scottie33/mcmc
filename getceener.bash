#!/bin/bash

if [ $# -lt 3 ];then
	echo " cmd energyfilename templistfilename coeff4contacts"
	exit -1
fi

rm tempREC.dat
touch tempREC.dat

cat $2 | while read line; do
	#getOE filename colid valTemp
	python getOE.py $1 2 $line
	cat OET.dat >> tempREC.dat
done

echo "coeff=$3" > temp.gpl

gnuplot < draw_ceener.gpl

exit 0