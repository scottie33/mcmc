#!/bin/bash

if [ $# -lt 2 ];then
	echo " cmd templistfilename disfilename"
	exit -1
fi

rm tempREC.dat
touch tempREC.dat

cat $1 | while read line; do
	#getOE filename colid valTemp
	./getOE.py $2 2 $line
	cat OET.dat >> tempREC.dat
done

gnuplot < draw_cedis.gpl

exit 0