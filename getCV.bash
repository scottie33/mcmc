#!/bin/bash

if [ $# -lt 1 ];then 
	echo " syntax: cmd templistfilename"
	exit -1
fi

if [ ! -f "entropy.dat" ];then
	echo " no [entropy.dat] found, please do wham first."
	exit -1
fi


cat entropy.dat | awk '{print $1,$1*$1}' > squaredE.dat

if [ -f MCEE.dat ]; then
	rm MCEE.dat
fi

touch MCEE.dat

cat $1 | while read line; do
	./getOE.py squaredE.dat 1 $line
	cat OET.dat >> MCEE.dat
done

if [ -f MCEE2.dat ]; then
	rm MCEE2.dat
fi

touch MCEE2.dat

cat $1 | while read line; do
	./getOE.py squaredE.dat 2 $line
	cat OET.dat >> MCEE2.dat
done

paste MCEE.dat MCEE2.dat | awk '{print $1,($4-$2*$2)/$1/$1}' > MCECV.dat

if [ -f CECV.dat ]; then
	rm CECV.dat
fi

touch CECV.dat
tempnum=1
cat _temperaturelist.pls | while read line; do
	let tempnum=$tempnum+1
	./getcvfromprob.py probability_each.dat $tempnum $line
	cat OET.dat >> CECV.dat
done

gnuplot < draw_cv.gpl



exit 0