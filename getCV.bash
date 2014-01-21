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
	python getOE.py squaredE.dat 1 $line
	cat OET.dat >> MCEE.dat
done

#if [ -f MCEE2.dat ]; then
#	rm MCEE2.dat
#fi

#touch MCEE2.dat

#cat $1 | while read line; do
#	python getOE.py squaredE.dat 2 $line
#	cat OET.dat >> MCEE2.dat
#done

#paste MCEE.dat MCEE2.dat | awk '{print $1,($4-$2*$2)/$1/$1}' > MCECV.dat

#if [ -f CECV.dat ]; then
#	rm CECV.dat
#fi

#touch CECV.dat
#tempnum=1
#cat _temperaturelist.pls | while read line; do
#	let tempnum=$tempnum+1
#	python getcvfromprob.py probability_each.dat $tempnum $line
#	cat OET.dat >> CECV.dat
#done

./derivative MCEE.dat deMCEE.dat

#echo "tmax=`cat $1 | head -n 1`" > temptrange.gpl
#echo "tmin=`cat $1 | tail -n 1`" >> temptrange.gpl

tempT=`cat deMCEE.dat | sort -k2,2n | tail -n 1 | awk '{print $1}'`
tempC=`cat deMCEE.dat | sort -k2,2n | tail -n 1 | awk '{print $2}'`

echo "set arrow 1 from $tempT, second 0.0 to $tempT, second $tempC nohead lt 1 lw 0.1 lc rgb 'gray'" > tempcv.gpl
echo "set label '{/Times-Italic T_{S-L}=$tempT}' at $tempT+0.5, second $tempC-0.5 left" >> tempcv.gpl

gnuplot < draw_cv.gpl

exit 0
