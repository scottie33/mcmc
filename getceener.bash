#!/bin/bash

if [ $# -lt 3 ];then
	echo " cmd energyfilename templistfilename coeff4contacts"
	exit -1
fi

rm tempREC.dat
touch tempREC.dat

cat $2 | while read line; do
	python getOE.py $1 2 $line
	cat OET.dat >> tempREC.dat
done

filename=`echo $1 | sed -e "s/\.dat//g"`

mv tempREC.dat ce${filename}.dat

./derivative ce${filename}.dat dece${filename}.dat

echo "fname='ce${filename}.dat'" > temp.gpl
echo "defname='dece${filename}.dat'" >> temp.gpl
echo "coeff=$3" >> temp.gpl

gnuplot < draw_ceener.gpl

mv ceener.eps ce${filename}.eps

exit 0