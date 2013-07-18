#!/bin/bash


if [ $# -lt 7 ]; then
	echo ' get2Dcemap.bash filename temperature cutoff coeff xangle zangle 1or2'
	echo '     ( xangle=0 zangle=0 for map )'
	exit -1
fi

python get2DFE.py $1 $2 $3

cp rg2maprange.gpl rangetemp.gpl
echo "coeff=$4" >> rangetemp.gpl
if [ $5 -eq 0 -a $6 -eq 0 ]; then
	echo "set pm3d map" >> rangetemp.gpl
	echo "set contour base" >> rangetemp.gpl
	echo "set cntrparam levels 50" >> rangetemp.gpl
else
	echo "set view $5,$6" >> rangetemp.gpl
fi

num=0

filename=`echo "inputfile.dat"`
cp $1.dat $filename

if [ $7 -eq 1 ]; then
	gnuplot < draw2Dcemap.gpl
fi
if [ $7 -eq 2 ]; then
	gnuplot < draw2Dcemap2.gpl
fi
mv oxmap.eps $1_2DFE.eps
convert -rotate 90 $1_2DFE.eps $1_2DFE.png
echo " check out your [ $1_2DFE.png ]. "

exit 0
