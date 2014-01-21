#!/bin/bash


if [ $# -lt 7 ]; then
	echo ' get2Dcemap.bash filename temperature cutoff coeff xangle zangle 1or2'
	echo '     ( xangle=0 zangle=0 for map )'
	exit -1
fi

python get2DFE.py $1 $2 $3

cp rg2maprange.gpl rangetemp.gpl
echo "coeff=$4" >> rangetemp.gpl
echo "set palette gray negative" >>rangetemp.gpl
#echo "set palette rgbformulae 23,28,3" >> rangetemp.gpl
#echo "set palette rgbformulae 33,13,10" >> rangetemp.gpl
#echo "set palette rgbformulae 7,5,15" >> rangetemp.gpl
#echo "set palette defined ( 0 '#9FAFDF', 1 'blue', 2 '#4B0082', 3 'white' ) " >> rangetemp.gpl
echo "set palette defined ( 0 '#9FAFDF', 1 'red', 2 '#4B0082', 3 'white' ) " >> rangetemp.gpl
#echo "set palette model HSV rgbformulae 7,5,15" >>rangetemp.gpl

echo "unset grid " >> rangetemp.gpl
if [ $5 -eq 0 -a $6 -eq 0 ]; then
	#echo "set pm3d map" >> rangetemp.gpl
	#echo "unset pm3d" >> rangetemp.gpl
	echo "unset surface" >> rangetemp.gpl
	#echo "set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover" >> rangetemp.gpl
	echo "set isosamples 50,50" >> rangetemp.gpl
	#echo "set dgrid3d" >> rangetemp.gpl
	echo "set contour base" >> rangetemp.gpl
	echo "set cntrparam bspline" >> rangetemp.gpl
	echo "set cntrparam levels 50 " >> rangetemp.gpl
	#echo "set cntrparam levels incremental 0.01,0.25,5" >> rangetemp.gpl
	echo "set cntrparam levels discrete 2.785,2.85,3.0,3.25,3.6,3,78,3.95,4.5,5.25" >> rangetemp.gpl
	echo "set cntrparam order 5" >> rangetemp.gpl
fi

echo "xangle=$5" >> rangetemp.gpl
echo "zangle=$6" >> rangetemp.gpl

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
