#!/bin/bash

if [ $# -lt 2 ]; then
	echo ' getdismap.bash xangle zangle'
	echo '           ( xangle=0 zangle=0 for map )'
	exit -1
fi

cp dismaprange.gpl rangetemp.gpl
if [ $1 -eq 0 -a $2 -eq 0 ]; then
	echo "set pm3d map" >> rangetemp.gpl
else
	echo "set view $1,$2" >> rangetemp.gpl
fi

for eafi in `ls dis*plot.dat`; do
	echo " working... please do not cut me, your majesty."
	#python mcenorm.py $eafi #op1
	#cp $eafi.dat inputfile.dat #op1
	cp $eafi inputfile.dat #op2
	#paste inputfile.dat entropy.dat | awk '{print $3*}'
	gnuplot < draw_map2.gpl
	mv oxmap.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
	convert -rotate 90 ${eafi}.eps ${eafi}.png
	echo " file: [ ${eafi}.png ] created. ;-)"
done

exit 0
