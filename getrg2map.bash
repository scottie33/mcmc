#!/bin/bash

if [ $# -lt 0 ]; then
	echo ' getrg2map.bash'
	#echo ' depending on: calcdis.py'
	exit -1
fi

cp rg2maprange.gpl range.gpl

for eafi in `ls rg2*plot.dat`; do
	echo " working... please do not cut me, your majesty."
	cp $eafi inputfile.dat
	#paste inputfile.dat entropy.dat | awk '{print $3*}'
	gnuplot < draw_map.gpl
	mv oxmap.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
	convert -rotate 90 ${eafi}.eps ${eafi}.png
	echo " file: [ ${eafi}.png ] created. ;-)"
done

exit 0
