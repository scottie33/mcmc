#!/bin/bash


if [ $# -lt 4 ]; then
	echo ' getrg2map.bash filename templistfile xangle zangle'
	echo '     ( xangle=0 zangle=0 for map )'
	exit -1
fi

python mce2cemap.py $1

cp rg2maprange.gpl rangetemp.gpl
if [ $3 -eq 0 -a $4 -eq 0 ]; then
	echo "set pm3d map" >> rangetemp.gpl
else
	echo "set view $1,$2" >> rangetemp.gpl
fi
echo "xmin=`tail -n 1 $2`" >> rangetemp.gpl
echo "xmax=`head -n 1 $2`" >> rangetemp.gpl


num=0

filename=`echo "inputfile.dat"`
rm $filename
touch $filename

cat $2 | while read temp; do
	num=0
	cat oxlist | while read ox; do
		echo " doing $temp $ox"
		let num=$num+1
		python getOE.py $1.dat $num $temp 1> templog 
		echo "$temp $ox `cut -d \" \" -f 2 OET.dat`" >> $filename
	done
	echo " " >> $filename
done

gnuplot < draw_cemap.gpl
mv oxmap.eps $1_CE.eps
convert -rotate 90 $1_CE.eps $1_CE.png
echo " check out your [ $1_CE.png ]. "

exit 0
