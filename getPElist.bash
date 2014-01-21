#!/bin/bash

if [ $# -lt 5 ]; then
	echo " cmd epsilon Natoms xmin xmax t1 [t2 [t3 [t4 [t5]]]] ... "
	exit -1
fi

rm PE*.dat

cp drawpelist.gpl drawpelistt.gpl



epsilon=$1
shift
Natoms=$1
shift

echo "set xrange [$1:$2]" >> drawpelistt.gpl
shift 
shift

num=1
for i in $@; do
	echo "t=$i"
	python getPE.py $i $epsilon $Natoms #> templog
	mv PE.dat PE_$num.dat
	
	if [ $num -eq 1 ]; then
		echo -n "plot 'PE_$num.dat' u (\$1/epsilon/numatoms):2  w l lt $num lw 1  title '{/Times-Italic T=$i}'" >> drawpelistt.gpl
	else
		echo -n " 'PE_$num.dat' u (\$1/epsilon/numatoms):2  w l lt $num lw 1  title '{/Times-Italic T=$i}'" >> drawpelistt.gpl
	fi
	if [ $num -ne $# ]; then
		echo ",\\" >> drawpelistt.gpl
	else 
		echo -e "\n" >> drawpelistt.gpl
	fi
	let num=$num+1
done

echo "unset multiplot" >> drawpelistt.gpl

gnuplot < drawpelistt.gpl
echo " now you can check out your [ pelist.eps ]."

exit 0