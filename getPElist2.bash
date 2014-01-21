#!/bin/bash

if [ $# -ne 6 ]; then
	echo " cmd epsilon Natoms xmin xmax fromid toid "
	exit -1
fi

rm PE*.dat

cp drawpelist.gpl drawpelistt.gpl



epsilon=$1
shift
Natoms=$1
shift

echo "set xrange [$1:$2]" >> drawpelistt.gpl
echo "set xtics 0.5" >> drawpelistt.gpl
shift 
shift

totnol=`cat _temperaturelist.pls | wc -l`
echo " there are $totnol lines in _temperaturelist.pls."
cutnol=$(($2-$1+1))
echo " we will get $cutnol lines from $1 to $2."
head -n $2 _temperaturelist.pls | tail -n $cutnol > newtemplist.pls

awk '{for(n=2;n<=NF;n++)t[n]+=$n}END{for(n=2;n<=NF;n++) printf t[n]"\n";}' probability_each.dat > temptotnum.dat

num=$1
ttnum=1
for i in `cat newtemplist.pls`; do
	echo "t=$i"
	python getPE.py $i $epsilon $Natoms #> templog
	mv PE.dat PE_$num.dat
	tempnum=0
	let tempnum=$num+1
	
	tttempn=0
	let tttempn=`cat temptotnum.dat | head -n $num | tail -n 1`
	echo " tttempncol=$tttempn"
	
	if [ $num -eq $1 ]; then
		echo -n "plot 'probability_each.dat' u (\$1/epsilon/numatoms):(\$${tempnum}/${tttempn}) w l lt 1 lw 4 lc rgb 'gray' title ''," >> drawpelistt.gpl
		echo -n " 'PE_$num.dat' u (\$1/epsilon/numatoms):2  w l lt $ttnum lw 1  title '{/Times-Italic T=$i}'" >> drawpelistt.gpl
	else
		echo -n " 'probability_each.dat' u (\$1/epsilon/numatoms):(\$${tempnum}/${tttempn}) w l lt 1 lw 4 lc rgb 'gray' title ''," >> drawpelistt.gpl
		echo -n " 'PE_$num.dat' u (\$1/epsilon/numatoms):2  w l lt $ttnum lw 1  title '{/Times-Italic T=$i}'" >> drawpelistt.gpl
	fi
	if [ $num -ne $2 ]; then
		echo ",\\" >> drawpelistt.gpl
	else 
		echo -e "\n" >> drawpelistt.gpl
	fi
	let ttnum=$ttnum+1
	let num=$num+1
done

echo "unset multiplot" >> drawpelistt.gpl

gnuplot < drawpelistt.gpl
echo " now you can check out your [ pelist.eps ]."

exit 0