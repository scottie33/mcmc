#!/bin/bash

if [ $# -lt 2 ]; then
	echo ' getDIS.bash halfboxsize coeff4contacts'
	echo ' depending on: calcdis.py'
	exit -1
fi

halflen=$1

echo " half boxsize is: $halflen"

echo "coeff=$2" > temp.gpl

for eafi in `ls COM*.dat`; do
	for eafj in `ls COM*.dat`; do  
		if [ $eafi != $eafj ]; then
			tempname=`echo dis$eafi$eafj | sed -e "s/\.dat//g" | sed -e "s/COM//g"`
			echo " creating $tempname"
			paste $eafi $eafj > tempCOM
			python calcdis.py tempCOM $halflen $halflen $halflen
			rm tempCOM
			mv tempCOM.dat ${tempname}.dat
			echo " file: [ $tempname.dat ] created. "
			cp ${tempname}.dat distance.dat 
			gnuplot < draw_distance.gpl
			mv distance.eps ${tempname}.eps
			echo " file: [ $tempname.eps ] created. ;-)"
			#cat tempCOM.dat | awk '{print $1, sqrt( ($2-$6)*($2-$6)+($3-$7)*($3-$7)+($4-$8)*($4-$8) )}' > $tempname
			#rm $tempname.dat
			#touch $tempname.dat
			#cat $tempname | while read line; do
			#	if [ $i > $halflen ]; then
			#		echo $i
			#		echo `expr $i-$halflen` >> $tempname.dat
			#	else
			#		echo $i >> $tempname.dat
			#	fi
			#done
			#rm $tempname
		fi
	done
done

for eafi in `ls RG2*.dat`; do
	cp $eafi tempRG2.dat
	gnuplot < draw_RG2.gpl
	mv tempRG2.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
done

for eafi in `ls EnerLJ*.dat`; do
	cp $eafi tempEnerLJ.dat
	gnuplot < draw_ELJ.gpl
	mv tempEnerLJ.eps ${eafi}.eps
	mv Contacts.eps contactN${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
	echo " file: [ contactN${eafi}.eps ] created. ;-)"
done

for eafi in `ls EnerAG*.dat`; do
	cp $eafi tempEnerAG.dat
	gnuplot < draw_EAG.gpl
	mv tempEnerAG.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
done

for eafi in `ls EnerDH*.dat`; do
	cp $eafi tempEnerDH.dat
	gnuplot < draw_EDH.gpl
	mv tempEnerDH.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
done

exit 0
