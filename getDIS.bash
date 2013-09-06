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
			mv tempCOM.dat ${tempname}_com.dat
			echo " file: [ ${tempname}_com.dat ] created. "
			#cp ${tempname}_com.dat distance.dat 
			#gnuplot < draw_distance.gpl
			#mv distance.eps ${tempname}_com.eps
			#echo " file: [ $tempname.eps ] created. ;-)"
		fi
	done
done

for eafi in `ls RG2*.dat`; do
	cp $eafi tempRG2.dat
	gnuplot < draw_RG2.gpl
	mv tempRG2.eps ${eafi}.eps
	mv tempFIC.eps ${eafi}_FIC.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
	echo " file: [ ${eafi}_FIC.eps ] created. ;-)"
done

for eafi in `ls EnerBF*.dat`; do
	cp $eafi tempEnerBF.dat
	gnuplot < draw_EBF.gpl
	mv tempEnerBF.eps ${eafi}.eps
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

for eafi in `ls DIS*.dat`; do
	if [ $eafi != "distance.dat" ]; then
		cp $eafi distance.dat 
		gnuplot < draw_distance2.gpl
		mv distance.eps ${eafi}.eps
		echo " file: [ ${eafi}.eps ] created. ;-)"
	fi
done

for eafi in `ls CN*.dat`; do
	cp $eafi contactnumber.dat 
	gnuplot < draw_CN.gpl
	mv contactnumber.eps ${eafi}.eps
	echo " file: [ ${eafi}.eps ] created. ;-)"
done

exit 0
