#!/bin/bash


if [ $# -lt 1 ];then
	echo "cmd sigma"
	exit -1
fi

echo " reading from file.lst then calc .dat, writing into data.lst."

#echo " log <-> gpl.log" > gpl.log
#for i in rg2 rh2; do 
for i in rg2 ; do 
	#echo $i
	if [ -f ${i}_file.lst ]; then rm ${i}_file.lst; fi
	if [ ! -f ${i}_file.lst ]; then
		echo " [ ${i}_file.lst ] does not exist, creating now...";
		if [ `ls ${i}*.dat | wc -l` -gt 0 ]; then 
			ls ${i}*.dat
			ls ${i}*.dat > ${i}_file.lst
		fi
	fi
	echo " calculation ..."
	./ana_rg2 ${i}_file.lst 1 2000
	echo " plotting into [ ${i}_data.eps ] ..."
	echo "fn='${i}_data'" > tempfn.gpl 
	echo "sigma=$1" >> tempfn.gpl 
 	echo "yname='{/Times-Italic R_{G}^2/{/Symbol s}^2}'" >> tempfn.gpl
	gnuplot rg2list.gpl  #| echo &> gpl.log
done


