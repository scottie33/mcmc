#!/bin/bash

if [ $# -lt 2 ]; then
	echo " cmd filename bsnum"
	exit -1
fi

echo $1 > rg2_file.lst

echo "#start log." > data_log.lst
for i in 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 40 50; do
	echo " ./ana_rg2 rg2_file.lst ${i} $2"
	./ana_rg2 rg2_file.lst ${i} $2
	cat rg2_data.lst >> data_log.lst
done

echo " please check [ data_log.lst ] "

gnuplot showbootstrap.gpl

echo " please check [ showbootstrap.eps ] "
#open showbootstrap.eps
