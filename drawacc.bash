#!/bin/bash

if [ $# -lt 4 ];then
   echo " cmd numrep-1 criterionT steplength scaleto"
   exit -1
fi
let allnum=$1+1
#echo "$allnum"

cat accept_ratio.dat | tail -n $1 | sed -e "s/accepted\ ratio: \[//g" | sed -e "s/\ ,//g" > acc_last.dat

cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed -e "s/\]//g" | awk '{print $1,$5}' | sort > accreptemp.dat
#cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed -e "s/\]//g" | awk '{print $1,$4}' | tail -n 1 | sed -e "s/succ/000 succ/g"> temp2.dat

#cat temp1.dat temp2.dat > accreptemp.dat
#cat _temperaturelist.pls | sort > temptemp.dat

echo "criterionT=$2" > accpara.gpl
#echo "steplength=$3" >> accpara.gpl
#echo "scaleto=$4" >> accpara.gpl

python power.py $1 $3 $4 > tempaccpow.dat


cat all.log | grep fail | tail -n $allnum | sed -e "s/@proc\[//g" | sed -e "s/\]//g" | awk '{print $1,$5}' | sort > filreptemp.dat
paste _temperaturelist.pls filreptemp.dat > filrep.dat


paste _temperaturelist.pls accreptemp.dat tempaccpow.dat filrep.dat > accrep.dat

#rm accreptemp.dat
#rm temp1.dat 
#rm temp2.dat 
echo "tmax=`cat _temperaturelist.pls | head -n 1`" > temptrange.gpl
echo "tmin=`cat _temperaturelist.pls | tail -n 1`" >> temptrange.gpl
#cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed "s/\]//g" | awk '{print $1,$5}' | sort

gnuplot < draw_acc.gpl

