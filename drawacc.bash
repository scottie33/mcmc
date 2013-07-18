#!/bin/bash

if [ $# -lt 1 ];then
   echo " cmd numrep-1"
   exit -1
fi
let allnum=$1+1
#echo "$allnum"

cat accept_ratio.dat | tail -n $1 | sed -e "s/accepted\ ratio: \[//g" | sed -e "s/\ ,//g" > acc_last.dat

cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed "s/\]//g" | awk '{print $1,$5}' | head -n $1 | sort > temp1.dat
cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed "s/\]//g" | awk '{print $1,$5}' | tail -n 1 | sed -e "s/succ/000 succ/g"> temp2.dat

cat temp1.dat temp2.dat > accreptemp.dat
paste _temperaturelist.pls accreptemp.dat > accrep.dat
rm accreptemp.dat
rm temp1.dat 
rm temp2.dat 

#cat all.log | grep succ | tail -n $allnum | sed -e "s/@proc\[//g" | sed "s/\]//g" | awk '{print $1,$5}' | sort

gnuplot < draw_acc.gpl

