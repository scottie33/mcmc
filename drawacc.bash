#!/bin/bash

if [ $# -lt 1 ];then
   echo " cmd numrep-1"
   exit -1
fi

cat accept_ratio.dat | tail -n $1 | sed -e "s/accepted\ ratio: \[//g" | sed -e "s/\ ,//g" > acc_last.dat

gnuplot < draw_acc.gpl

