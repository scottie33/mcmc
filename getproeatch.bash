#!/bin/bash
if [ $# -lt 1 ]; then
	echo " cmd numofrep"
	exit -1
fi

for((i=1;i<=$1;i++));do
	let ncol=$i+1
	cat probability_each.dat | awk '{print ${$ncol}}' > pro$i.dat
done

