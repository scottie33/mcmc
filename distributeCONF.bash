#!/bin/bash

if [ $# -lt 2 ]; then
	echo " usage: cmd totnum disnum "
	echo " e.g.: distributeCONF 64 4 "
	exit -1
fi

if [ $2 -gt $1 ]; then
	echo " disnum >= $1, error!"
	exit -1
fi

j=1
for((i=1;i<=$1;i++));do
	let j=$j+1
	if [ $j -gt $2 ]; then
		let j=1
	fi
	if [ $i -lt 10 ]; then
		if [ $j -lt 10 ]; then
			cp confT00${j}.pdb confT00${i}.pdb
		elif [ $j -lt 100 ]; then
			cp confT0${j}.pdb  confT00${i}.pdb
		else 
			cp confT${j}.pdb   confT00${i}.pdb
		fi
	elif [ $i -lt 100 ]; then
		if [ $j -lt 10 ]; then
			cp confT00${j}.pdb confT0${i}.pdb
		elif [ $j -lt 100 ]; then
			cp confT0${j}.pdb  confT0${i}.pdb
		else 
			cp confT${j}.pdb   confT0${i}.pdb
		fi
	else
		if [ $j -lt 10 ]; then
			cp confT00${j}.pdb confT${i}.pdb
		elif [ $j -lt 100 ]; then
			cp confT0${j}.pdb  confT${i}.pdb
		else 
			cp confT${j}.pdb   confT${i}.pdb
		fi
	fi
done

exit 0