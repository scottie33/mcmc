#!/bin/bash

if [ $# -lt 2 ]; then
	echo " cmd: fromid nprocs"
	exit -1
fi

for((i=0;i<$2;i++)); do
	let tid=$i+$1
	let ttd=$i+1
	echo " from $tid to $ttd ... "
	if [ $tid -lt 10 ]; then
		if [ $ttd -lt 10 ]; then
			cp confT00$tid.pdb confT00$ttd.pdb
			cp confT00$tid.xyz confT00$ttd.xyz
		elif [ $ttd -lt 100 ]; then
			cp confT00$tid.pdb confT0$ttd.pdb
			cp confT00$tid.xyz confT0$ttd.xyz
		else
			cp confT00$tid.pdb confT$ttd.pdb
			cp confT00$tid.xyz confT$ttd.xyz
		fi
	elif [ $tid -lt 100 ]; then
		if [ $ttd -lt 10 ]; then
			cp confT0$tid.pdb confT00$ttd.pdb
			cp confT0$tid.xyz confT00$ttd.xyz
		elif [ $ttd -lt 100 ]; then
			cp confT0$tid.pdb confT0$ttd.pdb
			cp confT0$tid.xyz confT0$ttd.xyz
		else
			cp confT0$tid.pdb confT$ttd.pdb
			cp confT0$tid.xyz confT$ttd.xyz
		fi
	else
		if [ $ttd -lt 10 ]; then
			cp confT$tid.pdb confT00$ttd.pdb
			cp confT$tid.xyz confT00$ttd.xyz
		elif [ $ttd -lt 100 ]; then
			cp confT$tid.pdb confT0$ttd.pdb
			cp confT$tid.xyz confT0$ttd.xyz
		else
			cp confT$tid.pdb confT$ttd.pdb
			cp confT$tid.xyz confT$ttd.xyz
		fi
	fi
done
