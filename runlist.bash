#!/bin/bash

runtimes=2
numprocs=64
host=hostfile

startid=1 
endid=1

echo " mcmc proc ... "
for((i=0;i<${runtimes};i++));do
   #mpirun -np $numprocs ./mcmc
   mpirun -np $numprocs --hostfile $host ./mcmc
done

if [ $runtimes -gt 0 ]; then
	echo " wham ... "
	./wham 25 1.0 1.0 1.0 55
fi

echo " muca proc ... "
for((i=${startid};i<=${endid};i++)); do
	if [ $i -eq 1 ]; then
		if [ -e dirmuca ]; then
			#then echo "yes"; fi
			rm -fr dirmuca
		fi
		echo " muca in ./dirmuca ... "
		mkdir dirmuca
		if [ -e dirmuca ]; then
			echo "yes, dirmuca built!";
		else
			echo "can not create ./dirmuca/, exit ... "
			exit -1
		fi
		cp beta.dat beta.00.dat
		cp entropy.dat entropy.00.dat
		cd dirmuca
		cp -f ../_config.in ./; cp -f ../pnc.* ./;       cp -f ../_temperaturelist.pls ./
		cp -f ../hostfile ./;   cp -f ../_ffpara.pls ./; cp -f ../beta.dat ./
		cp -f ../muca ./;       cp -f ../confT???.* ./;   cp -f ../parabeta.gpl ./
		cp -f ../draw_entropy.gpl ./; cp -f ../draw_beta.gpl ./;  cp -f ../draw.gpl ./
		cp -f ../draw_prob.gpl ./; cp -f ../pdb2psf.py ./; cp -f ../smoothdata ./
		cp -f ../paraentropy.gpl ./;
		cd ..
	fi
	cd dirmuca
	# mpirun -np $numprocs ./muca
	mpirun -np $numprocs --hostfile $host ./muca
	if [ $i -lt 10 ]; then
		cp beta.dat beta.0$i.dat
		cp entropy.dat entropy.0$i.dat
	else
		cp beta.dat beta.$i.dat
		cp entropy.dat entropy.$i.dat
	fi
	cd ..
done



