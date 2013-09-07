#!/bin/bash

runtimes=2
numprocs=64
host=hostfile
for((i=0;i<${runtimes};i++));do
   #mpirun -np $numprocs ./mcmc 
   mpirun -np $numprocs --hostfile $host ./mcmc 
done

./wham 25 1.0 1.0 1.0 16

