#!/bin/bash

runtimes=2
numprocs=32
host=hostfile2
for((i=0;i<${runtimes};i++));do
   #mpirun -np $numprocs ./mcmc 
   mpirun -np $numprocs --hostfile $host ./mcmc 
done

