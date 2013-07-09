#!/bin/bash

runtimes=2
numprocs=32

for((i=0;i<${runtimes};i++));do
   mpirun -np $numprocs ./mcmc 
done

