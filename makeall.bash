#!/bin/bash

make clean

make
make -f Makefile_muca
make -f Makefile_wham
make -f Makefile_newpdb

g++ smoothdata.cpp -o smoothdata
g++ derivative.cpp -o derivative
g++ integration.cpp -o integration 

make clean

