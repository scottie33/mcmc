#!/bin/bash

make clean

make
make -f Makefile_wham
make -f Makefile_newpdb

make clean

