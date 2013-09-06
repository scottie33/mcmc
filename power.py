#!/usr/bin/python

from math import pow
from os import sys#

for i in range(int(sys.argv[1])+1):
	print i,float(sys.argv[2])*pow(float(sys.argv[3]), float(i)/float(sys.argv[1])) 