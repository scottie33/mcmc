#!/usr/bin/python

import sys,os

if len(sys.argv)<3:
	print " cmd Thigh Tlow NREM[10]"
	exit(-1)

Thigh=float(sys.argv[1])
Tlow=float(sys.argv[2])
if Thigh<Tlow:
	print " ERROR: [ Thigh < Tlow ]"
	exit(-1)

NREM=10
if len(sys.argv)>3:
	NREM=int(sys.argv[3])

print " T[%f:%f] into %d replicas" % (Tlow, Thigh, NREM)

for i in range(0,NREM):
	print Thigh/pow(Thigh/Tlow, float(i)/float(NREM-1))
Boltzmann=1.38e-23*0.23901*0.001*6.02e23 # kcal/k/mol
print "\n T for real:"
for i in range(0,NREM):
	print Thigh/Boltzmann/pow(Thigh/Tlow, float(i)/float(NREM-1))
exit(0)
