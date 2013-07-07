#!/usr/bin/python

import sys,os
from math import exp,log

nradi=9.2
Nn=8

massload=0.154
rho=0.567

pmass=56.0
pradi=4.7
nmass=pow(nradi/pradi,3.0)*pmass

Mass=float(Nn)*nmass/massload
pMass=Mass*(1.0-massload)
pnum=pMass/float(pmass)

volume=Mass/rho
print " \n ===== \n"
print " there are %f NP in the system with mass %f each \n" % (Nn,nmass)
print " there should be [ %f ] monomers in the system <- mload %f & NumN %f" % (pnum, massload, Nn)
print " "
print " and the volume should be [ %f ] given the Rho %f fixed." % (volume, rho)
print " "
print " and the length of each dimension should be [ %f ]." % (pow(volume,1.0/3.0))
print " \n =====\n"
exit(0)