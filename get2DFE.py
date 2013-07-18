#!/usr/bin/python

import sys,os
from math import exp,log,pow

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))
cutoff=1000.0

if len(sys.argv) < 3: 
	print " Syntax: get2DFE.py filename temperature cutoff" # [shift]"
	print " to get: "
	print "    Ei: rg0 rg1 rg2 ... rgN [\enter] ... "
  	exit(-1)

filename=sys.argv[1]
valTemp=float(sys.argv[2]) # k_BT
beta=1.0/valTemp
cutoff=float(sys.argv[3])
print " temp=",valTemp,"beta=",beta

try:
	fp=open("entropy.dat", 'r')
	print " loading information of column from: [ entropy.dat ]"
except IOError:
	print " can not open file: [ entropy.dat ]"
	exit(-1)

ener=[]
entropy=[]

while True:
	line=fp.readline().rstrip()
	if line:
		elements=line.split()
		if len(elements)==2:
			ener.append(float(elements[0]))
			entropy.append(float(elements[1]))
	else:
		break

print " entropy information loaded, now create the normalized file ... "
fp.close()

listlen=len(ener)
print " listlen=",listlen
PFZ=-1.7e+308
for i in range(0,listlen):
	PFZ=log_cc(PFZ, entropy[i]-ener[i]*beta)
print " log(Z)=%f" % (PFZ)

pelist=[]
for i in range(0,listlen):
	pelist.append(exp(entropy[i]-ener[i]*beta-PFZ))
	
print " pelist created for the later calculation ... " # % (PFZ)
print " listlen=",listlen


try:
	oefp=open(filename, 'r')
	print " loading information of column from: [ "+filename+" ]"
except IOError:
	print " can not open file: [",filename,"]"
	exit(-1)

tfname=filename+".dat"
try:
	outfp=open(tfname, 'w')
	print " writing information into: [ "+tfname+" ]"
except IOError:
	print " can not open file: [ "+tfname+" ]"
	exit(-1)



i=0
x=0
maxfe=0.0
while True:
	line=oefp.readline()
	if line:
		line=line.rstrip()
		elements=line.split()

		if len(elements)==3:
			tempnum=float(elements[2])
			if tempnum > 0.0 :
				if pelist[i] > 0.0:
					print >>outfp, "%s %s %f" % ( elements[0], elements[1], -log(tempnum*pelist[i])*valTemp )
				else:
					print >>outfp, "%s %s %f" % ( elements[0], elements[1], cutoff )
			else:
				print >>outfp, "%s %s %f" % ( elements[0], elements[1], cutoff )
		if len(elements)==0:
			#print i
			print >>outfp, ""
			#i=0
			#x=x+1
			i=i+1
	else:
		break
print "",i,"data series loaded." 

ener=[]
pelist=[]
entropy=[]

oefp.close()
outfp.close()

exit(0)


