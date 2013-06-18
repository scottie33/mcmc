#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 4: 
	print " Syntax: getOE.py filename colid valTemp" # [shift]"
  	exit(-1)

filename=sys.argv[1]
colid=int(sys.argv[2])
valTemp=float(sys.argv[3])
try:
	tempfp=open("paraentropy.gpl", 'r')
	print " loading information of column from: [paraentropy.gpl]"
except IOError:
	print " can not open file: [paraentropy.gpl]"
	exit(-1)

epsilon=1.0
while True:
	line=tempfp.readline().rstrip()
	if line:
		elements=line.split('=')
		print "",elements
		if len(elements)==2:
			if elements[0].lstrip()=="epsilon":
				epsilon=float(elements[1])
				break
	else:
		break
print " epsilon=%f" % (epsilon)
tempfp.close()

try:
	fp=open("entropy.dat", 'r')
	print " loading information of column from: [entropy.dat]"
except IOError:
	print " can not open file: [entropy.dat]"
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

print " entropy information loaded, now create the OE(T) ... "
fp.close()

listlen=len(ener)

PFZ=-1.7e+308

for i in range(0,listlen):
	#print entropy[i]-ener[i]/valTemp
	PFZ=log_cc(PFZ, entropy[i]-ener[i]/valTemp)

print " log(Z)=%f" % (PFZ)

try:
	oefp=open(filename, 'r')
	print " loading information of column from: [",filename,"]"
except IOError:
	print " can not open file: [",filename,"]"
	exit(-1)

oelist=[]

while True:
	line=oefp.readline().rstrip()
	if line:
		elements=line.split()
		if len(elements)>=colid:
			oelist.append(float(elements[colid-1]))
	else:
		break
oefp.close()

if len(oelist) != listlen:
	print "there are %d and %d in oe and entropy respectively, exit." % (len(oelist), listlen)
	exit(-1)

OET_p=-1.7e+308
OET_n=-1.7e+308

for i in range(0,listlen):
	if oelist[i]>0.0:
		OET_p=log_cc( OET_p, log(oelist[i])+entropy[i]-ener[i]/valTemp ) 
	elif oelist[i]<0.0:
		OET_n=log_cc( OET_n, log(-oelist[i])+entropy[i]-ener[i]/valTemp ) 
	else:
		continue
	#print oelist[i]
oelist=[]
ener=[]
entropy=[]

print " oet_p=",OET_p
print " oet_n=",OET_n
OET=exp(OET_p-PFZ)-exp(OET_n-PFZ)

tempstr="echo "+str(valTemp)+" "+str(OET)+" > OET.dat"

print " T=",valTemp,"OE=",OET
os.system(tempstr)

exit(0)


