#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 2: 
	print " Syntax: getPE.py valTemp [idxCE]" # [shift]"
	print " in _temperaturelist.pls T[idxCE]=valTemp, this option is for comparing CE v.s. MCE"
  	exit(-1)

valTemp=float(sys.argv[1])

idxCE=0
if len(sys.argv) > 2: 
	idxCE=int(sys.argv[2])


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
	pefp2=open("probability_each.dat", 'r')
	print " loading information of column from: [probability_each.dat]"
except IOError:
	print " can not open file: [probability_each.dat]"
	exit(-1)

pelist2=[]
pemax2=0.0
while True:
	line=pefp2.readline().rstrip()
	if line:
		elements=line.split()
		if idxCE>0 and len(elements)>(idxCE):
			petemp=float(elements[idxCE])
			pelist2.append(petemp)
			if petemp>pemax2:
				pemax2=petemp
	else:
		break
if len(pelist2)!=0: 
	print " probabilty_each information loaded. "
pefp2.close()


pemax1=0.0
pelist1=[]
for i in range(0,listlen):
	petemp=exp(entropy[i]-ener[i]/valTemp-PFZ)
	pelist1.append( petemp )
	if petemp>pemax1:
		pemax1=petemp

try:
	pefp=open("PE.dat", 'w')
	print " loading information of column from: [PE.dat]"
except IOError:
	print " can not open file: [PE.dat]"
	exit(-1)

coeff=pemax2/pemax1

for i in range(0,listlen):
	if len(pelist2)==0:
		print >> pefp, "%f %f %d %f" % (ener[i], pelist1[i], 0, ener[i]-valTemp*entropy[i])
	else:
		print >> pefp, "%f %f %f %f" % (ener[i], pelist1[i], pelist2[i]/coeff, ener[i]-valTemp*entropy[i])

pefp.close()

print " y range: [0.0,%f]" % (pemax1)
pemax1=1.05*pemax1

print " y range: [0.0,%f]" % (pemax1)
os.system("cp paraentropy.gpl parape.gpl")
tempstr="echo valTemp="+str(valTemp)+" >> parape.gpl"
os.system(tempstr)
tempstr="echo smin=0.0 >> parape.gpl"
os.system(tempstr)
tempstr="echo smax="+str(pemax1)+" >> parape.gpl"
os.system(tempstr)
tempstr="echo idxCE="+str(idxCE)+" >> parape.gpl"
os.system(tempstr)
os.system("gnuplot < draw_PE.gpl")
tempstr="mv PE.eps PE_"+str(valTemp)+".eps"
os.system(tempstr)

print " check [PE_%f.eps]" % (valTemp)

