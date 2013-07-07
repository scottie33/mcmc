#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 2: 
	print " Syntax: mcenorm.py filename" # [shift]"
  	exit(-1)

filename=sys.argv[1]


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

#PFZ=-1.7e+308
#for i in range(0,listlen):
#	PFZ=log_cc(PFZ, entropy[i]-ener[i]/valTemp)
#print " log(Z)=%f" % (PFZ)

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

i=0;
listrg=[]
while True:
	line=oefp.readline()
	if line:
		line=line.rstrip()
		elements=line.split()
		if len(elements)==3:
			if i==0:
				listrg.append(elements[1])
			print >>outfp, elements[2],
		if len(elements)==0:
			print >>outfp, ""
			i=i+1
	else:
		break

print "",i,"data series loaded."
print " the ox list is: ",listrg

lenrg=len(listrg)

oxfp=open("oxlist", 'w')
for i in range(0,lenrg):
	print >> oxfp, listrg[i]
oxfp.close()


#for i in range(0,templen):
	#valTemp=templist[i]
	#for j in range(0,lenrg):
		#valrg=listrg[j]
		#tempstr="python getOE.py "+tfname+" "+str(valTemp)
		#sprint tempstr
		#system("python getOE.py tfname valTemp")


ener=[]
listrg=[]
entropy=[]

oefp.close()
outfp.close()

exit(0)


