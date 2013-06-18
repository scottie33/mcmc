#!/usr/bin/python

import sys,os
from math import exp,log

if len(sys.argv) < 2: 
	print " Syntax: getPE.py valTemp" # [shift]"
  	exit(-1)
valTemp=float(sys.argv[1])

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
			ener.append(float(elements[0])/epsilon)
			entropy.append(float(elements[1]))
	else:
		break

print " entropy information loaded, now create the PE(T) ... "
fp.close()

listlen=len(ener)

try:
	pefp=open("PE.dat", 'w')
	print " loading information of column from: [PE.dat]"
except IOError:
	print " can not open file: [PE.dat]"
	exit(-1)

pelist=[]
pemin=0.0
pemax=0.0
for i in range(0,listlen):
	#print entropy[i]
	#print -ener[i]/valTemp
	if entropy[i]<-1.0e308:
		temppe=""
	else:
		temppe=entropy[i]-ener[i]/valTemp
	pelist.append(temppe)
	if temppe!="" and temppe<pemin:
		pemin=temppe
print " log(pe) min = %f" % (pemin)

for i in range(0,listlen):
	if pelist[i]!="": # and pelist[i]>0.0:
		pelist[i]=pelist[i]-pemin
		if pelist[i]>pemax:
			pemax=pelist[i]
pemin=0
print " y range: [%f,%f]" % (pemin, pemax)

threshold=16
tempminus=pemax-threshold
if tempminus>0.0:
	for i in range(0,listlen):
		if pelist[i]!="": # and pelist[i]>0.0:
			pelist[i]=pelist[i]-tempminus
print " y range: [%f,%f]" % (pemin, pemax)

for i in range(0,listlen):
	if pelist[i]!="": # and pelist[i]>0.0:
		#print " log(pe)=%f" % (temppe)
		print >> pefp, "%f %f" % (ener[i], exp(pelist[i]))
	else:
		print >> pefp, "%f %s" % (ener[i], pelist[i])
	
pefp.close()

print " y range: [%f,%f]" % (pemin, pemax)
pemax=exp(threshold)
distc=0.05*pemax
pemax=distc+pemax
pemin=pemin-distc
print " y range: [%f,%f]" % (pemin, pemax)
os.system("cp paraentropy.gpl parape.gpl")
tempstr="echo valTemp="+str(valTemp)+" >> parape.gpl"
os.system(tempstr)
tempstr="echo smin="+str(pemin)+" >> parape.gpl"
os.system(tempstr)
tempstr="echo smax="+str(pemax)+" >> parape.gpl"
os.system(tempstr)
os.system("gnuplot < draw_PE.gpl")
tempstr="mv PE.eps PE_"+str(valTemp)+".eps"
os.system(tempstr)

print " check [PE_%f.eps]" % (valTemp)





