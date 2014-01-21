#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 4: 
	print " Syntax: getlocalcv.py index1 index2 templist [energyfile]" # [shift]"
  	exit(-1)

index1=int(sys.argv[1])
index2=int(sys.argv[2])

filename=''
if len(sys.argv) > 4: 
	filename=sys.argv[4]


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

print " entropy information loaded."
fp.close()

if len(sys.argv) > 4 : 
	try:
		fp_ener=open(filename, 'r')
		print " loading information of column from: [",filename,"]"
	except IOError:
		print " can not open file: [",filename,"]"
		exit(-1)

	ener=[]

	while True:
		line=fp_ener.readline().rstrip()
		if line:
			elements=line.split()
			if len(elements)==2:
				ener.append(float(elements[1]))
		else:
			break

	print " energy from: [",filename,"] loaded."
	fp_ener.close()

Ener1=ener[index1]
Ener2=ener[index2]
Entropy1=entropy[index1]
Entropy2=entropy[index2]
try:
	fp2=open(sys.argv[3], 'r')
	print " loading information of column from: [",sys.argv[3],"]."
except IOError:
	print " can not open file: [",sys.argv[3],"."
	exit(-1)

templist=[]

while True:
	line=fp2.readline().rstrip()
	if line:
		elements=line.split()
		if len(elements)!=0:
			templist.append(float(elements[0]))
	else:
		break

print " templist loaded."
fp2.close()

listlen=len(ener)
templen=len(templist)

PFZ1list=[]
PFZ2list=[]

for i in range(0,templen):
	# print "Entropy1=",Entropy1
	# print "Ener1=",Ener1
	# print "Temp=",templist[i]
	# print "pfz=",Entropy1-Ener1/templist[i]
	# if Entropy1<-1.7e+308:
	# 	PFZ1list.append(0.0)
	# else:
	# 	PFZ1list.append(exp(Entropy1-Ener1/templist[i]))
	PFZ1list.append((Entropy1-Ener1/templist[i]))
	# print "Entropy2=",Entropy2
	# print "Ener2=",Ener2
	# print "Temp=",templist[i]
	# print "pfz=",Entropy2-Ener2/templist[i]
	# if Entropy2<-1.7e+308:
	# 	PFZ2list.append(0.0)
	# else:
	# 	PFZ2list.append(exp(Entropy2-Ener2/templist[i]))
	PFZ2list.append((Entropy2-Ener2/templist[i]))

EME2=(Ener1-Ener2)
EME2=EME2*EME2

CVlist=[]

for i in range(0,templen):
	tempPFZ1=exp(PFZ1list[i]-PFZ2list[i])
	tempPFZ2=exp(PFZ2list[i]-PFZ1list[i])
	CVlist.append(EME2/templist[i]/templist[i]/(tempPFZ1+tempPFZ2+2.0))

try:
	fp3=open(filename[0:len(filename)-4]+'localcv.dat', 'w')
	print " writing information into : [ "+filename[0:len(filename)-4]+'localcv.dat'+"]."
except IOError:
	print " can not open file: [ "+filename[0:len(filename)-4]+'localcv.dat'+"]."
	exit(-1)

for i in range(0,templen):
	print >> fp3, "%f %f" % (templist[i],CVlist[i])

print " output done."
fp3.close()


CVlist=[]
templist=[]
ener=[]
entropy=[]
PFZ1list=[]
PFZ2list=[]
os.system("paste deMCEE.dat "+filename[0:len(filename)-4]+'localcv.dat'+" | awk '{print $1,$2-$4}' > newcv.dat")
#os.system("paste deMCEE.dat deceEnerLJ001001.dat | awk '{print $1,$2-$4}' > minusBFcv.dat")
os.system("gnuplot < draw_newcv.gpl")

exit(0)


