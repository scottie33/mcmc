#!/usr/bin/python

import sys,os
from math import exp,log

if len(sys.argv) < 4: 
	print " Syntax: getcvfromprob.py filename colid valTemp" # [shift]"
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
	oefp=open(filename, 'r')
	print " loading information of column from: [",filename,"]"
except IOError:
	print " can not open file: [",filename,"]"
	exit(-1)

nume1=0.0
nume2=0.0
numer=0.0
tmpn=0
while True:
	line=oefp.readline().rstrip()
	if line:
		elements=line.split()
		if len(elements)>=colid:
			ener=float(elements[0])
			enum=float(elements[colid-1])
			numer=numer+enum
			nume1=nume1+enum*ener
			nume2=nume2+enum*ener*ener
			tmpn=tmpn+1
	else:
		break
oefp.close()
print " there are",tmpn,"lines information loaded."
print " len_row=",len(elements),"colnum=",colid
tempcv=0.0
tempe=0.0
if numer!=0:
	tempe=nume1/numer
	tempcv=((nume2/numer-nume1*nume1/numer/numer)/valTemp/valTemp)
else:
	tempcv=0

tempstr="echo "+str(valTemp)+" "+str(tempe)+" "+str(tempcv)+" > OET.dat"

print " T=",valTemp,"E=",tempe,"CECV=",tempcv
os.system(tempstr)

exit(0)


