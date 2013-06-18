#!/usr/bin/python

import os,sys
import math

if len(sys.argv) < 5 :
	print " Syntax: calcdis.py filename halflenx halfleny halflenz" # [shift]"
  	exit(-1)

try:
	fp=open(sys.argv[1], 'r')
	print " loading information of column from:",sys.argv[1]
except IOError:
	print " can not open file:",sys.argv[1]
	exit(-1)

halflenx=float(sys.argv[2])
halfleny=float(sys.argv[3])
halflenz=float(sys.argv[4])

enelist=[]
dislist=[]
while True:
	line=fp.readline().rstrip()
	if line:
		elements=line.split()
		if len(elements)==8:
			distXCoor=abs( float(elements[1])-float(elements[5]) )
			if distXCoor>halflenx:
				distXCoor=(halflenx*2.0-distXCoor)
			distYCoor=abs( float(elements[2])-float(elements[6]) )
			if distYCoor>halfleny:
				distYCoor=(halfleny*2.0-distYCoor)
			distZCoor=abs( float(elements[3])-float(elements[7]) )
			if distZCoor>halflenz:
				distZCoor=(halflenz*2.0-distZCoor)
			disttemp=math.sqrt(distXCoor*distXCoor+distYCoor*distYCoor+distZCoor*distZCoor)
			enelist.append(float(elements[0]))
			dislist.append(disttemp)
	else:
		break
fp.close()

try:
	fpw=open(sys.argv[1]+".dat", 'w')
	print " writing information of column into: "+sys.argv[1]+".dat"
except IOError:
	print " can not open file: "+sys.argv[1]+".dat"
	exit(-1)

lenlist=len(enelist)

for i in range(0,lenlist):
	print >> fpw, "%10.3f %10.3f" % (enelist[i], dislist[i])

fpw.close()




