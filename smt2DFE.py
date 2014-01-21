#!/usr/bin/python

import sys,os
from math import exp,log,pow,sqrt


if len(sys.argv) < 2: 
	print " Syntax: smt2DFE.py ifilename" # [shift]"
  	exit(-1)

filename=sys.argv[1]

try:
	oefp=open(filename, 'r')
	print " loading information of column from: [ "+filename+" ]"
except IOError:
	print " can not open file: [",filename,"]"
	exit(-1)

tfname=filename+".smoothed.dat"
try:
	outfp=open(tfname, 'w')
	print " writing information into: [ "+tfname+" ]"
except IOError:
	print " can not open file: [ "+tfname+" ]"
	exit(-1)



i=0
x=0
leftval1=0.0
leftval2=0.0
while True:
	line=oefp.readline()
	if line:
		line=line.rstrip()
		elements=line.split()
		if len(elements)==3:
			tempnum=float(elements[2])
			if x==0 :
				print >>outfp, "%s %s %s" % ( elements[0], elements[1], tempnum )
			elif x==1 :
				print >>outfp, "%s %s %f" % ( elements[0], elements[1], (tempnum+leftval2)/2.0 )
			else:
				print >>outfp, "%s %s %f" % ( elements[0], elements[1], (tempnum+leftval2+leftval1)/3.0 )
			leftval1=leftval2 #leftval1 leftval2 tempnum
			leftval2=tempnum #update
			x=x+1
		if len(elements)==0:
			#print i
			print >>outfp, ""
			#i=0
			x=0
			i=i+1
	else:
		break
print "",i,"data series loaded." 

oefp.close()
outfp.close()

exit(0)


