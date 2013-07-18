#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 4: 
	print " Syntax: getMCEFE.py lowtemp hightemp nTbins" # [shift]"
	print " e.g.: getFE.py 0.01 10.0 100" # [shift]"
  	exit(-1)

Tlow=float(sys.argv[1])
Thigh=float(sys.argv[2])
NREM=int(sys.argv[3])

try:
	tfp=open("paraentropy.gpl", 'r')
	print " loading information of column from: [paraentropy.gpl]"
except IOError:
	print " can not open file: [paraentropy.gpl]"
	exit(-1)

epsilon=1.0
while True:
	line=tfp.readline().rstrip()
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
tfp.close()

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
try:
	hfp=open("MCEFE.dat", 'w')
	print " writing free energy information into: [MCEFE.dat]"
except IOError:
	print " can not open file: [MCEFE.dat]"
	exit(-1)

templist=[]
for i in range(0,NREM):
	templist.append(Thigh/pow(Thigh/Tlow, float(i)/float(NREM-1)))


for j in range(0,listlen):
	for i in range(0,NREM):
	#PFZ=-1.7e+308
	#meanEn=-1.7e+308
	#meanEp=-1.7e+308
	
		##print entropy[i]-ener[i]/valTemp
		#PFZ=log_cc(PFZ, entropy[i]-ener[i]/valTemp)
		#if ener[i]<0.0: 
		#	meanEn=log_cc(meanEn, log(-ener[i])+entropy[i]-ener[i]/valTemp)
		#elif ener[i]>0.0:
		#	meanEp=log_cc(meanEp, log(ener[i])+entropy[i]-ener[i]/valTemp)
	#meanEp=exp(meanEp-PFZ)-exp(meanEn-PFZ)
		if entropy[j]<-1.0e308:
			#print >> hfp, "%f %f %f" % (ener[i], valTemp, entropy[j])
			print >> hfp, "%f %f %f" % (ener[j], templist[i], entropy[j])
		else:
			#print >> hfp, "%f %f %f" % (ener[i], valTemp,ener[j]-valTemp*entropy[j])
			print >> hfp, "%f %f %f" % (ener[j], templist[i], ener[j]-templist[i]*entropy[j])

	print >> hfp, " "

hfp.close()
print " please check your [ MCEFE.dat ]."
print " format: E1 T1 fe11"
print "         E1 T2 fe12"
print "         E ... ..."
print "         EN TM feNM"

#os.system("echo \"set xrange ["+str(ener[0])+":"+str(ener[listlen-1])+"]\" > range.gpl")
os.system("echo \"set yrange ["+str(Tlow/epsilon)+":"+str(Thigh/epsilon)+"]\" >> range.gpl")
os.system("gnuplot < draw_mcefe.gpl")
print " please check [ mcefe.eps ]."
os.system("convert -rotate 90 mcefe.eps mcefe.png")
print " please check [ mcefe.png ]."

exit(0)



