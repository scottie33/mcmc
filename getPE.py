#!/usr/bin/python

import sys,os
from math import exp,log

def log_cc(log_aa, log_bb): #cc=aa+bb; give log_aa, log_bb, return log_cc;
	if log_aa>log_bb: 
		return log_aa+log(1+exp(log_bb-log_aa))
	else:
		return log_bb+log(1+exp(log_aa-log_bb))

if len(sys.argv) < 4: 
	print " Syntax: getPE.py valTemp epsilon Natoms [emin] [emax] [idxCE]" # [shift]"
	print " in _temperaturelist.pls T[idxCE]=valTemp, this option is for comparing CE v.s. MCE"
  	exit(-1)

valTemp=float(sys.argv[1])

epsilon=float(sys.argv[2])

natoms=float(sys.argv[3])

emin=0
if len(sys.argv) > 4: 
	emin=float(sys.argv[4])

emax=0
if len(sys.argv) > 5: 
	emax=float(sys.argv[5])

idxCE=0
if len(sys.argv) > 6: 
	idxCE=int(sys.argv[6])
print " emin=",emin
print " emax=",emax

emin=emin*epsilon*natoms
emax=emax*epsilon*natoms

print " emin=",emin
print " emax=",emax
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
			pemax2=pemax2+petemp
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

delmin=1e+308
delmax=1e+308
idemin=0
idemax=0



for i in range(0,listlen):
	if abs(emin-ener[i])<delmin:
		delmin=abs(emin-ener[i])
		idemin=i
	if abs(emax-ener[i])<delmax:
		delmax=abs(emax-ener[i])
		idemax=i
	if len(pelist2)==0:
		print >> pefp, "%f %f %d %f" % (ener[i], pelist1[i], 0, ener[i]-valTemp*(entropy[i]))
		#print >> pefp, "%f %f %d %f" % (ener[i], pelist1[i], 0, ener[i]-valTemp*(entropy[i]-PFZ))
	else:
		print >> pefp, "%f %f %f %f" % (ener[i], pelist1[i], pelist2[i]/pemax2, ener[i]-valTemp*(entropy[i]))
		#print >> pefp, "%f %f %f %f" % (ener[i], pelist1[i], pelist2[i]/pemax2, ener[i]-valTemp*(entropy[i]-PFZ))

ymin=ener[idemin]-valTemp*(entropy[idemin])
#ymin=ener[idemin]-valTemp*(entropy[idemin]-PFZ)
ymax=ymin

for i in range(0,idemax-idemin):
	
	temp=ener[i+idemin]-valTemp*entropy[i+idemin]
	if ymin>temp:
		ymin=temp
	if ymax<temp:
		ymax=temp
		
print " emin=",emin,"i_min=",idemin,"e_i_min=",ener[idemin]
print " emax=",emax,"i_max=",idemax,"e_i_max=",ener[idemax]

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
if abs(emin-0)>1e-12 and abs(emax-0)>1e-12:
	tempstr="echo emin="+str(ener[idemin])+" >> parape.gpl"
	os.system(tempstr)
	tempstr="echo emax="+str(ener[idemax])+" >> parape.gpl"
	os.system(tempstr)

# tempstr="echo ymin="+str(ener[idemin]-valTemp*entropy[idemin])+" >> parape.gpl"
# os.system(tempstr)
# tempstr="echo ymax="+str(ener[idemax]-valTemp*entropy[idemax])+" >> parape.gpl"
# os.system(tempstr)
tempstr="echo ymin="+str(ymin)+" >> parape.gpl"
os.system(tempstr)
#tempstr="echo ymax="+str(ener[idemax]-valTemp*entropy[idemax])+" >> parape.gpl"
tempstr="echo ymax="+str(ymax)+" >> parape.gpl"
os.system(tempstr)
tempstr="echo idxCE="+str(idxCE)+" >> parape.gpl"
os.system(tempstr)
os.system("./smoothdata PE.dat smPE.dat 10")
os.system("gnuplot < draw_PE.gpl")
tempstr="mv PE.eps PE_"+str(valTemp)+".eps"
os.system(tempstr)

print " check [PE_%f.eps]" % (valTemp)

