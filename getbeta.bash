#!/bin/bash

if [ $# -lt 6 ]; then
	echo " cmd: xmin xmax y1min y1max y2min y2max Natoms"
	exit -1
fi 


./smoothdata beta.dat smbeta.dat 20
./derivative smbeta.dat desmbeta.dat 
./smoothdata desmbeta.dat smdesmbeta.dat 20

cat smdesmbeta.dat | awk '{if($2<0) {print $1,-$2} else {print $1,$2}}' > absbeta.dat

echo "set xrange [$1:$2]" > newrangebeta.gpl
echo "set yrange [$3:$4]" >> newrangebeta.gpl
echo "set y2range [$5:$6]" >> newrangebeta.gpl
echo "set arrow 1 from $1, second 0.0 to $2, second 0.0 nohead lt 4 lw 1 lc rgb 'black'" >> newrangebeta.gpl

#cat smdesmbeta.dat | awk '{if($1<$) {print $1,-$2} else {print $1,$2}}' > absbeta.dat
#tempE=`cat deMCEE.dat | sort -k2,2n | tail -n 1 | awk '{print $1}'`
#echo "set arrow 2 from $tempE, second $5 to $tempE, second $6 nohead lt 2 lw 1 lc rgb 'black'" >> newrangebeta.gpl
gnuplot < draw_smdebeta.gpl

exit 1