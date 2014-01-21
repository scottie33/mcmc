#!/bin/bash

if [ $# -ne 4 ]; then
	#              $1      $2  $3       $4
	echo " cli: cmd inputfn col outputfn smoothintvl"
	echo " dependency: [smoothdata] "
	exit -1
fi

inputfn=$1
col=$2
outputfn=$3
sintvl=$4
shift
shift
shift
shift

cp $inputfn anotherinput

totalcol=`head -n 1 anotherinput | wc -w`
echo "there are $totalcol columns"
xindex=$(($totalcol+2))
echo "the col-index in pasted file is: $xindex"

echo "cat anotherinput | awk '{print \$1,\$$col}' > tempcol.dat" > temp1.bash
source temp1.bash
rm temp1.bash
./smoothdata tempcol.dat tempoutput $sintvl
rm tempcol.dat
echo -n "paste anotherinput tempoutput | awk '{print " > temp2.bash
for((i=1;i<=$totalcol;i++));do
	if [ $i -ne $col ];then
		if [ $i -ne $totalcol ];then
			echo -n "\$$i," >> temp2.bash
		else
			echo -n "\$$i" >> temp2.bash
		fi
	else
		if [ $i -ne $totalcol ];then
			echo -n "\$$xindex," >> temp2.bash
		else 
			echo -n "\$$xindex" >> temp2.bash
		fi
	fi
done
echo "}' > $outputfn" >> temp2.bash

source temp2.bash
rm temp2.bash

rm tempoutput
rm anotherinput

exit 0