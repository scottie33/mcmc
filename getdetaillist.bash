#!/bin/bash

if [ $# -lt 1 ]; then
	echo " cmd index1 [index2 [index3 [index4 [index5]]]] ... "
	exit -1
fi

totnum=$#
echo " there are $totnum index need to be done."


cp draw_bondfluclist.gpl draw_bondfluclistt.gpl
echo "set yrange [0.5:3.0]" >> draw_bondfluclistt.gpl
echo "set xtics 5" >> draw_bondfluclistt.gpl
echo "set mxtics 5" >> draw_bondfluclistt.gpl
echo -n "plot " >> draw_bondfluclistt.gpl
num=1

rm tempindex.dat
touch tempindex.dat
for i in $@; do
	echo $i >> tempindex.dat
done

#for i in `cat tempindex.dat | sort -r`; do

for i in $@; do
	if [ $num -eq 1 ]; then
		cid='#708090'
		iid='Ground State'
	elif [ $num -eq 2 ]; then
		cid='#191970'
		iid='S State'
	elif [ $num -eq 3 ]; then
		cid='#FF4500'
		iid='L State'
	elif [ $num -eq 4 ]; then
		cid='#C71585'
		iid='{/Symbol Q}_{min}'
	elif [ $num -eq 5 ]; then
		cid='#006400'
		iid='{/Symbol Q}_{max}'
	else 
		cid='#708090'
		iid='abnormal'
	fi
	if [ $num -lt $totnum ]; then
		echo "'dis_${i}.dat' u 1:2 w l lt 1 lw 5 lc rgb '$cid' title '',\\" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2 w l lt 1 lw 1 lc rgb '$cid' title '{/Times-Italic $i}',\\" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2:3:4 with errorbars lt 1 lw 3 lc rgb '$cid' title '',\\" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2:3:4 with errorbars lt 1 lw 1 lc rgb '$cid' title '{/Times-Italic $iid}',\\" >> draw_bondfluclistt.gpl
	else
		echo "'dis_${i}.dat' u 1:2 w l lt 1 lw 5 lc rgb '$cid' title ''" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2 w l lt 1 lw 1 lc rgb '$cid' title '{/Times-Italic $i}',\\" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2:3:4 with errorbars lt 1 lw 3 lc rgb '$cid' title ''" >> draw_bondfluclistt.gpl
		#echo "'dis_${i}.dat' u 1:2:3:4 with errorbars lt 1 lw 1 lc rgb '$cid' title '{/Times-Italic $iid}'" >> draw_bondfluclistt.gpl
	fi
	let num=$num+1
done
echo "unset multiplot" >> draw_bondfluclistt.gpl
gnuplot < draw_bondfluclistt.gpl
echo " now you can check out your [ bondfluclist.eps ]."



cp draw_bonddislist.gpl draw_bonddislistt.gpl
echo "set xrange [0.8:1.5]" >> draw_bonddislistt.gpl
echo "set ytics 0.05" >> draw_bonddislistt.gpl
echo "set mytics 5" >> draw_bonddislistt.gpl
echo -n "plot " >> draw_bonddislistt.gpl
num=1
for i in $@; do
	if [ $num -eq 1 ]; then
		cid='#708090'
		iid='Ground State'
	elif [ $num -eq 2 ]; then
		cid='#191970'
		iid='S State'
	elif [ $num -eq 3 ]; then
		cid='#FF4500'
		iid='L State'
	elif [ $num -eq 4 ]; then
		cid='#C71585'
		iid='{/Symbol Q}_{min}'
	elif [ $num -eq 5 ]; then
		cid='#006400'
		iid='{/Symbol Q}_{max}'
	else 
		cid='#708090'
		iid='abnormal'
	fi
	if [ $num -lt $totnum ]; then
		echo "'disdis_${i}.dat' u 1:2 w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}',\\" >> draw_bonddislistt.gpl
	else
		echo "'disdis_${i}.dat' u 1:2 w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}'" >> draw_bonddislistt.gpl
	fi
	let num=$num+1
done
echo "unset multiplot" >> draw_bonddislistt.gpl
gnuplot < draw_bonddislistt.gpl
echo " now you can check out your [ bonddistlist.eps ]."



cp draw_rdflist.gpl draw_rdflistt.gpl
echo "set arrow 1 from 0,1.0 to 4,1.0 nohead lt 1 lw 4 lc rgb rgb '#708090'" >> draw_rdflistt.gpl
echo "set xtics 0.5" >> draw_rdflistt.gpl
echo "set mxtics 5" >> draw_rdflistt.gpl
echo "set xrange [0:4]" >> draw_rdflistt.gpl
echo -n "plot " >> draw_rdflistt.gpl
num=1
for i in $@; do
	if [ $num -eq 1 ]; then
		cid='#708090'
		iid='Ground State'
	elif [ $num -eq 2 ]; then
		cid='#191970'
		iid='S State'
	elif [ $num -eq 3 ]; then
		cid='#FF4500'
		iid='L State'
	elif [ $num -eq 4 ]; then
		cid='#C71585'
		iid='{/Symbol Q}_{min}'
	elif [ $num -eq 5 ]; then
		cid='#006400'
		iid='{/Symbol Q}_{max}'
	else 
		cid='#708090'
		iid='abnormal'
	fi
	if [ $num -eq 1 ]; then
		echo "'grdf_${i}.dat' u 1:(\$4*0.06) w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}',\\" >> draw_rdflistt.gpl
	elif [ $num -eq 2 ]; then
		echo "'grdf_${i}.dat' u 1:(\$4) w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}',\\" >> draw_rdflistt.gpl
	elif [ $num -lt $totnum ]; then
		echo "'grdf_${i}.dat' u 1:(\$4) smooth unique w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}',\\" >> draw_rdflistt.gpl
	else
		echo "'grdf_${i}.dat' u 1:(\$4*0.01) w l lt 1 lw 4 lc rgb '$cid' title '{/Times-Italic $iid}'" >> draw_rdflistt.gpl
	fi
	let num=$num+1
done
echo "unset multiplot" >> draw_rdflistt.gpl
gnuplot < draw_rdflistt.gpl
echo " now you can check out your [ rdflist.eps ]."

exit 0