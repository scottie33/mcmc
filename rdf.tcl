

proc getrdf {seltext1 seltext2 deltar radiimax PBCorNot index} {

	set outputfile1 "grdf_${index}.dat"

	set sel1 [atomselect top "$seltext1"]

	set sel2 [atomselect top "$seltext2"]

	puts " now RDF calculation start, please wait ..."

	set outfp1 [open $outputfile1 w]

	set gr0 [measure gofr $sel1 $sel2 delta $deltar rmax $radiimax usepbc $PBCorNot selupdate 1 first 0 last -1 step 1]

	set r [lindex $gr0 0]
	set gr [lindex $gr0 1]
	set igr [lindex $gr0 2]
	set isto [lindex $gr0 3]
	foreach j $r k $gr l $igr m $isto {
	  puts $outfp1 [format "%.4f\t%.4f\t%.4f\t%.4f" $j $k $l $m]
	}

	close $outfp1

	puts " you have your RDF \[$outputfile1\] now, we then shall do the contacting number parts."
}



proc getrdfcom {seltext1 seltext2 deltar radiimax PBCorNot index comrlower} {

	set outputfile1 "grdf_${index}.dat"
	source aligncom.tcl
	aligncom "$seltext1"
	set sel1 [atomselect top "$seltext1"]
	set com1 [measure center $sel1]
	foreach {a b c} $com1 { 
		#coding nothing;
	}

	set comrlower [expr $comrlower*$comrlower]

	set sel1 [atomselect top "(((x-$a)*(x-$a)+(y-$b)*(y-$b)+(z-$c)*(z-$c))<=$comrlower)"]
	set sel2 [atomselect top "(((x-$a)*(x-$a)+(y-$b)*(y-$b)+(z-$c)*(z-$c))>$comrlower)"]

	puts " now RDF calculation start, please wait ..."
	set outfp1 [open $outputfile1 w]

	set gr0 [measure gofr $sel1 $sel2 delta $deltar rmax $radiimax usepbc $PBCorNot selupdate 1 first 0 last -1 step 1]

	set r [lindex $gr0 0]
	set gr [lindex $gr0 1]
	set igr [lindex $gr0 2]
	set isto [lindex $gr0 3]
	foreach j $r k $gr l $igr m $isto {
	  puts $outfp1 [format "%.4f\t%.4f\t%.4f\t%.4f" $j $k $l $m]
	}

	close $outfp1

	puts " you have your RDF \[$outputfile1\] now, we then shall do the contacting number parts."
}



