proc aligncom {seltext {outputornot 0}} {
	set outDataFile [open "rmsd.dat" w]
	#set sela [atomselect top all]
	set selb [atomselect top all]
	set alignsela [atomselect top "$seltext"]
	set alignselb [atomselect top "$seltext"]
	set nf [molinfo top get numframes]
	$alignsela frame 0
	if { $outputornot > 0 } {
		puts $outDataFile "0 0.0"
	}
	for {set fb 1} {$fb<$nf} {incr fb} {
		puts " aligning $fb-th frame ... "
		$selb frame $fb
		$alignselb frame $fb
		#display update
		#set val 
		#set resid $r
		set trans_mat [measure fit $alignselb $alignsela]
		$selb move $trans_mat
		if { $outputornot > 0 } {
			puts $outDataFile "$fb [measure rmsd $sela $selb]"
		} ;#puts $outDataFile "$fa $fb [measure rmsd $sela $selb weight mass]"
	}
	#puts $outDataFile " "
	close $outDataFile
}