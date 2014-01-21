source aligncom.tcl

proc bonddistance {seltext f_r_out f_d_out r_min r_max N_d} {
  set sel [atomselect top "$seltext"]
  set natoms [$sel num]
  set sel [$sel get index]
  puts $sel
  puts " $natoms atoms loaded..."
  
  set itern [expr $natoms-1]
  for {set m 0} {$m < $itern} {incr m} {
    set mindata($m.r) 1000000.0
    set simdata($m.r) 0.0
    set maxdata($m.r) 0.0
  }
  set nf [molinfo top get numframes]

  set dr 0.0
  if { $N_d > 1 } { 
    set dr [expr ($r_max - $r_min) /($N_d - 1)]
  }
  for {set k 0} {$k < $N_d} {incr k} {
    set distribution($k) 0
  }
  ##################################################
  # Loop over all frames.                          #
  ##################################################
  set outfile [open $f_r_out w]
  for {set i 0} {$i < $nf} {incr i} {
    puts "frame $i of $nf" 
    for {set m 0} {$m < $itern} {incr m} {
      set sel1 [atomselect top "index $m" frame $i]
      set sel2 [atomselect top "index [expr $m+1]" frame $i]
      set com1 [measure center $sel1]
      set com2 [measure center $sel2]
      #set com1 [measure center $sel1 weight mass]
      #set com2 [measure center $sel2 weight mass]
      set blen [veclength [vecsub $com1 $com2]]
      if {$blen<$mindata($m.r)} {
        set mindata($m.r) $blen
      }
      if {$blen>$maxdata($m.r)} {
        set maxdata($m.r) $blen
      }
      set simdata($m.r) [expr $simdata($m.r)+$blen]
      if { $N_d > 1 } {
        set k [expr int(($blen - $r_min) / $dr)]
        incr distribution($k)
      }
    }
  }
  set sum 0.0
  set mean 0.0
  set deldata 0.0
  if { $N_d > 1 } { 
    for {set k 0} {$k < $N_d} {incr k} {
      set sum [expr $sum+$distribution($k)]
    }
  }
  for {set m 0} {$m < $itern} {incr m} {
    set simdata($m.r) [vecscale $simdata($m.r) [expr 1.0/$nf]]
  }

  if { $N_d > 1 } { 
    for {set m 0} {$m < $itern} {incr m} {
      set mean [expr $mean+$simdata($m.r)]
    }
    set mean [expr $mean/$itern]
  }
  puts "mean=$mean"
  if { $N_d > 1 } { 
    for {set k 0} {$k < $N_d} {incr k} {
      set tempr [expr $r_min+$k*$dr-$mean]
      set tempr [expr $tempr*$tempr]
      set deldata [expr $deldata+$distribution($k)*$tempr ]
    }
    set deldata [expr $deldata/$sum]
  }
  puts "delta=[expr sqrt($deldata)]"
  for {set m 0} {$m < $itern} {incr m} {
    puts $outfile "[expr $m+1] $simdata($m.r) $mindata($m.r) $maxdata($m.r) $mean $deldata"
  }
  close $outfile
  ##################################################
  ##################################################

  if { $N_d > 1 } { 
    set outfile [open $f_d_out w]
    for {set k 0} {$k < $N_d} {incr k} {
      puts $outfile "[expr $r_min + $k * $dr] [expr $distribution($k)/$sum]"
    }
    close $outfile
  } else {
    puts " N of bins is $N_d, no distribution can be got, quit."
  }
}


