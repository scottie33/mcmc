proc bonddistance {seltext f_r_out f_d_out r_min r_max N_d} {
  
  source aligncom.tcl
  aligncom "$seltext" ;#all frame alighed

  set sel [atomselect top "$seltext"]
  set natoms [$sel num]
  #set sel [$sel get index]
  #puts $sel
  puts " $natoms atoms loaded..."
  set comall [measure center $sel]
  foreach {a b c} $comall { 
    #coding nothing;
  }
  
  set itern [expr $natoms]
  set istor [expr $natoms*$natoms]
  for {set m 0} {$m < $istor} {incr m} {
    #set mindata($m.r) 1000000.0
    set simdata($m.r) 0.0
    set simdata2($m.r) 0.0
    set numdata($m.r) 0.0
    #set maxdata($m.r) 0.0
  }
  set ncuttot 0.0

  set nf [molinfo top get numframes]

  set dr 0.0
  if { $N_d > 1 } { 
    set dr [expr ($r_max - $r_min) /($N_d - 1)]
  }
  for {set k 0} {$k < $N_d} {incr k} {
    set distribution($k) 0
  }
  set maxrange [expr $r_max*$r_max]
  ##################################################
  # Loop over all frames.                          #
  ##################################################
  set outfile [open $f_r_out w]
  set realnf 0.0
  for {set i 0} {$i < $nf} {incr i} {
    set selcut [atomselect top "(((x-$a)*(x-$a)+(y-$b)*(y-$b)+(z-$c)*(z-$c))<=$maxrange)" frame $i]
    set ncuttot [expr $ncuttot+[$selcut num]]
    if { $ncuttot > 1.0 } {
      set realnf [expr $realnf+1]
    } else {
      continue
    }

    puts "frame $i of $nf selection made" 
    set selindex [$selcut get index]
    set selindext $selindex
    #puts $selt
    foreach m $selindex {
      set sel1 [atomselect top "index $m" frame $i] 
      set com1 [measure center $sel1]
      set selindext [lreplace $selindext 0 0]
      foreach n $selindext {
        set sel2 [atomselect top "index $n" frame $i]
        set com2 [measure center $sel2]
        $sel2 delete
        #set com1 [measure center $sel1 weight mass]
        #set com2 [measure center $sel2 weight mass]
        set blen [veclength [vecsub $com1 $com2]]
        if { $m < $n } {
          set tempm [expr $m*$itern+$n]
        } else {
          set tempm [expr $n*$itern+$m]
        }
        # if {$blen<$mindata($tempm.r)} {
        #   set mindata($tempm.r) $blen
        # }
        # if {$blen>$maxdata($tempm.r)} {
        #   set maxdata($tempm.r) $blen
        # }
        #if { $blen < $r_max } {
        set numdata($tempm.r) [expr $numdata($tempm.r)+1.0]
        set simdata($tempm.r) [expr $simdata($tempm.r)+$blen]
        set simdata2($tempm.r) [expr $simdata2($tempm.r)+$blen*$blen]
        if { $N_d > 1 } {
          set k [expr int(($blen - $r_min) / $dr)]
          incr distribution($k)
        }
        #}
      }
      $sel1 delete
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

  set ncuttot [expr $ncuttot/$realnf]
  puts "ncuttot=$ncuttot"
  for {set m 0} {$m < $itern} {incr m} {
    for {set n [expr $m+1]} {$n < $itern} {incr n} {
      set tempm [expr $m*$itern+$n]
      #set simdata($tempm.r) [vecscale $simdata($tempm.r) [expr 1.0/$numdata($tempm.r)]] ;# <r_{i,j}>
      if { $numdata($tempm.r) > 0.0 } {
        set simdata($tempm.r) [vecscale $simdata($tempm.r) [expr 1.0/$numdata($tempm.r)]] ;# <r_{i,j}>
      }
      #set simdata2($tempm.r) [vecscale $simdata2($tempm.r) [expr 1.0/$numdata($tempm.r)]] ;# <r_{i,j}>
      if { $numdata($tempm.r) > 0.0 } {
        set simdata2($tempm.r) [vecscale $simdata2($tempm.r) [expr 1.0/$numdata($tempm.r)]] ;# <r^2_{i,j}>
      }
      #set numdata($tempm.r) [expr $numdata($tempm.r)/$nf] ;# (N_cut-1)*N_cut/2
    }
  }
  puts " all index data averaged over!"
  # if { $N_d > 1 } { 
  #   for {set m 0} {$m < $istor} {incr m} {
  #     set mean [expr $mean+$simdata($m.r)]
  #   }
  #   set mean [expr $mean/$istor] ;# <r>
  # }
  # puts "mean=$mean"
  # if { $N_d > 1 } { 
  #   for {set k 0} {$k < $N_d} {incr k} { 
  #     set tempr [expr $r_min+$k*$dr-$mean] ;# r-<r>  
  #     set tempr [expr $tempr*$tempr]       ;# (r-<r>)^2 
  #     set deldata [expr $deldata+$distribution($k)*$tempr ] 
  #   } 
  #   set deldata [expr $deldata/$sum]
  # }
  # puts "delta=[expr sqrt($deldata)]"
  set outfile2 [open "bondeddis.dat" w]
  puts " creating data files .... "
  # for {set m 0} {$m < $itern} {incr m} {
  #   set numcut($m.r) 0.0
  #   #puts $m
  #   for {set n 0} {$n < $itern} {incr n} {
  #     #puts $n
  #     if { $m != $n } {
  #       if { $n > $m } {
  #         set tempm [expr $m*$itern+$n]
  #       } elseif { $n < $m } {
  #         set tempm [expr $n*$itern+$m]
  #       }
  #       set numcut($m.r) [expr $numcut($m.r)+$numdata($tempm.r)]
  #     }
  #   }
  #   #puts "$m $numcut($m.r)"
  # }
  puts " creating data files .... 2 "
  set factorBERRYallnone 0.0
  set factorBERRYallbond 0.0
  for {set m 0} {$m < $itern} {incr m} {
    #set factorBERRY 0.0
    #set factorBERRYbond 0.0
    for {set n [expr $m+1]} {$n < $itern} {incr n} {
      #if { $n != $m } {
        #if { $n > $m } {
        set tempm [expr $m*$itern+$n]
        #} elseif { $n < $m } {
        #  set tempm [expr $n*$itern+$m]
        #}
        ;# index <r_{i,j}> min(r_{i,j}) max(r_{i,j}) <r^2_{i,j}>
        ;# puts $outfile "[expr $tempm+1] $simdata($tempm.r) $mindata($tempm.r) $maxdata($tempm.r) $mean $deldata"
        ;# index <r_{i,j}> <r^2_{i,j}>
        ;# set factorLINDEMANN [expr $factorLINDEMANN+$simdata2($tempm.r)-$simdata($tempm.r)*$simdata($tempm.r)]
        if { $simdata($tempm.r) > 0.0 } {
          if { [expr $n-$m] == 1 } {
            #puts "$m $n $numdata($tempm.r) $simdata($tempm.r) $simdata2($tempm.r)"
            set tempnums [ expr $simdata2($tempm.r)-$simdata($tempm.r)*$simdata($tempm.r) ]
            if { $tempnums > 0.0 } {
              set factorBERRYallbond [expr $factorBERRYallbond+ sqrt($tempnums) / $simdata($tempm.r)]
            }
            puts $outfile2 "[expr $tempm] $simdata($tempm.r) $simdata2($tempm.r)"
          } else {
            #puts "$m $n $numdata($tempm.r) $simdata($tempm.r) $simdata2($tempm.r)"
            set tempnums [ expr $simdata2($tempm.r)-$simdata($tempm.r)*$simdata($tempm.r) ]
            if { $tempnums > 0.0 } {
              set factorBERRYallnone [expr $factorBERRYallnone+ sqrt($tempnums) / $simdata($tempm.r)]
            }
            puts $outfile "[expr $tempm] $simdata($tempm.r) $simdata2($tempm.r)"
          }
        }
      #}
    }
    #set factorBERRYallnone [expr $factorBERRYallnone+$factorBERRY]
    #set factorBERRYallbond [expr $factorBERRYallbond+$factorBERRYbond]
  }
  puts " files closed.... "
  #set factorLINDEMANN [expr sqrt($factorLINDEMANN*2/$itern/($itern-1))]
  set factorBERRYallnone [expr $factorBERRYallnone*2.0/$ncuttot/($ncuttot-1)]
  set factorBERRYallbond [expr $factorBERRYallbond*2.0/$ncuttot/($ncuttot-1)]
  close $outfile
  close $outfile2
  ##################################################
  puts " creating factor files..."
  set outfile [open "factors.dat" w]
  puts $outfile "factorBERRY=$factorBERRYallnone"
  puts $outfile "factorBERRYbond=$factorBERRYallbond"
  puts $outfile "factorBERRYall=[expr $factorBERRYallnone+$factorBERRYallbond]"
  close $outfile
  puts " factor files closed..."
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


