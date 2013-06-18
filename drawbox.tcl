

set celltx 125.0
set cellty 125.0
set celltz 125.0

set cellox 0.0
set celloy 0.0
set celloz 0.0

set cellx [expr $celltx-$cellox]
set celly [expr $cellty-$celloy]
set cellz [expr $celltz-$celloz]

set centx [expr ($celltx+$cellox)/2.0]
set centy [expr ($cellty+$celloy)/2.0]
set centz [expr ($celltz+$celloz)/2.0]

puts "{$cellx $celly $cellz}"
set lenvec {}
lappend lenvec $cellx
lappend lenvec $celly
lappend lenvec $cellz
pbc set $lenvec -all -molid top

set sel [atomselect top "resid 2"]
pbc box -centersel $sel