#!/bin/bash


#725 750 1263 1511 2950 6776

#./getrdf.bash pnc.psf ./pdbfiles 6776 0.01 5.0 0 1 0.5 # in rdf.tcl change sel1 and sel2 to index 27 and not index 27 respectively;
#./getrdf.bash pnc.psf ./pdbfiles 2950 0.01 5.0 0 1 0.6 &> doall2950.log &
#./getrdf.bash pnc.psf ./pdbfiles 1511 0.01 5.0 0 1 0.6 &> doall1511.log &
#./getrdf.bash pnc.psf ./pdbfiles 1263 0.01 5.0 0 1 0.5 &> doall1263.log &
#./getrdf.bash pnc.psf ./pdbfiles 0750 0.01 5.0 0 1 0.5 
#./getrdf.bash pnc.psf ./pdbfiles 0725 0.025 1.0 0 1 0.25
#./showdisfluc.bash pnc.psf ./pdbfiles 0725 0.00 1.50 350 0 10000
./showdisfluc.bash pnc.psf ./pdbfiles 0750 0.00 1.0 350 0 10000
#./showdisfluc.bash pnc.psf ./pdbfiles 1263 0.00 1.50 350 0 50000
#./showdisfluc.bash pnc.psf ./pdbfiles 1511 0.00 1.50 350 0 50000
#./showdisfluc.bash pnc.psf ./pdbfiles 2950 0.00 1.50 350 0 50000 
#./showdisfluc.bash pnc.psf ./pdbfiles 6776 0.00 1.50 350 0 50000


exit