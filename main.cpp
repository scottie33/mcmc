#include "cmc.h"
#include <ctime>


int main(int argc, char **argv) {
  cmc mymcsystem;
  //mymcsystem.init_conformation_a(65,1,4.7,"pnc.pdb");
  double start=clock();
  mymcsystem.init_mpi(argc, argv);
  mymcsystem.initialization();
  mymcsystem.run();
  mymcsystem.ensembler();
  mymcsystem.output_statistic();
  mymcsystem.term_mpi();
  double end=clock();
  cout<<" [time]="<<(end-start)/CLOCKS_PER_SEC<<"seconds. "<<endl;//16~18s
  return 0;
}

//end of file;