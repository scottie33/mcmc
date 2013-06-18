
#include "cmc.h"
#include <ctime>

int main(int argc, char **argv) {
  if(argc<4) {
    cout<<" cmd seglen chainlen nchain"<<endl;
    exit(-1);
  }
  double seglen=atof(argv[1]);
  double chainlen=atof(argv[2]);
  double nchain=atof(argv[3]);
  cmc mymcsystem;
  mymcsystem.init_conformation_a(chainlen,nchain,seglen,"pnc.pdb");
  return 0;
}
/*#define SQUARE(X) X*X
template<class T> inline T Energy_LJ1(T EPSILON, T SIGMA, T d_2) {//pay attention d_2=d*d; 
  return 4.0*(fabs(EPSILON)*HEXIAL(SQUARE(SIGMA)/d_2)-EPSILON*TRIPLE(SQUARE(SIGMA)/d_2)); 
}*/
int main(int argc, char **argv) {
  for(int i=4; i<512; i*=2) {
    cmc mymcsystem;
    mymcsystem.init_mpi(argc, argv);
    //mymcsystem;
    //mymcsystem.load("test.pdb",1.0,false);
    //tmol.writepdbinfo("result.pdb",false);
    mymcsystem.init_conformation_a(i,1,4.7,"test.pdb");
    mymcsystem.initialization();
    double start=clock();
    mymcsystem.run_with_stepnumber(i*100000);
    mymcsystem.rand_update(0);
    mymcsystem.run();
    double end=clock();
    cout<<" [time]="<<(end-start)/CLOCKS_PER_SEC<<"seconds. "<<endl;//16~18s
    mymcsystem.fout_conformation(false);
    mymcsystem.memo_free();
  }
  /*cmc mymcsystem;
  //mymcsystem.init_conformation_a(21,1,1.0,"test.pdb");
  mymcsystem.initialization();
  double start=clock();
  mymcsystem.run();
  double end=clock();
  cout<<" [time]="<<(end-start)/CLOCKS_PER_SEC<<"seconds. "<<endl;//16~18s
  */
  return 0;
}


int main(int argc, char **argv) {
  cmc mymcsystem;
  //mymcsystem.init_conformation_a(21,1,4.7,"pnc.pdb");
  mymcsystem.init_mpi(argc, argv);
  mymcsystem.initialization();
  mymcsystem.run();
  mymcsystem.ensembler();
  mymcsystem.term_mpi();
  return 0;
}