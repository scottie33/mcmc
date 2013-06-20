#include "cmc.h"
#include <ctime>


int main(int argc, char **argv) {
  /*CMyArray3<double> tarray;
  tarray.Build(3,2,100);
  tarray.SetZero();
  int i, j, k;
  double temp=1.0;
  for(i=0; i<3; i++) {
    for(j=0; j<2; j++) {
      for(k=0; k<100; k++) {
        cout<<" i="<<i<<" j="<<j<<" k="<<k<<endl;
        tarray.pArray[i][j][k]=temp;
        temp+=1.0;
      }
    }
  }
  for(i=0; i<3; i++) {
    for(j=0; j<2; j++) {
      for(k=0; k<100; k++) {
        cout<<" "<<tarray.pArray[i][j][k]<<endl;
      }
    }
  }*/
  cmc mymcsystem;
  double start=clock();
  mymcsystem.init_mpi(argc, argv);
  mymcsystem.initialization();
  mymcsystem.run();
  mymcsystem.ensembler();
  mymcsystem.output_statistic();
  mymcsystem.term_mpi();
  double end=clock();
  cout<<" [time]="<<(end-start)/CLOCKS_PER_SEC<<"seconds. "<<endl;//16~18s/**/
  return 0;
}

//end of file;