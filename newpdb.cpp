
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
