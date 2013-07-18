
#include "cmc.h"
#include <ctime>

int main(int argc, char **argv) {
  if(argc<7) {
    cout<<" cmd seglen chainlen nchain xoff yoff zoff"<<endl;
    exit(-1);
  }
  double seglen=atof(argv[1]);
  double chainlen=atof(argv[2]);
  double nchain=atof(argv[3]);
  double xoffset=atof(argv[4]);
  double yoffset=atof(argv[5]);
  double zoffset=atof(argv[6]);
  cmc mymcsystem;
  mymcsystem.init_conformation_a(chainlen,nchain,seglen,"pnc.pdb",xoffset,yoffset,zoffset);
  return 0;
}
