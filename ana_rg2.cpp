

#include "ssbf.h"
#include "aglcalc.h"
#include "bootstrap_deviation.h"
#include <cstdlib>
#include <algorithm>
#include <fstream>
using namespace std;

int main(int argc, char** argv) {
  if(argc<4) {
    cout<<" cmd filelist bsgap bsnum [numofcol]"<<endl;
    cout<<" please specify the 'filelist' and bsgap (d'better thinking before use)."<<endl; 
    cout<<" warning: if 'numofcol' not specified, all data in a row will be added."<<endl;
    return -1;
  }
  int bsgap=atoi(argv[2]);
  int bsnum=atoi(argv[3]);
  cout<<" bsgap="<<bsgap<<" bsnum="<<bsnum<<endl;
  int numofcol=0;
  if(argc>=5) {
    numofcol=atoi(argv[4]);
  }
  int i=0; int j=0;
  int tempsize=0;
  string tempstr;
  vector<string> tempvec;
  vector<double> data_x; 
  double data_x_sum=0.0;
  double data_x2_sum=0.0;
  int temp_sz;
  double tempnum;
  string filename;
  double dev[2];
  double adg[3];
  //string rm_filelist=string("rm -fr ")+string(argv[1]);
  string fileprefix=string(argv[1]).substr(0,3);
  /*string createfilelist=string("if [ `ls ")+fileprefix //move to script
                       +string("*.dat | wc -l` -gt 0 ]; then ls ")
                       +fileprefix+string("*.dat > ")
                       +fileprefix+string("_file.lst; fi");*/
  /////////////////////////////////////
  //system(rm_filelist.c_str());
  //string filelistname=fileprefix+string("_file.lst");
  ifstream rg2filelist( argv[1] );
  if(rg2filelist==NULL) {
    cout<<" no [ "<<string(argv[1])<<" ], exit... " <<endl;
    rg2filelist.close();
    return -1;
  }
  //rg2filelist.open( filelistname.c_str() );
  cout<<endl<<" Loading filename from file [ "<<string(argv[1])<<" ] ..."<<endl;
  string outputdatafile=fileprefix+string("_data.lst");
  ofstream rg2datalist(outputdatafile.c_str());
  if(rg2datalist==NULL) {
    cout<<" Error: Can not open file: [ "<<outputdatafile<<" ]"<<endl;
    exit(IOERROR);
  }
  rg2datalist<<"#N1M m2an lo3 4igh de5^2 6mean bde7^2 bsg8p bsn9m g10 bd11v^2*g n_eff12"<<endl;
  cout<<endl<<" writing filename into file [ "<<outputdatafile<<" ] ..."<<endl;
  while( getline(rg2filelist, tempstr) ) {
    tempvec=Split(BFilter(tempstr)); // a file
    if(tempvec.size()!=0) {
      filename=tempvec[0];
      ifstream rg2file(filename.c_str());
      if(rg2file==NULL) {
        cout<<" Error: Can not open file: [ "<<filename<<" ]"<<endl;
        exit(IOERROR);
      }
      cout<<" loading data from:[ "<<filename<<" ] ..."<<endl;
      data_x_sum=0.0;
      data_x2_sum=0.0;
      data_x.clear();
      while( getline(rg2file, tempstr) ) {
        tempvec=Split(BFilter(tempstr));
        tempsize=tempvec.size();
        if(tempsize!=0) {
          tempnum=0.0;
          if(numofcol>tempsize) {
            cout<<" column > datacol... error. "<<endl;
            return -1;
          } else if(numofcol==0){
            for(i=0;i<tempsize;i++) tempnum+=atof(tempvec[i].c_str());
          } else {
            for(i=0;i<numofcol;i++) tempnum+=atof(tempvec[i].c_str());
          }
          data_x.push_back(tempnum);
          //data_x_sum+=tempnum;
          //data_x2_sum+=tempnum*tempnum;
          //tempnum=sqrt(tempnum);
        }
      }//all data loaded;
      rg2file.close();
      temp_sz=data_x.size();
      //tempnum=data_x2_sum/double(temp_sz)-( data_x_sum/double(temp_sz) )*( data_x_sum/double(temp_sz) );//dev^2
      cout<<" "<<temp_sz<<" data loaded... "<<endl;
      int realarraynum=int(temp_sz/bsgap);
      cout<<" arraylen="<<realarraynum<<endl;
      double* realarray;
      realarray=new double[realarraynum];
      for(i=0, j=0; j<temp_sz; j+=bsgap) {
        realarray[i++]=data_x[j];
      }
      bootstrap_deviation(realarray,realarraynum,bsnum,dev);
      autocorrelation(realarray,realarraynum,adg);
      delete[] realarray;;
      if(filename.size()>9) {
        rg2datalist<<atoi(filename.substr(9,5).c_str())<<" "
                   <<adg[0]<<" "
                   <<*min_element(data_x.begin(),data_x.end())<<" "
                   <<*max_element(data_x.begin(),data_x.end())<<" "
                   <<adg[1]<<" "
                   <<dev[0]<<" "<<dev[1]<<" "
                   <<bsgap<<" "<<bsnum<<" "<<adg[2]<<" "
                   <<dev[1]*adg[2]/double(realarraynum)<<" "<<double(realarraynum)/adg[2]<<endl;
      } else {
        rg2datalist<<(filename)<<" " //1
                   <<adg[0]<<" "     //2
                   <<*min_element(data_x.begin(),data_x.end())<<" "
                   <<*max_element(data_x.begin(),data_x.end())<<" "
                   <<adg[1]<<" " //mean
                   <<dev[0]<<" "<<dev[1]<<" " //6mean //7bdev
                   <<bsgap<<" "<<bsnum<<" "<<adg[2]<<" "
                   <<adg[1]*adg[2]/double(realarraynum)<<" "<<double(realarraynum)/adg[2]<<endl;
      }
      //mean, low, high, deviation, bs_mean, bs_dev;
      //cout<<" passed."<<endl;
    }
  } 
  rg2datalist.close();
  rg2filelist.close();
  //system( (string("echo \"fn=\'")+fileprefix+string("_data\'\" > tempfn.gpl")).c_str() );   
  //system( (string("echo \"yname=\'")+fileprefix+string(" [0.1nm]\'\" >> tempfn.gpl")).c_str() );  
  //system("gnuplot < rg2list.gpl"); //move to script
  //system("./creat_html.bash"); //move to script
  return 0;
}
