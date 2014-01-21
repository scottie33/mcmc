
#include <iostream>
#include "ssbf.cpp"
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;


int main(int argc, char** argv) {
    if(argc<3) {
        cout<<" cmd: [PWD]integration inputfn outputfn"<<endl;
        cout<<" warning: pay attention to the input data's order!"<<endl;
        exit(-1);
    }

    char* filename;
    filename=argv[1];

    char* outputfn;
    outputfn=argv[2];

    string FILENAME_para;
    FILENAME_para=string(filename);
    ifstream para_ifstream(FILENAME_para.c_str());
    if(para_ifstream==NULL) {
        cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
        exit(-1);
    }
    cout<<" Loading data from file [ "<<FILENAME_para<<" ] ..."<<endl;
    
    vector<double> xaxis;
    vector<double> alldatavec;
    vector<string> tempvec;
    int tempsz=0;
    string tempstr;
    while( getline(para_ifstream, tempstr) ) {
        tempvec=Split(BFilter(tempstr));
        tempsz=tempvec.size();
        if( tempsz<2 ) {
            cout<<" error, tempvec.size()="<<tempvec.size()<<endl;
            exit(-1);
        } else {
            xaxis.push_back(atof(tempvec[0].c_str()));
            alldatavec.push_back(atof(tempvec[1].c_str()));
        }
    }
    para_ifstream.close();
    tempsz=alldatavec.size();
    int i=0;
    double* alldata;
    double* alldata_new;
    alldata=new double[tempsz];
    alldata_new=new double[tempsz];
    for(i=0; i<tempsz; i++) {
        alldata[i]=alldatavec[i];
        alldata_new[i]=0.0;
    }
    cout<<" all data loaded."<<endl;

    //integration procedure;
    double sum=0.0;
    for(i=1; i<tempsz; i++ ) {
       sum+=(alldata[i]+alldata[i-1])/2.0*(xaxis[i]-xaxis[i-1]);
       alldata_new[i]=sum;
    }
    alldata_new[0]==alldata_new[1];
    //integration procedure ends here;

    FILENAME_para=string(outputfn);
    ofstream para_ofstream(FILENAME_para.c_str());
    if(para_ofstream==NULL) {
        cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
        exit(-1);
    }
    cout<<" Writing data into file [ "<<FILENAME_para<<" ] ..."<<endl;

    for(i=0;i<tempsz;i++) {
        para_ofstream<<xaxis[i]<<" "<<alldata_new[i]<<endl; 
    } 
    for(i=1;i<tempsz;i++) {
        if(alldata_new[i]*alldata_new[i-1]<0.0) {
            if( fabs(alldata_new[i])<fabs(alldata_new[i-1]) ) {
                cout<<i<<" "<<xaxis[i]<<" "<<alldata[i]<<" "<<alldata_new[i]<<endl; 
            } else {
                cout<<i<<" "<<xaxis[i-1]<<" "<<alldata[i-1]<<" "<<alldata_new[i-1]<<endl; 
            }
            i+=1;
        }
    }

    para_ofstream.close();

    cout<<" all data written."<<endl;

    xaxis.clear();
    alldatavec.clear();
    delete[] alldata_new;
    delete[] alldata;
    return 1;
}

