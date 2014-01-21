
#include <iostream>
#include "ssbf.cpp"
#include <fstream>
#include <cstdlib>
using namespace std;

void smoothArray(const double *in, double *out, long N, int n) {
    long k, n1, n2;
    int i;
    for(k=0; k<N; k++) {
        n1=k-n;
        if(n1<0) n1=0;
        n2=k+n;
        if(n2>N-1) n2=N-1; 
        out[k]=0;
        for(i=n1; i<=n2; i++) {
            out[k]+=in[i]/(n2-n1+1);
        }
    }
}

int main(int argc, char** argv) {
    if(argc<4) {
        cout<<" cmd: [PWD]smoothdata inputfn outputfn smoothstep"<<endl;
        exit(-1);
    }

    char* filename;
    filename=argv[1];

    char* outputfn;
    outputfn=argv[2];

    int smoothstep;
    smoothstep=atoi(argv[3]);

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
    for(i=0;i<tempsz;i++) {
        alldata[i]=alldatavec[i];
    }
    cout<<" all data loaded."<<endl;

    smoothArray(alldata, alldata_new, tempsz, smoothstep);

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

    para_ofstream.close();

    cout<<" all data written."<<endl;

    xaxis.clear();
    alldatavec.clear();
    delete[] alldata_new;
    delete[] alldata;
    return 1;
}

