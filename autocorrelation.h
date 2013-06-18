
#include <cstdlib>
#include <algorithm>
#include "aglcalc.h"
int smooth(double* arrayp, int arrayn);
int autocorrelation(double* arrayp, int numN, double* ret);
int autocorrelation(vector<double> arrayp, int numN, double* ret);

int smooth(double* arrayp, int arrayn) {
  int i=0;
  //int maxi=0;
  //double maxgap=0.0;
  //int j=1;
  bool flag=true;
  while(flag) {
    flag=false;
    for(i=1;i<arrayn-1;i++) {
      if(arrayp[i]<0.0) {
        arrayp[i]=(arrayp[i-1]+arrayp[i+1])/2.1;
        flag=true;
      }
    }
  }
  cout<<" pass negtive test. but it may be trapped ..."<<endl
      <<" if you wait too long, please press <Ctrl+C> ..."<<endl;
  i=1;
  flag=true;
  while(true) {
    for(;i<arrayn-1;i++) {
      if( (arrayp[i]+arrayp[i-1]) < (arrayp[i]+arrayp[i+1]) ) { //pairsum not decreasing...
        break;
      }
    }
    if( i==(arrayn-1) ) {//all pairs sum decreasing?... 
      if(flag) {//when overlooped, 
        break; 
      } else {
        //flag=true; 
        i=1;
      }
    } else { //smooth;
      /*if(flag) { 
        cout<<" not all sum of pairs are decreasing... break..."<<endl; 
        flag=false;
      }
      cout<<" i="<<i<<" array_i-1="<<arrayp[i-1];
      cout<<" array_i="<<arrayp[i];
      cout<<" array_i+1="<<arrayp[i+1]<<endl;
      getchar();*/
      arrayp[i]=(arrayp[i-1]+arrayp[i+1])/2.1;
      if(i==arrayn-2) {
        i=1;
        //flag=true; //overloop;
      } else { 
        //arrayp[i+1]=(arrayp[i+2]+arrayp[i])/2.0;
        i=i+1; 
        //flag=false; //not overloop;
      }
      //break;
    }
  }
  return 0;
}

int autocorrelation(double* arrayp, int numN, double* ret) { //matrix=MxN;
  
  //double data_sum=0.0;
  //double data2_sum=0.0;
  double data_avg=0.0;
  double data2_avg=0.0;
  double dev2=0.0;
  double g_const=0.0;
  int i,k;

  bool testinfo=false;

  for(i=0; i<numN; i++) {
  	data_avg+=arrayp[i];
    data2_avg+=arrayp[i]*arrayp[i];
  }
  data_avg=data_avg/double(numN);
  ret[0]=data_avg;
  data2_avg=data2_avg/double(numN);
  dev2=data2_avg-data_avg*data_avg;
  ret[1]=dev2;

  double* A_k; // the normalized autocorrelation function;
  A_k=new double[numN];
  for(k=0;k<numN;k++) {
    A_k[k]=0.0;
  }

  A_k[0]=1.0;
  int end_k=0;
  double tau=0.0;
  for(k=1;k<numN;k++) {
    for(i=0;i<numN-k;i++) {
      A_k[k]+=(arrayp[i]-data_avg)*(arrayp[i+k]-data_avg);
    }
    A_k[k]=A_k[k]/double(numN-k); // co-variance;
    A_k[k]=A_k[k]/dev2; // co-variance/variance;
    A_k[k]=A_k[k]*( 1-double(k)/double(numN) ); // co-eff.
    if ( (A_k[k]+A_k[k-1])<=0.0 ) { // && (k>20) ) 
      end_k=k;
      A_k[end_k]=0.0;
      if(testinfo) cout<<A_k[k]<<endl;
      break;
    } else {
      if(testinfo) cout<<A_k[k]<<" ";
    }
  }

  smooth(A_k,end_k+1); 
  //need to be improved ... from A_k[0]=1.0 to A_k[end_k] smoothed and convex all the time;

  for(k=1;k<end_k;k++) { //since A_k[end_k]=0.0;
    if(testinfo) cout<<A_k[k]<<" ";
    tau=tau+A_k[k];  
  }
  if(testinfo) cout<<A_k[k]<<endl;
  g_const=1.0+2.0*tau;
  delete[] A_k;
  if(g_const<1.0) {
    ret[2]=1.0; 
  } else { 
    ret[2]=g_const;
  }
  return 0;
}
int autocorrelation(vector<double> arrayvec, int numN, double* ret) {
  double* sample;
  sample=new double[numN];
  for(int i=0; i<numN; i++) {
      sample[i]=arrayvec[i];
      //cout<<arrayvec[i]<<endl;
      //getchar();
    }
    //double ret[2];
    autocorrelation(sample, numN, ret);
    delete[] sample;
    return 0;
}
