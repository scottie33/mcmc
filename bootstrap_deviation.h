
#include <cstdlib>
#include <algorithm>
#include "autocorrelation.h"

int bootstrap_deviation(double* arrayp, int arrayn, int numB, double* ret);
int bootstrap_deviation(vector<double> arrayvec, int arrayn, int numB, double* ret);

int bootstrap_deviation(double* arrayp, int arrayn, int numB, double* ret) { // better choose numB>1000;
  /*if(numB<1000) {
  	return NULL;
  }*/
  //double ret[2];//0 mean value, 1 deviation;
  //cout<<arrayp[1999]<<endl;
  //getchar();
  double* resample;
  resample=new double[numB]; // sum of resumple suffix i
  int i;
  //cout<<arrayp[1999]<<endl;
  for(i=0; i<numB; i++) {
  	resample[i]=0.0;
  }
  //cout<<arrayp[1999]<<endl;
  int j;

  //cout<<" realarray[last]="<<realarray[--i]<<endl;
  //CMyArray<double> newarray;
  //newarray.Build(numB,arrayn);
  //newarray.SetZero();

  srand(numB*arrayn);
  for(i=0; i<numB; i++) {
  	for(j=0; j<arrayn; j++) {
  		/*int temp=rand()%arrayn;
  		//resample[i]+=arrayp[temp];
  		cout<<temp<<" "<<arrayp[temp]<<endl;*/
  		//newarray.pArray[i][j]=arrayp[rand()%arrayn];
  		//resample[i]+=newarray.pArray[i][j];
      resample[i]+=arrayp[rand()%arrayn];
  	}
  	resample[i]/=double(arrayn);
  	//cout<<resample[i]<<endl;
  	//getchar();
  }
  double meanvalue=0.0;
  double deviation=0.0;
  for(i=0; i<numB; i++) {
  	meanvalue+=resample[i];
  	deviation+=resample[i]*resample[i];
  }
  delete[] resample;
  ret[0]=meanvalue/double(numB);
  ret[1]=deviation/double(numB)-ret[0]*ret[0];
  //cout<<newarray.pArray[3][150];
  //delete[] realarray;
  //newarray.Release();
  return 0;
}

int bootstrap_deviation(vector<double> arrayvec, int arrayn, int numB, double* ret) {
	double* sample;
	sample=new double[arrayn];
	for(int i=0; i<arrayn; i++) {
  		sample[i]=arrayvec[i];
  		//cout<<arrayvec[i]<<endl;
  		//getchar();
  	}
  	//double ret[2];
  	bootstrap_deviation(sample, arrayn, numB, ret);
  	delete[] sample;
  	return 0;
}