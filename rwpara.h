#ifndef _RWPARA_H 
#define _RWPARA_H

#include "ssbf.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

template<typename T>
void writeparameter(ofstream &OFSName, const string ParaName, T Para) {
	OFSName<<setw(20)<<setiosflags(ios::left)<<ParaName.c_str()<<" "
		   <<setw(20)<<Para<<resetiosflags(ios::left)<<endl;
}
////////////////////////////////////
void readparameter(const string TempString, const string ParaName, bool &BPara);
void readparameter(const string TempString, const string ParaName, int &IPara);
void readparameter(const string TempString, const string ParaName, float &FPara);
void readparameter(const string TempString, const string ParaName, double &DPara);
void readparameter(const string TempString, const string ParaName, string &SPara);

/*void writeparameter(ofstream &OFSName, const string ParaName, bool BPrara);
void writeparameter(ofstream &OFSName, const string ParaName, int IPrara);
void writeparameter(ofstream &OFSName, const string ParaName, float FPrara);
void writeparameter(ofstream &OFSName, const string ParaName, double DPrara);
void writeparameter(ofstream &OFSName, const string ParaName, string SPrara);
*/

//////////
#endif