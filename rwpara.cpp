#include "rwpara.h"


/////////////////////////////////////////////////////
void readparameter(const string TempString, const string ParaName, bool &BPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		if( TempVector[1]==string("true") || TempVector[1]==string("TRUE") ) {
			BPara=true;
        } else if( TempVector[1]==string("false") || TempVector[1]==string("FALSE") ) {
			BPara=false;
        } else {
			cout<<"Error, bool variabes can only be true|TRUE|false|FALSE ... "<<endl;
			exit(-1);
        }
	}
}
void readparameter(const string TempString, const string ParaName, int &IPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		IPara=atoi(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, float &FPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		FPara=atof(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, double &DPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		DPara=atof(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, string &SPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		SPara=TempVector[1];
	}
}





