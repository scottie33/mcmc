
//ssbf.cpp

#include "ssbf.h"
///////////////////////////////////////////////////////
vector<string> Split(const string TEMP) {
	vector<string> RET;
	int i=0;
	int j=0;
	int SizeStr=TEMP.size();
	while( ( i!=SizeStr ) ) {
		while( (TEMP[i]!=' ') && (TEMP[i]!='\t') && (i!=SizeStr) ) {
			i++;
		}
		if(i!=j) {
			RET.push_back( TEMP.substr(j, i-j) );
		}
		while( (TEMP[i]==' ') || (TEMP[i]=='\t') ) {
			i++;
		}
		j=i;
	}
	return RET;
}
///////////////////////////////////////////////////////
string Split_By_Underline(const string TEMP) {
	string RET=string("");
	int i=0;
	//int j=0;
	int SizeStr=TEMP.size();
	while( i!=SizeStr && TEMP.substr(i, 1)!=string("_") ) {
		i++;
	}
	if( i<=SizeStr-2 ) {
		i++;
		RET=TEMP.substr(i, SizeStr-i);
	}
	return RET;
}
///////////////////////////////////////////////////////
string Split_By_Underline_Header(const string TEMP) {
	string RET=string("");
	int i=0;
	//int j=0;
	int SizeStr=TEMP.size();
	while( i!=SizeStr && TEMP.substr(i, 1)!=string("_") ) {
		i++;
	}
	if( i<=SizeStr-1 ) {
		RET=TEMP.substr(0, i);
	}
	return RET;
}
///////////////////////////////////////////////////////
string LBFilter(const string InStr) {
	//string Ret;
	int i=0;
	int Size_Str=InStr.size();
	while( ( (InStr.c_str()[i]==' ') || (InStr.c_str()[i]=='\t') ) && (i!=Size_Str) ) {
		i++;
	}
	//Ret=InStr.substr(i, Size_Str-i);
	//return Ret;
	return InStr.substr(i, Size_Str-i);
}
///////////////////////////////////////////////////////////////////
string RBFilter(const string InStr) {
	//string Ret;
	int Size_Str=InStr.size();
	int i=Size_Str-1;
	while( (i>=0) && 
		   ( (InStr.c_str()[i]==' ') || (InStr.c_str()[i]=='\t') || 
		   	 (InStr.c_str()[i]=='\r') || (InStr.c_str()[i]=='\n') ) ) {
		i--;
	}
	//Ret=InStr.substr(0, i+1);
	//return Ret;
	return InStr.substr(0, i+1);
}
/////////////////////////////////////////////////
string BFilter(const string InStr) {
	return RBFilter( LBFilter(InStr) );
}
///////////////////////////////////////////////////////


