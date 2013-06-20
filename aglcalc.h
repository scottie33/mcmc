
#ifndef _AGLCALC_H
#define _AGLCALC_H

#include <cmath>
#include <iostream>
// gcc over 4.3 requires explicitly including the following headers.
#include <cstring> 
#include <cstdlib>
#include <cstring>
#include <algorithm>
using namespace std;

#define _MIN_DOUBLE     -1.7e308
#define _MAX_DOUBLE      1.7e308
#define PRECISION        1e-6
#define IOERROR         -1
#define LOGICERROR      -2
#define SIZEERROR       -3
#define NODEFINEERROR   -4
#define MEMOSETZERO(POINTER, SIZEOFDATA) memset(POINTER, 0, SIZEOFDATA)
#define MEMOCOPY(TPOINTER, RPOINTER, SIZEOFDATA) memcpy(TPOINTER, RPOINTER, SIZEOFDATA)
template<class T> 
inline T MEMOSIZE(T SizeNUM) { 
	return SizeNUM+1; 
}
/////////////////////////////////////////////////////////////////////
#define SQUARE(X) X*X
#define TRIPLE(X) X*X*X
//template<class T> inline T QUADRI(T x) { return (x)*(x)*(x)*(x); }
//template<class T> inline T QUIQUE(T x) { return (x)*(x)*(x)*(x)*(x); }
#define HEXIAL(X) X*X*X*X*X*X
//template<class T> inline T MaxNum(T A, T B) { return A>B?A:B; }
//template<class T> inline T MinNum(T A, T B) { return A<B?A:B; }
/////////////////////////////////////
template<class T> class CMyArray {
private:
	bool released;
	int Dim1_Size, Dim2_Size;
public:
	T** pArray;
	CMyArray(void);
	~CMyArray(void);
	void Build(int Size1, int Size2);
	void SetZero();
	void Release(void);
};
/////////////////////////////////////
//////////////////////////
template<class T> CMyArray<T>::CMyArray(void) {
	pArray=NULL;
	released=true;
	Dim1_Size=0;
	Dim2_Size=0;
}
//////////////////////////
template<class T> CMyArray<T>::~CMyArray(void) {
	if(!released) {
		this->Release();
	}
}
//////////////////////////
template<class T> void CMyArray<T>::Build(int Size1, int Size2) {
	pArray=NULL;
	Dim1_Size=Size1;
	Dim2_Size=Size2;
	int i=0;
	pArray=new T*[Dim1_Size];
	i=Dim1_Size*Dim2_Size;
	pArray[0]=new T[i];
	MEMOSETZERO(pArray[0], sizeof(T)*i);
	for(i=1; i<Dim1_Size; i++) {
		pArray[i]=pArray[i-1]+Dim2_Size;
	}
	released=false;
}
template<class T> void CMyArray<T>::SetZero() {
	int totalsize=Dim1_Size*Dim2_Size;
	if( totalsize==0 || released ) {
		cout<<" can not setzero: CMyArray<T>::SetZero"<<endl;
		exit(LOGICERROR);
	}
	MEMOSETZERO(pArray[0], sizeof(T)*totalsize);
}
template<class T> void CMyArray<T>::Release(void) {
	if(pArray) {
		if(pArray[0]) {
			delete[] pArray[0];
		}
		delete[] pArray;
		pArray=NULL;
	}
	Dim1_Size=Dim2_Size=0;
	released=true;
}
///////////////////////
/////////////////////////////////////
template<class T> class CMyArray3 {
private:
	bool released;
	int Dim1_Size, Dim2_Size, Dim3_Size;
public:
	T*** pArray;
	CMyArray3(void);
	~CMyArray3(void);
	void Build(int Size1, int Size2, int Size3);
	void SetZero();
	void Release(void);
};
/////////////////////////////////////
//////////////////////////
template<class T> CMyArray3<T>::CMyArray3(void) {
	pArray=NULL;
	released=true;
	Dim1_Size=0;
	Dim2_Size=0;
	Dim3_Size=0;
}
//////////////////////////
template<class T> CMyArray3<T>::~CMyArray3(void) {
	if(!released) {
		this->Release();
	}
}
//////////////////////////
template<class T> void CMyArray3<T>::Build(int Size1, int Size2, int Size3) {
	pArray=NULL;
	Dim1_Size=Size1;
	Dim2_Size=Size2;
	Dim3_Size=Size3;
	int i=0;
	int j=0;
	pArray=new T**[Dim1_Size]; // x
	for(i=0; i<Dim1_Size; i++) {
		pArray[i]=new T*[Dim2_Size]; // y
	}
	i=Dim1_Size*Dim2_Size*Dim3_Size;
	pArray[0][0]=new T[i];
	MEMOSETZERO(pArray[0][0], sizeof(T)*i);
	for(j=1; j<Dim2_Size; j++) {
		pArray[0][j]=pArray[0][j-1]+Dim3_Size;
	}
	for(i=1; i<Dim1_Size; i++) {
		for(j=0; j<Dim2_Size; j++) {
			pArray[i][j]=pArray[0][0]+Dim3_Size*(i*Dim2_Size+j);
		}
	}
	released=false;
}
//////////////////////////
//////////////////////////
template<class T> void CMyArray3<T>::SetZero() {
	int totalsize=Dim1_Size*Dim2_Size*Dim3_Size;
	if( totalsize==0 || released ) {
		cout<<" can not setzero: CMyArray<T>::SetZero"<<endl;
		exit(LOGICERROR);
	}
	MEMOSETZERO(pArray[0][0], sizeof(T)*totalsize);
}
//////////////////////////
template<class T> void CMyArray3<T>::Release(void) {
	if(pArray) {
		if(pArray[0]) {
			for(int i=Dim1_Size; i>0; i--) {
				delete[] pArray[i-1];
			}
		}
		delete[] pArray;
		pArray=NULL;
	}
	Dim1_Size=Dim2_Size=Dim3_Size=0;
	released=true;
}
///////////////////////
//////////////////////////////////////////////
template<class T>
inline T log_cc(const T &log_aa, const T &log_bb) { //cc=aa+bb; give log_aa, log_bb, return log_cc;
	if(log_aa>log_bb) {
		return log_aa+log(1+exp(log_bb-log_aa));
	} else {
		return log_bb+log(1+exp(log_aa-log_bb));
	}
}
/////////////////////////////////////////////////////////
/*template<class T>
inline T distance_2(T x1, T y1, T z1, 
					T x2, T y2, T z2) {
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}
////////////////////////
template<class T>
inline void T_display(T** temp_T) {
	cout<<"  "<<temp_T[0][0]<<" "<<temp_T[0][1]<<" "<<temp_T[0][2]<<endl;
	cout<<"  "<<temp_T[1][0]<<" "<<temp_T[1][1]<<" "<<temp_T[1][2]<<endl;
	cout<<"  "<<temp_T[2][0]<<" "<<temp_T[2][1]<<" "<<temp_T[2][2]<<endl;
}
//////////////////////////////////////////////
template<class T>
inline void P_display(T* temp_P) {
	cout<<"  "<<temp_P[0]<<endl;
	cout<<"  "<<temp_P[1]<<endl;
	cout<<"  "<<temp_P[2]<<endl;
}
//////////////////////////////////////////////
template<class T>
inline T Omega_Det(T sin_num, T cos_num) {
	if(sin_num>=0.0) {
		return acos(cos_num);
	} else {
		return (-1.0)*acos(cos_num);
	}
}
//////////////////////////////////////////////
template<class T>
inline T Test_Value_of_Angle(T somevalue) {
	string flag="";
	if( somevalue>1.00 ) {
		cout<<" this value is over limit: "<<somevalue<<endl;
		cout<<" change to 1.00 ?(Yes/No)"<<endl;
		cin>>flag;
		cout<<endl;
		if( flag=="Yes" || flag=="Y" ) {
			return 1.00;
		} else {
			cout<<" there's some problem in ur calculation, check them and try again! (Test_Value_of_Angle)"<<endl;
			exit(LOGICERROR);
		}
	}
	if( somevalue<-1.00 ) {
		cout<<" this value is over limit:"<<somevalue<<endl;
		cout<<" change to -1.00 ?(Yes/No)"<<endl;
		cin>>flag;
		cout<<endl;
		if( flag=="Yes" || flag=="Y" ) {
			return -1.00;
		} else {
			cout<<" there's some problem in ur calculation, check them and try again! (Test_Value_of_Angle)"<<endl;
			exit(LOGICERROR);
		}
	}
	return somevalue;
}
//////////////////////////////////////////////
template<class T>
inline T** T_x(T someangle) {
	T** Tret;
	Tret=new T*[3];
	for(int i=0; i!=3; i++) {
		Tret[i]=new T[3];
	}
	Tret[0][0]=1.0;
	Tret[0][1]=Tret[0][2]=Tret[1][0]=Tret[2][0]=0.0;
	Tret[1][1]=Tret[2][2]=cos(someangle);
	Tret[1][2]=-sin(someangle);
	Tret[2][1]=-Tret[1][2];
	return Tret;
}
///////////////////////////////////////
template<class T>
inline T** T_y(T someangle) {
	T** Tret;
	Tret=new T*[3];
	for(int i=0; i!=3; i++) {
		Tret[i]=new T[3];
	}
	Tret[0][0]=Tret[2][2]=cos(someangle);
	Tret[0][2]=sin(someangle);
	Tret[2][0]=-Tret[0][2];
	Tret[0][1]=Tret[1][2]=Tret[1][0]=Tret[2][1]=0.0;
	Tret[1][1]=1.0;
	return Tret;
}
///////////////////////////////////////
template<class T>
inline T** T_z(T someangle) {
	T** Tret;
	Tret=new T*[3];
	for(int i=0; i!=3; i++) {
		Tret[i]=new T[3];
	}
	Tret[0][0]=Tret[1][1]=cos(someangle);
	Tret[0][1]=-sin(someangle);
	Tret[1][0]=-Tret[0][1];
	Tret[0][2]=Tret[1][2]=Tret[2][0]=Tret[2][1]=0.0;
	Tret[2][2]=1.0;
	return Tret;
}
///////////////////////////////////////
template<class T>
inline T* Add_p(T* p1, T *p2) {
	T* retp;
	retp=new T[3];
	retp[0]=p1[0]+p2[0];
	retp[1]=p1[1]+p2[1];
	retp[2]=p1[2]+p2[2];
	return retp;
}
///////////////////////////////////////
template<class T>
inline T* Substract_p(T* p1, T *p2) {
	T* retp;
	retp=new T[3];
	retp[0]=p1[0]-p2[0];
	retp[1]=p1[1]-p2[1];
	retp[2]=p1[2]-p2[2];
	return retp;
}
//////////////////////////////////////////////
template<class T>
inline T* C_p(T TConst, T* pi) {
	T* qret;
	qret=new T[3];
	for(int i=0; i!=3; i++)
	{
		qret[i]=(T_angle[i][0]*pi[0]+T_angle[i][1]*pi[1]+T_angle[i][2]*pi[2]);
	}
	qret[0]=TConst*pi[0];
	qret[1]=TConst*pi[1];
	qret[2]=TConst*pi[2];
	return qret;
}
//////////////////////////////////////////////
template<class T>
inline T* C_p_r(T TConst, T* pi) {
	T* qret;
	qret=new T[3];
	for(int i=0; i!=3; i++)
	{
		qret[i]=(T_angle[i][0]*pi[0]+T_angle[i][1]*pi[1]+T_angle[i][2]*pi[2]);
	}
	qret[0]=pi[0]/TConst;
	qret[1]=pi[1]/TConst;
	qret[2]=pi[2]/TConst;
	return qret;
}
//////////////////////////////////////////////
template<class T>
inline T* T_p(T** T_angle, T* pi) {
	T* qret;
	qret=new T[3];
	for(int i=0; i!=3; i++)
	{
		qret[i]=(T_angle[i][0]*pi[0]+T_angle[i][1]*pi[1]+T_angle[i][2]*pi[2]);
	}
	qret[0]=(T_angle[0][0]*pi[0]+T_angle[0][1]*pi[1]+T_angle[0][2]*pi[2]);
	qret[1]=(T_angle[1][0]*pi[0]+T_angle[1][1]*pi[1]+T_angle[1][2]*pi[2]);
	qret[2]=(T_angle[2][0]*pi[0]+T_angle[2][1]*pi[1]+T_angle[2][2]*pi[2]);
	return qret;
}
//////////////////////////////////////////////
template<class T>
inline T* p_X_p(T* p1, T* p2) {
	T* qret;
	qret=new T[3];
	qret[0]=(p1[1]*p2[2]-p1[2]*p2[1]);
	qret[1]=(p1[2]*p2[0]-p1[0]*p2[2]);
	qret[2]=(p1[0]*p2[1]-p1[1]*p2[0]);
	return qret;
}
//////////////////////////////////////////////////////
template<class T>
inline T cos_angle(T x1, T y1, T z1, 
				T x2, T y2, T z2,
				T x3, T y3, T z3) {
	T x21=x2-x1;	T y21=y2-y1;	T z21=z2-z1;
	T x32=x3-x2;	T y32=y3-y2;	T z32=z3-z2;
	T r21=sqrt( x21*x21+y21*y21+z21*z21 );
	T r32=sqrt( x32*x32+y32*y32+z32*z32 );
	//following R21.R32/r21/r32 
	T ANGLE_TEMP=( x21*x32+y21*y32+z21*z32 )/r21/r32;
	//cout<<DIHEDRAL_TEMP<<endl;
	//cout<<acos(DIHEDRAL_TEMP)<<endl;
	///////////////////////important!//////////////////////
	if(ANGLE_TEMP>1.0) {
		return 1.0;
	} else if(ANGLE_TEMP<-1.0) {
		return -1.0;
	} else {
		return ANGLE_TEMP;
	}
}
//////////////////////////////////////////////////////
template<class T>
inline T cos_angle(T x1, T y1, T z1, 
			T x2, T y2, T z2,
			T x3, T y3, T z3,
			T x4, T y4, T z4) {
	T x21=x2-x1;	T y21=y2-y1;	T z21=z2-z1;
	T x43=x4-x3;	T y43=y4-y3;	T z43=z4-z3;
	T r21=sqrt( x21*x21+y21*y21+z21*z21 );
	T r43=sqrt( x43*x43+y43*y43+z43*z43 );
	//following R21.R32/r21/r32 
	T ANGLE_TEMP=( x21*x43+y21*y43+z21*z43 )/r21/r43;
	//cout<<DIHEDRAL_TEMP<<endl;
	//cout<<acos(DIHEDRAL_TEMP)<<endl;
	///////////////////////important!//////////////////////
	if(ANGLE_TEMP>1.0) {
		return 1.0;
	} else if(ANGLE_TEMP<-1.0) {
		return -1.0;
	} else {
		return ANGLE_TEMP;
	}
}
//////////////////////////////////////////////////////
template<class T>
inline T angle(T x1, T y1, T z1, 
			T x2, T y2, T z2,
			T x3, T y3, T z3) {
	T x21=x2-x1;	T y21=y2-y1;	T z21=z2-z1;
	T x32=x3-x2;	T y32=y3-y2;	T z32=z3-z2;
	T r21=sqrt( x21*x21+y21*y21+z21*z21 );
	T r32=sqrt( x32*x32+y32*y32+z32*z32 );
	//following R21.R32/r21/r32 
	T ANGLE_TEMP=( x21*x32+y21*y32+z21*z32 )/r21/r32;
	//cout<<DIHEDRAL_TEMP<<endl;
	//cout<<acos(DIHEDRAL_TEMP)<<endl;
	///////////////////////important!//////////////////////
	if(ANGLE_TEMP>1.0) {
		ANGLE_TEMP=1.0;
	} else if(ANGLE_TEMP<-1.0) {
		ANGLE_TEMP=-1.0;
	}
	return acos( ANGLE_TEMP );//acos(theta)
}*/
//////////////////////////////////////////////////////
template<class T>
inline T dihedral(T x1, T y1, T z1, 
				T x2, T y2, T z2,
				T x3, T y3, T z3,
				T x4, T y4, T z4) {	//see ref. hw00_dihedral.pdf
	T x21=x2-x1;	T y21=y2-y1;	T z21=z2-z1;
	T x32=x3-x2;	T y32=y3-y2;	T z32=z3-z2;
	T x43=x4-x3;	T y43=y4-y3;	T z43=z4-z3;
	T N1_x=y21*z32-y32*z21;//N1=R21xR32;
	T N1_y=z21*x32-z32*x21;
	T N1_z=x21*y32-x32*y21;
	T N2_x=y32*z43-y43*z32;//N2=R32xR43;
	T N2_y=z32*x43-z43*x32;
	T N2_z=x32*y43-x43*y32;	
	/*T n1=sqrt( N1_x*N1_x+N1_y*N1_y+N1_z*N1_z );
	T n2=sqrt( N2_x*N2_x+N2_y*N2_y+N2_z*N2_z );
	if( fabs(n1*n2)<1e-12 ) {
		if( fabs(n1)<1e-12 ) {
			cout<<" Error: Normal_1 can not be 0.000 !"<<endl;
			exit(LOGICERROR);
		} else {
			cout<<" Error: Normal_2 can not be 0.000 !"<<endl;
			exit(LOGICERROR);
		}
	}*/
	//T DIHEDRAL_cos=( N1_x*N2_x+N1_y*N2_y+N1_z*N2_z )/n1/n2;
	T DIHEDRAL_cos=N1_x*N2_x+N1_y*N2_y+N1_z*N2_z;
	//cout<<DIHEDRAL_TEMP<<endl;
	//cout<<acos(DIHEDRAL_TEMP)<<endl;
	///////////////////////important!//////////////////////
	/*if(DIHEDRAL_cos>1.0) {
		DIHEDRAL_cos=1.0;
	} else if(DIHEDRAL_cos<-1.0) {
		DIHEDRAL_cos=-1.0;
	}*/
	///////////////////////////////////////////////////////
	/*T N3_x=N1_y*N2_z-N1_z*N2_y; //y21*z32-y32*z21;//N3=N1xN2;
	T N3_y=N1_z*N2_x-N1_x*N2_z; //z21*x32-z32*x21;
	T N3_z=N1_x*N2_y-N1_y*N2_x; //x21*y32-x32*y21;
	T DIHEDRAL_sin=(sqrt(N3_x*N3_x+N3_y*N3_y+N3_z*N3_z))/n1/n2;*/
	T DIHEDRAL_sin=sqrt(x32*x32+y32*y32+z32*z32)*(x43*N1_x+y43*N1_y+z43*N1_z);
	//cout<<DIHEDRAL_TEMP<<endl;
	////////// rotate n1 to map on n2, if ..... ////////////////////////////
	/*T N1_y_again=N1_y*cos(DIHEDRAL_TEMP)-N1_z*sin(DIHEDRAL_TEMP);
	T N1_z_again=N1_y*sin(DIHEDRAL_TEMP)+N1_z*cos(DIHEDRAL_TEMP);
	n1=sqrt( N1_x*N1_x+N1_y_again*N1_y_again+N1_z_again*N1_z_again );
	T DELTA=( N1_x*N2_x+N1_y_again*N2_y+N1_z_again*N2_z )/n1/n2;
	if( fabs(DELTA-1.0)<1e-12 )
	{
		return DIHEDRAL_TEMP;
	}
	return -DIHEDRAL_TEMP;*/
	////////// much easier .... ////////////////////////////
	/*T Angle=x21*N2_x+y21*N2_y+z21*N2_z;//cosine of angle between R21 and N2;
	------------------------
	   Angle<0.0, N1->N2 + 
	   Angle=0.0, N1==N2
	   Angle>0.0, N1->N2 - 
	--------------------------
	if( Angle<=0.0 )
	{
		return -DIHEDRAL_TEMP;
	}*/
	return atan2(DIHEDRAL_sin,DIHEDRAL_cos);
}
////////////////////////////

template<class T>
inline T dihedral2(T x21, T y21, T z21, 
			       T x32, T y32, T z32,
			       T x43, T y43, T z43)	{	//see ref. hw00_dihedral.pdf
	T N1_x=y21*z32-y32*z21;//N1=R21xR32;
	T N1_y=z21*x32-z32*x21;
	T N1_z=x21*y32-x32*y21;
	T N2_x=y32*z43-y43*z32;//N2=R32xR43;
	T N2_y=z32*x43-z43*x32;
	T N2_z=x32*y43-x43*y32;	
	/*T n1=sqrt( N1_x*N1_x+N1_y*N1_y+N1_z*N1_z );
	T n2=sqrt( N2_x*N2_x+N2_y*N2_y+N2_z*N2_z );
	if( fabs(n1*n2)<1e-12 ) {
		if( fabs(n1)<1e-12 ) {
			cout<<" Error: Normal_1 can not be 0.000 !"<<endl;
			exit(LOGICERROR);
		} else {
			cout<<" Error: Normal_2 can not be 0.000 !"<<endl;
			exit(LOGICERROR);
		}
	}*/
	//T DIHEDRAL_cos=( N1_x*N2_x+N1_y*N2_y+N1_z*N2_z )/n1/n2;
	T DIHEDRAL_cos=N1_x*N2_x+N1_y*N2_y+N1_z*N2_z;
	//cout<<DIHEDRAL_TEMP<<endl;
	//cout<<acos(DIHEDRAL_TEMP)<<endl;
	///////////////////////important!//////////////////////
	/*if(DIHEDRAL_cos>1.0) {
		DIHEDRAL_cos=1.0;
	} else if(DIHEDRAL_cos<-1.0) {
		DIHEDRAL_cos=-1.0;
	}*/
	///////////////////////////////////////////////////////
	/*T N3_x=N1_y*N2_z-N1_z*N2_y; //y21*z32-y32*z21;//N3=N1xN2;
	T N3_y=N1_z*N2_x-N1_x*N2_z; //z21*x32-z32*x21;
	T N3_z=N1_x*N2_y-N1_y*N2_x; //x21*y32-x32*y21;
	T DIHEDRAL_sin=(sqrt(N3_x*N3_x+N3_y*N3_y+N3_z*N3_z))/n1/n2;*/
	T DIHEDRAL_sin=sqrt(x32*x32+y32*y32+z32*z32)*(x43*N1_x+y43*N1_y+z43*N1_z);
	//cout<<DIHEDRAL_TEMP<<endl;
	////////// rotate n1 to map on n2, if ..... ////////////////////////////
	/*T N1_y_again=N1_y*cos(DIHEDRAL_TEMP)-N1_z*sin(DIHEDRAL_TEMP);
	T N1_z_again=N1_y*sin(DIHEDRAL_TEMP)+N1_z*cos(DIHEDRAL_TEMP);
	n1=sqrt( N1_x*N1_x+N1_y_again*N1_y_again+N1_z_again*N1_z_again );
	T DELTA=( N1_x*N2_x+N1_y_again*N2_y+N1_z_again*N2_z )/n1/n2;
	if( fabs(DELTA-1.0)<1e-12 )
	{
		return DIHEDRAL_TEMP;
	}
	return -DIHEDRAL_TEMP;*/
	////////// much easier .... ////////////////////////////
	/*T Angle=x21*N2_x+y21*N2_y+z21*N2_z;//cosine of angle between R21 and N2;
	------------------------
	   Angle<0.0, N1->N2 + 
	   Angle=0.0, N1==N2
	   Angle>0.0, N1->N2 - 
	--------------------------
	if( Angle<=0.0 )
	{
		return -DIHEDRAL_TEMP;
	}*/
	return atan2(DIHEDRAL_sin,DIHEDRAL_cos);
}
////////////////////////////

#endif

