//energy.h

#ifndef _ENERGY_H
#define _ENERGY_H


#define _MIN_DOUBLE     -1.7e308
#define _MAX_DOUBLE      1.7e308
#define _EMIN_DOUBLE     -1.7e5
#define _EMAX_DOUBLE      1.7e5
#define KBT         1.00
#define e		   -1.00
#define E_x			1.00
#define E_y			0.00
#define E_z			0.00
//static double tempdis24;
static double tempdis12;
static double tempdis6;
static double tempdis2;
//#define Energy_BF(paraK,dMd2) paraK*dMd2/2.0
inline double  Energy_BF(const double&  paraK, const double& dMd, const double& bondlen, const double& bondlen2, 
	const double&  EPSILON, const double&  SIGMAWELL2, const double&  SIGMA6, const double&  lambda, const double&  SIGMA) {
	/*tempdis6=dMd-bondlen;
	return paraK*tempdis6*tempdis6/2.0;*/
	/*tempdis2=dMd/bondlen2;
	tempdis12=1.0-tempdis2;
	tempdis12=tempdis12<1e-3?1e-3:tempdis12;
	if(dMd<SIGMAWELL2) {
		tempdis6=dMd*dMd*dMd;
		tempdis6=SIGMA6/tempdis6;
		return -paraK/2.0*bondlen2*log(tempdis12)+4.0*EPSILON*tempdis6*(tempdis6-1.0)+EPSILON;
		//return 4.0*EPSILON*tempdis6*(tempdis6-1.0)+EPSILON;
	} else {
		return -paraK/2.0*bondlen2*log(tempdis12);
		//return 0.0;
	}*/
	/*tempdis12=sqrt(dMd)/bondlen;
	tempdis2=dMd/bondlen2;sfsf
	tempdis6=dMd*dMd*dMd;
	tempdis6=SIGMA6/tempdis6;
	return EPSILON*((-207.12+342.88*tempdis12-163.52*tempdis2+24.32*tempdis12*tempdis2)+tempdis6*(tempdis6-2.0));*/
	/*tempdis6=dMd*dMd*dMd;
	tempdis6=SIGMA6/tempdis6;
	tempdis2=sqrt(dMd)-bondlen;
	return paraK*tempdis2*tempdis2+EPSILON*tempdis6*(tempdis6-2.0);*/
	tempdis2=sqrt(dMd);
	tempdis6=1.0-exp(-lambda*(tempdis2-SIGMA));
	tempdis2=tempdis2-bondlen;
	tempdis2=tempdis2*tempdis2/0.25;
	if(tempdis2<1.0) {
		return -0.125*paraK*log(1.0-tempdis2)+EPSILON*(tempdis6*tempdis6-1.0);
	} else {
		return _MAX_DOUBLE;
	}
}
//#define Energy_AG(paraK,aMa2) paraK*aMa2/2.0
/*inline const double&  Energy_AG(const double&  paraK, const double&  cosaMa2, const double&  sinaMa2, const double&  costhe0, const double&  sinthe0) {
	return paraK*(1-cosaMa2*costhe0-sinaMa2*sinthe0);
}
inline const double&  Energy_DH(const double&  paraK, const double&  cosaMa2, const double&  sinaMa2, const double&  costhe0, const double&  sinthe0) {
	return paraK*(1+cosaMa2*costhe0+sinaMa2*sinthe0);//gromacs
}*/
inline double  Energy_AG(const double&  paraK, const double&  aMa2) {
	return paraK*(1-cos(aMa2));
}
inline double  Energy_DH(const double&  paraK, const double&  aMa2) {
	return paraK*(1+cos(aMa2));//gromacs
}
//{paraK*SQUARE(d-d0)/2.0} dMd2=square(d-d0)

//#define Energy_LJ_cut(EPSILON,SIGMA_TWELVE,SIGMA_HEXIAL,d_2,E_cut) {4.0*(fabs(EPSILON)*SIGMA_TWELVE/HEXIAL(d_2)-EPSILON*SIGMA_HEXIAL/TRIPLE(d_2))-E_cut}
//{//pay attention d_2=d*d; 
	//if( d_2 <= R_cut_2 ) { //if moved to the main energy calculation process;
		//; //LJ-1
	//} else {
	//	return 0.0;
	//}
	//return 4.0*( HEXIAL(SQUARE(SIGMA)/d_2) - EPSILON*TRIPLE(SQUARE(SIGMA)/d_2) ); //LJ-2
	//return d_2>SQUARE(SIGMA)?_MAX_DOUBLE:0.0; //hard-sphere;
//}
//#template<class T> inline T
//#define Energy_LJ(EPSILON,SIGMA12,SIGMA6,d12,d6) 4.0*(fabs(EPSILON)*SIGMA12/d12-EPSILON*SIGMA6/d6)

inline double  Energy_LJ(const int&  type, const double&  lambda, const double&  EPSILON,
	const double&  SIGMA12, const double&  SIGMA6, const double& SIGMA2, const double& SIGMA, 
	const double& d2) {
	//return 4.0*(fabs(EPSILON)*SIGMA12/d12-EPSILON*SIGMA6/d6);
	if(type==1) { // hard sphere
		if(d2<SIGMA2) {
			return _EMAX_DOUBLE;
		} else if(d2<lambda*lambda*SIGMA2) {
			return -EPSILON;
		} else {
			return 0.0;
		}
	} else if(type==2) { // normal lj
		tempdis6=d2*d2*d2;
		//tempdis12=tempdis6*tempdis6;
		tempdis6=SIGMA6/tempdis6;
		return EPSILON*tempdis6*(tempdis6-2.0); 
		//cout<<" tdis6="<<tempdis6<<" eps="<<EPSILON<<endl;
		//return 4.0*tempdis6*(tempdis6-EPSILON); //for PRL benchmark;
	} else if(type==3) { // morse potential
		tempdis6=1.0-exp(-lambda*(sqrt(d2)-SIGMA));
		return EPSILON*(tempdis6*tempdis6-1.0);
	} else if(type==4) { // 24-12 LJ
		tempdis6=d2*d2*d2;
		tempdis12=SIGMA12/tempdis6/tempdis6;
		return EPSILON*tempdis12*(tempdis12-2.0); 
	} else {
		return 0.0;
	}
}
//pay attention, only value can be used in it, not any proc or function;
/*inline const double&  Energy_LJ(const double&  ep, const double&  s12, const double&  s6, const double&  d12, const double&  d6) {
	return 4.0*( fabs(ep)*s12/d12 - ep*s6/d6 );
}*/
//{//pay attention d_2=d*d; 
	//return ; //LJ-1
	//return 4.0*( HEXIAL(SQUARE(SIGMA)/d_2) - EPSILON*TRIPLE(SQUARE(SIGMA)/d_2) ); //LJ-2
	//return d_2>SQUARE(SIGMA)?_MAX_DOUBLE:0.0; //hard-sphere;
//}
/*inline T Energy_LJ(T EPSILON, T SIGMA, T d_2)//pay attention d_2=d*d;
{
	//return 4.0*(fabs(EPSILON)*HEXIAL(SQUARE(SIGMA)/d_2)-EPSILON*TRIPLE(SQUARE(SIGMA)/d_2)); //LJ-1
	return 4.0*( HEXIAL(SQUARE(SIGMA)/d_2) - EPSILON*TRIPLE(SQUARE(SIGMA)/d_2) ); //LJ-2
	//return d_2>SQUARE(SIGMA)?_MAX_DOUBLE:0.0; //hard-sphere;
}
*/
template<class T>
inline T Energy_LJ_SURF(T EPSILON, T SIGMA, T d) {
	//return 4.0*( fabs(EPSILON)*HEXIAL(SQUARE(SIGMA/d)) - EPSILON*TRIPLE(SQUARE(SIGMA/d)) ); //LJ-2
	return fabs(EPSILON)*TRIPLE(TRIPLE(SIGMA/d))/7.5 - EPSILON*TRIPLE(SIGMA/d); //LJ-3
}
/*template<class T>
inline T Energy_LJ_SURF(T EPSILON, T SIGMA, T d)
{
	//return 4.0*( fabs(EPSILON)*QUIQUE(SQUARE(SIGMA/d)) - EPSILON*QUADRI(SIGMA/d) ); //LJ-1
	return 4.0*( QUIQUE(SQUARE(SIGMA/d)) - EPSILON*QUADRI(SIGMA/d) ); //LJ-2
	//return d>SIGMA?_MAX_DOUBLE:0.0; //hard-sphere;
}*/

template<class T> inline T Energy_ELE(T Q1, T Q2, T d) { return Q1*Q2/d; }


#endif


