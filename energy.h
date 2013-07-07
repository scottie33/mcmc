//energy.h

#ifndef _ENERGY_H
#define _ENERGY_H

#define KBT         1.00
#define e		   -1.00
#define E_x			1.00
#define E_y			0.00
#define E_z			0.00

//#define Energy_BF(paraK,dMd2) paraK*dMd2/2.0
inline double Energy_BF(double paraK, double dMd2) {
	return paraK*dMd2/2.0;
}
//#define Energy_AG(paraK,aMa2) paraK*aMa2/2.0
/*inline double Energy_AG(double paraK, double cosaMa2, double sinaMa2, double costhe0, double sinthe0) {
	return paraK*(1-cosaMa2*costhe0-sinaMa2*sinthe0);
}
inline double Energy_DH(double paraK, double cosaMa2, double sinaMa2, double costhe0, double sinthe0) {
	return paraK*(1+cosaMa2*costhe0+sinaMa2*sinthe0);//gromacs
}*/
inline double Energy_AG(double paraK, double aMa2) {
	return paraK*(1-cos(aMa2));
}
inline double Energy_DH(double paraK, double aMa2) {
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
inline double Energy_LJ(double EPSILON, double SIGMA12, double SIGMA6, double d12, double d6) {
	//return 4.0*(fabs(EPSILON)*SIGMA12/d12-EPSILON*SIGMA6/d6);
	return 4.0*(fabs(EPSILON)*SIGMA12/d12-EPSILON*SIGMA6/d6); //prl
}
//pay attention, only value can be used in it, not any proc or function;
/*inline double Energy_LJ(double ep, double s12, double s6, double d12, double d6) {
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


