//cmc.h

#ifndef _CMC_H 
#define _CMC_H

//#include <omp.h>
#include <mpi.h>
#include "cmolecule.h"
#include "energy.h"
#include <sstream>
#include "rwpara.h"
using namespace MPI;

#define _OUTPUT_LEN  50
#define _CONVERGENCE 1e-12
#define _CONF_DIR {string("_conformations")}
#define _CN_coeff 1.35

class cmc {
public:
	cmc();
	~cmc();
	void load_parameters(const string FILENAME_para);
	void write_parameters(const string FILENAME_para);// and bcast;
	void broadcast_parameters();
	//void make_dir();
	void init_conformation_a(const int NUM_atoms, const int NUM_chains, const double LEN_bond, const string FILENAME_conf,const double xorigin, const double yorigin, const double zorigin);
	void initialization();
	///////////////////////////////////////////////////////////////////////////
	void load_conformation(); //'a' or 'f'	
	void load_temperatures();
	void scatter_temperatures();
	///////////////////////////////////////////////////////////////////////////
	void tell_procid();
	void memo_allocation();
	bool _Memory_freed;
	void memo_setzero();    
	void memo_evaluation_fnode();
	void memo_evaluation_bcast();
	void memo_evaluation();
	inline void check_memo();
	void init_epsilonsigma();
	void broadcast_epsilonsigma();
	///////////////////////////////////////////////////////////////////////////
    double get_energy();                       //should = _ener_total
	double calc_energy();                      //the total energy; 
	void check_energy();
	void   init_energy();                      //doing before move;
	bool _Enery_initialized;
	void   calc_energy_each_atom();      
	inline void   calc_energy_chosen_atom();          //before & after each atom change; <1+<2
	inline void init_statistic();
	inline void close_statistic();
	bool _Statistic_over;
	void ensembler();
	void output_statistic();
	inline void make_choice(const int INDEX_coor);
	//inline void prep_change(); //now embedded in make_choice;
	inline bool make_change();

	inline void  make_judge();
	inline void  make_judge_mc();

	inline void make_accept();
	inline void make_reject_all();	
	inline void make_reject_only_move();	

	inline int  make_mcmove(const int INDEX_coor);
	inline bool ener_dispatch();
	inline int  make_mcmove_mc(const int INDEX_coor);

	int _RUNTIMES_recording;
	char* _RecordingNAME;
	inline void recording();
	bool _RESTARTFlag;
	double _ogboxlx;
	double _ogboxly;
	double _ogboxlz;
	inline void loadrestart(double ogboxlx, double ogboxly, double ogboxlz);
	//inline void trajectory_rec();
	void make_mapping(); // from _XX, _YY, _ZZ to the vector of atoms.
	void make_mapping_wrap(); // from _XX, _YY, _ZZ to the vector of atoms.
	bool chck_bond_len();
    void fout_conformation(int CHCK_BOND_LEN_OR_NOT);
	void memo_free();
	//////////////////////////////////////////////////////////////////////////
	CMolecule _system_;
	string _FILENAME_conf;
	bool _IFVERBOSE;
	
	int _ACCINDEX;
	//int _ACCINDEX_bf;
	//int _ACCINDEX_ag;
	//CMyArray3<double> _EAG_ACC;
	//CMyArray3<double> _EBF_ACC;
	CMyArray3<double> _ELJ_ACC;
	//CMyArray3<double> _EDH_ACC;
	//double _AG_interval;
	//CMyArray<double> _BF_interval;
	CMyArray<double> _LJ_interval;
	int index_lj;
	//int bf_index;
	//double _DH_interval;

	int _NUM_atoms;
	int _NUM_residues;
	int _NUM_chains;
	//int** _Series_eachain;//for statistics;
	int _NUM_chains_real;
	int _INDEX_chosen;
	double _ENER_total;
	int _NUM_replicas;
	double* _T_rep_eachrep;      // each replica
	double* _T_eachrep;          // each replica
	int*    _Index_T_rep;        // each replica, only on node0, index of which rep the T is in;
	int*    _Index_rep_T;        // each replica, every node,    index of which T the rep is;
	double* _E_rep;      // each replica
	int*    _NUM_rem;
	int*    _NUM_rem_acc;
	double _B_delta;
	int _INDEX_TEMPERATURE;

	double _E_lowest;    //       -20
	double _E_interval;  //        0.04
	double _E_highest;   // ....
	int _E_totalnum;  //        5000

	double _RG_lowest;    //       -20
	double _RG_interval;  //        0.04
	double _RG_highest;   // ....
	int _RG_totalnum;  //        5000

	CMyArray<double> _Probability;
	CMyArray<double> _Probability_tot;
	double* _Probability_all;
private:
	///////////////for test////////////////////////
	//int enumer;
	///////////////////////////////////////////////
	int _SIZE_memo;
	///////////////for bond fluctuation ///////////
	bool _BF_flag;
	double _BF_coeff;
	double _BF_scale;
	bool _AG_flag;
	bool _DH_flag;
	CMyArray<double> _BOND_length;
	CMyArray<double> _BOND_length2;
	CMyArray<double> _PARA_KB;
	CMyArray<double> _BOND_delta;
	CMyArray<double> _RMIN;
	CMyArray<double> _RMAX;
	CMyArray<double> _RMIN2;
	CMyArray<double> _RMAX2;
	CMyArray<double> _DDELTA;
	CMyArray<double> _DDELTA_ea;
	double tempddelta;
	int temprandom;
	CMyArray<double> _PARA_KA;
	CMyArray<double> _THETAZERO;
	CMyArray<double> _THETAZERO_cos;
	CMyArray<double> _THETAZERO_sin;
	CMyArray<double> _PARA_KD;
	CMyArray<double> _PHIZERO;
	CMyArray<double> _PHIZERO_cos;
	CMyArray<double> _PHIZERO_sin;
	CMyArray<double> _EPSILON_eachres;
	CMyArray<double> _LAMBDA_eachres;
	CMyArray<int> _PPTYPE_eachres;
	CMyArray<double> _SIGMA_eachres;
	CMyArray<double> _SIGMA2_eachres;
	CMyArray<double> _SIGMADIS_eachres;
	CMyArray<double> _SIGMAWELL_eachres;
	CMyArray<double> _SIGMAWELL2_eachres;
	CMyArray<double> _SIGMA3_eachres;
	CMyArray<double> _SIGMA6_eachres;
	CMyArray<double> _SIGMA9_eachres;
	CMyArray<double> _SIGMA12_eachres;
	CMyArray<double> _SIGMA24_eachres;
	CMyArray<double> _E_cut_RR_eachres;
	CMyArray<double> _R_cut_RR_eachres;
	//CMyArray<int>    _IF_RR_eachatom;
	//int*             _IF_RR_calc_or_not;

	string _FFPARA_FN;

	//for use in calc_ener_each_atom;
	double* cos_angle;
	double* sin_angle;
	int tempind;
	int tempind_x;
	int tempind_y;
	//for use in calc_ener_each_atom;

	//for judge
	//int tempindex_judge_bak;
	int tempindex_judge;
	int tempnumer_ret;
	double TempEner;
	//for judge

	int     _NEIGHBOR_lj;

	//double  _epsilon_surf;
	//double  _sigma_surf;
	//double  _sigma_surf_1_6;
	
	//bool _ITERATION_OR_NOT;
	//bool _WHAM_EACH_OR_NOT;
	//bool _ENSEMBLE_AVERAGE;
	double  _CoorX1;
	double  _CoorX2;
	double  _CoorY1;
	double  _CoorY2;
	double  _CoorZ1;
	double  _CoorZ2;
	double  _PBL_X;
	double  _PBL_Y;
	double  _PBL_Z;
	double  _PBL_X_2;
	double  _PBL_Y_2;
	double  _PBL_Z_2;
	int     _PBC_Dim;
	double  _DIS_x;
	double  _DIS_y;
	double  _DIS_z;
	////////////////////////////////////////
	template<class T>
	inline T Coor_PBC_X(T coor_x) {
		return coor_x>_CoorX1?( coor_x-_PBL_X*int((coor_x-_CoorX1)/_PBL_X) ) : ( coor_x-_PBL_X*int((coor_x-_CoorX2)/_PBL_X) );
	}
	template<class T>
	inline T Coor_PBC_Y(T coor_y) {
		return coor_y>_CoorY1?( coor_y-_PBL_Y*int((coor_y-_CoorY1)/_PBL_Y) ) : ( coor_y-_PBL_Y*int((coor_y-_CoorY2)/_PBL_Y) );
	}
	template<class T>
	inline T Coor_PBC_Z(T coor_z) {
		return coor_z>_CoorZ1?( coor_z-_PBL_Z*int((coor_z-_CoorZ1)/_PBL_Z) ) : ( coor_z-_PBL_Z*int((coor_z-_CoorZ2)/_PBL_Z) );
	}
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_X(T C1, T C2) {//c2-c1;
		T coor2=Coor_PBC_X(C2);
		T coor1=Coor_PBC_X(C1);
		T tempDIS=coor2-coor1;
		if( fabs(tempDIS) < (_PBL_X_2) ) {
			return tempDIS;
		} else {
			return coor2>coor1?(tempDIS-_PBL_X):(tempDIS+_PBL_X);
		}
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Y(T C1, T C2) {//c2-c1;
		T coor2=Coor_PBC_Y(C2);
		T coor1=Coor_PBC_Y(C1);
		T tempDIS=coor2-coor1;
		if( fabs(tempDIS) < (_PBL_Y_2) ) {
			return tempDIS;
		} else {
			return coor2>coor1?(tempDIS-_PBL_Y):(tempDIS+_PBL_Y);
		}
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Z(T C1, T C2) {//c2-c1;
		T coor2=Coor_PBC_Z(C2);
		T coor1=Coor_PBC_Z(C1);
		T tempDIS=coor2-coor1;
		if( fabs(tempDIS) < (_PBL_Z_2) ) {
			return tempDIS;
		} else {
			return coor2>coor1?(tempDIS-_PBL_Z):(tempDIS+_PBL_Z);
		}
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_X(T tempDIS) {//c2-c1;
		tempDIS=tempDIS-int(tempDIS/_PBL_X)*_PBL_X;
		if( fabs(tempDIS)>(_PBL_X_2) ) {
			tempDIS=tempDIS>0.0?(tempDIS-_PBL_X):(tempDIS+_PBL_X);
		}
		return tempDIS;
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Y(T tempDIS) {//c2-c1;
		tempDIS=tempDIS-int(tempDIS/_PBL_Y)*_PBL_Y;
		if( fabs(tempDIS)>(_PBL_Y_2) ) {
			tempDIS=tempDIS>0.0?(tempDIS-_PBL_Y):(tempDIS+_PBL_Y);
		}
		return tempDIS;
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Z(T tempDIS) {//c2-c1;
		tempDIS=tempDIS-int(tempDIS/_PBL_Z)*_PBL_Z;
		if( fabs(tempDIS)>(_PBL_Z_2) ) {
			tempDIS=tempDIS>0.0?(tempDIS-_PBL_Z):(tempDIS+_PBL_Z);
		}
		return tempDIS;
	}
	/////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void translate2box();

	//int    _INDEX_TEMPERATURE;
	double _TEMPERATURE;
	double _TEMPERATURE_REP;

	double _tempDH;
	double _tempphi0;
	double _tempkd;
	inline void dihedral2H1(int index) {
		_tempDH=dihedral2(_DIS_x_eachatom.pArray[index][index+1], 
			              _DIS_y_eachatom.pArray[index][index+1],
			              _DIS_z_eachatom.pArray[index][index+1],
			             _DIS_x_eachatom.pArray[index+1][index+2], 
			             _DIS_y_eachatom.pArray[index+1][index+2],
			             _DIS_z_eachatom.pArray[index+1][index+2],
			             _DIS_x_eachatom.pArray[index+2][index+3], 
			             _DIS_y_eachatom.pArray[index+2][index+3],
			             _DIS_z_eachatom.pArray[index+2][index+3] );
		_tempphi0=(_PHIZERO.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;

		_tempkd=( _PARA_KD.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;
		//cout<<" phi="<<_tempDH<<" e="<<_tempkd*(1+cos(_tempDH-_tempphi0))<<endl;
	}
	inline void dihedral2H2(int index) {
		_tempDH=dihedral2(-_DIS_x_eachatom.pArray[index+1][index], 
			              -_DIS_y_eachatom.pArray[index+1][index],
			              -_DIS_z_eachatom.pArray[index+1][index],
			             _DIS_x_eachatom.pArray[index+1][index+2], 
			             _DIS_y_eachatom.pArray[index+1][index+2],
			             _DIS_z_eachatom.pArray[index+1][index+2],
			             _DIS_x_eachatom.pArray[index+2][index+3], 
			             _DIS_y_eachatom.pArray[index+2][index+3],
			             _DIS_z_eachatom.pArray[index+2][index+3] );
		_tempphi0=(_PHIZERO.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;

		_tempkd=( _PARA_KD.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;
		//cout<<" phi="<<_tempDH<<" e="<<_tempkd*(1+cos(_tempDH-_tempphi0))<<endl;
	}
	inline void dihedral2H3(int index) {
		_tempDH=dihedral2(_DIS_x_eachatom.pArray[index][index+1], 
			              _DIS_y_eachatom.pArray[index][index+1],
			              _DIS_z_eachatom.pArray[index][index+1],
			             -_DIS_x_eachatom.pArray[index+2][index+1], 
			             -_DIS_y_eachatom.pArray[index+2][index+1],
			             -_DIS_z_eachatom.pArray[index+2][index+1],
			             _DIS_x_eachatom.pArray[index+2][index+3], 
			             _DIS_y_eachatom.pArray[index+2][index+3],
			             _DIS_z_eachatom.pArray[index+2][index+3] );
		_tempphi0=(_PHIZERO.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;

		_tempkd=( _PARA_KD.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;
		//cout<<" phi="<<_tempDH<<" e="<<_tempkd*(1+cos(_tempDH-_tempphi0))<<endl;
	}
	inline void dihedral2H4(int index) {
		_tempDH=dihedral2(_DIS_x_eachatom.pArray[index][index+1], 
			              _DIS_y_eachatom.pArray[index][index+1],
			              _DIS_z_eachatom.pArray[index][index+1],
			             _DIS_x_eachatom.pArray[index+1][index+2], 
			             _DIS_y_eachatom.pArray[index+1][index+2],
			             _DIS_z_eachatom.pArray[index+1][index+2],
			            -_DIS_x_eachatom.pArray[index+3][index+2], 
			            -_DIS_y_eachatom.pArray[index+3][index+2],
			            -_DIS_z_eachatom.pArray[index+3][index+2] );
		_tempphi0=(_PHIZERO.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				 +_PHIZERO.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;

		_tempkd=( _PARA_KD.pArray[_INDEX_RES_ATM[index]][_INDEX_RES_ATM[index+1]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+1]][_INDEX_RES_ATM[index+2]]
				+_PARA_KD.pArray[_INDEX_RES_ATM[index+2]][_INDEX_RES_ATM[index+3]])/3.0;
		//cout<<" phi="<<_tempDH<<" e="<<_tempkd*(1+cos(_tempDH-_tempphi0))<<endl;
	}

	//char* _OUTPUTDIR; 
	/* when rem, the following may need change;
	FILE* _OUTfile;
	char* _FNoutput;
	vector<string> _SAVElist;
	int _NUM_SAVElist;*/
	
	/*double* _PARA_beta;
	double* _PARA_alpha;
	double* _PARA_beta_new;
	double* _PARA_alpha_new;*/
////////////////////////////////
	CMyArray<double> _ENER_LJ_eachatom; //L-J potential;
	double* _ENER_LJ_chosenatom_backup;
	double* _ENER_BF_eachatom; //bond fluctuation;
	double  _ENER_BF_chosenatom_backup_l;
	double  _ENER_BF_chosenatom_backup_r;
	double* _ENER_AG_eachatom; //bond fluctuation;
	double  _ENER_AG_chosenatom_backup_l;
	double  _ENER_AG_chosenatom_backup_c;
	double  _ENER_AG_chosenatom_backup_r;
	double* _ENER_DH_eachatom; //bond fluctuation;
	double  _ENER_DH_chosenatom_backup_1;
	double  _ENER_DH_chosenatom_backup_2;
	double  _ENER_DH_chosenatom_backup_3;
	double  _ENER_DH_chosenatom_backup_4;
	/*CMyArray<double> _ENER_AG_eachatom; //bond fluctuation;
	double*  _ENER_AG_chosenatom_backup;
	CMyArray<double> _ENER_DH_eachatom; //dihedral;
	double*  _ENER_DH_chosenatom_backup;*/
	CMyArray<double> _DIS2_eachatom; // in one box
	double* _DIS2_chosenatom_backup; // in one box
	CMyArray<double> _DIS_x_eachatom; // in one box
	double* _DIS_x_chosenatom_backup; // in one box
	CMyArray<double> _DIS_y_eachatom; // in one box
	double* _DIS_y_chosenatom_backup; // in one box
	CMyArray<double> _DIS_z_eachatom; // in one box
	double* _DIS_z_chosenatom_backup; // in one box
	double* _COM_x;
	double* _COM_y;
	double* _COM_z;
	double* _DISSTAT;
	double _DISSTAT_lowest;
	double _DISSTAT_highest;
	double _DISSTAT_interval;
	double* _RG2_x;
	double* _RG2_y;
	double* _RG2_z;
	double* _RG2_ec;
	CMyArray<double> _COM_x_stat;
	CMyArray<double> _COM_y_stat;
	CMyArray<double> _COM_z_stat;
	CMyArray3<double> _DIS_stat;
	CMyArray<double> _DIS_stat_actual;
	CMyArray<double> _DIS_stat_actual_tot;

	int** _CINDEXMAP;
	CMyArray<double> _RG2_actual_x;
	CMyArray<double> _RG2_actual_y;
	CMyArray<double> _RG2_actual_z;
	CMyArray<double> _RG2_actual_xtot;
	CMyArray<double> _RG2_actual_ytot;
	CMyArray<double> _RG2_actual_ztot;
	CMyArray3<double> _RG2_stat;
	CMyArray<double> _EINT_stat;
	CMyArray<double> _CN_backup;
	double _CNSTAT;
	double* _CN_stat;
	double* _CN_stat_tot;
	int _EINT_totalnum;
	double _EINT_interval;
	double _EINT_lowest;
	double _EINT_highest;
	int* _EINTlist;
    int _len_EINT;
    int _EINT_index;
	int* _indexRGSTAT;
	int tempindexstat;
	int tempindexstat_ij;
	double tempnumerstat_x;
	double tempnumerstat_y;
	double tempnumerstat_z;
	CMyArray<double> _COM_x_stat_tot;
	CMyArray<double> _COM_y_stat_tot;
	CMyArray<double> _COM_z_stat_tot;

	CMyArray3<double> _DIS_stat_tot;
	/*CMyArray3<double> _RG2_x_stat_tot;
	CMyArray3<double> _RG2_y_stat_tot;
	CMyArray3<double> _RG2_z_stat_tot;*/
	CMyArray3<double> _RG2_stat_tot;
	CMyArray<double> _EINT_stat_tot;
	//CMyArray<double> _ENERAA_stat;
	CMyArray<double> _ENERLJ_stat;
	CMyArray<double> _ENERBF_stat;
	CMyArray<double> _ENERAG_stat;
	CMyArray<double> _ENERDH_stat;


	CMyArray<double> _ENERLJ_stat_tot;
	CMyArray<double> _ENERBF_stat_tot;
	CMyArray<double> _ENERAG_stat_tot;
	CMyArray<double> _ENERDH_stat_tot;
    char* tcharn;

    int* _OPlist;
    int _len_OP;
    vector<int> _Runninglist;
    int _realSize;
    //int _tfrom;
    //int _tto;
	bool _FLAG_rg2;
	bool _FLAG_dis;
	//double* _rgyration2; //for each selection...
	//ofstream* _rg2_stream;
	//////////////
	bool _FLAG_ener;
	double _FLAG_pivot;
	double _numerpivot;
	//double* _rhead2;     //for each selection...
	//ofstream* _rh2_stream;
	//ofstream _ener_stream;
	inline void statistic();
	inline void getrealrate();
	
	int stat_i;
	int stat_j;
	int stat_k;
	int stat_l;
	int stat_size;
	int stat_head;
	int stat_tail; 
	int stat_head_j;
	int stat_tail_j; 
	double* stat_com_x;
	double* stat_com_y;
	double* stat_com_z;
	double tempnum;
	int* _INDEX_CHN_HEAD; //for each chain...
	int* _INDEX_CHN_TAIL; //for each chain...
	int _INDEX_chnhead;
	int _INDEX_chntail;
	int* _INDEX_CHN_HEAD_real; //for each chain...
	int* _INDEX_CHN_TAIL_real; //for each chain...
	int _INDEX_chnhead_real;
	int _INDEX_chntail_real;
	int* _INDEX_LNEIGHBOR;
	int* _INDEX_RNEIGHBOR;
	int* _INDEX_LNEIGHBOR_real;
	int* _INDEX_RNEIGHBOR_real;
	int _Index_lneighbor_ind;
	int _Index_rneighbor_ind;
	int _Index_lneighbor_ind_real;
	int _Index_rneighbor_ind_real;
	int _Res_chosen;

	double _DIS;
	double _DIS2;
	double _DIS3;
	double _DIS6;
	double _DIS9;
	double _DIS12;

	double*  tempElj;
	double*  tempEbf;
	double*  tempEag;
	double*  tempEdh;
	/////
	double  _ENER_delta;
	double* _ENER_delta_LJ_inter;
	double  _ENER_delta_LJ_intra;
	double  _ENER_delta_LJ;
	double  _ENER_delta_BF_l;
	double  _ENER_delta_BF_r;
	double  _ENER_delta_AG_l;
	double  _ENER_delta_AG_c;
	double  _ENER_delta_AG_r;
	double  _ENER_delta_DH_1;
	double  _ENER_delta_DH_2;
	double  _ENER_delta_DH_3;
	double  _ENER_delta_DH_4;
	// the following 4 is for the chosen atom;
	int   _INDEX_chn_ind;           // for the chosen one;  a
	int   _SIZE_of_chosen_chn;      // for the chosen one;  b
	int   _SIZE_of_chosen_chn_real;      // for the chosen one;  b
	int   _INDEX_chn_or_not;        // for the chosen one;  c
	int   _INDEX_atm_in_chn_ind;    // for the chosen one;  d
	int   _INDEX_atm_in_chn_ind_real;    // for the chosen one;  d
	int   _TYPE_atom_ind;           // for the chosen one;  e
	int   _TYPE_atom_ind_real;      // for the chosen one;  f chain move only.
	//void  memoset_fnode();
	int*  _INDEX_chn_atom;         // a
	int*  _SIZE_of_chn;            // b
	int*  _SIZE_of_chn_real;            // b
	int*  _INDEX_chn_or_not_atom;  // c.1
	int*  _INDEX_chn_or_not_chain; // c.2
	int*  _INDEX_atm_in_chn;  // d
	int*  _INDEX_atm_in_chn_real;  // d
	int*  _TYPE_atom;         // e
	int*  _TYPE_atom_real;         // e
	int*  _INDEX_RES_ATM;  // f.1
	//int*  _INDEX_RES_RES;  // f.2
	////////////////////////////////////////////////////////
	///////////////for rotation calculation ////////////////
	/*double _rmin_l;
	double _rmax_l;
	double _ddelta_l;
	double _rmin_r;
	double _rmax_r;
	double _ddelta_r;*/
	double r2min;
	double r2max;
	double ee;
	double ff;
	double dd;
	double dd_sqrt;
	double* ec;
	//double* fc;
	double* dc;
	double* vecDC;
	double len_vecDC;//for _BF
	double* coorRES;
	double factor;
	double* tc;
	double tl;
	double planeD;
	double _PI_single;
	double _PI_double;
	double _PI_half;
	double Omega;
	double THETA;
	double uc;
	double vc;
	double wc;
	double TempBX;
	double TempBY;
	double TempBZ;

	double sin_PHI;
	double cos_PHI;
	//double sin_THETA_B;
	//double cos_THETA_B;
	double Len_CB;
	double Len_CB_xy;
	double Len_CB_xy_2;
	double sin_Omega;
	double cos_Omega;
	double sin_THETA;
	double cos_THETA;
	/////////////// ends here //////////////////////////////
public:
	//void make_preparation();
	//void make_conformations();
	////////////////// mpisub //////////////////
	void init_mpi(int pargc, char **pargv);              // init mpi or not;
	void term_mpi();
	//void init_parameters();                            // init and bcast parameters;
	//void init_conformations();  //initialization;
	//void load_conformations();  //just load conformations;
	void run_with_stepnumber(const int stepnumber);
	void run();
	void run_rem();
	inline void   rand_update(const int procid);

private:

	double* _XX;
	double* _YY;
	double* _ZZ;
	CMyArray<double> _XX_allrep;
	CMyArray<double> _YY_allrep;
	CMyArray<double> _ZZ_allrep;
	double* _XX_rec;
	double* _YY_rec;
	double* _ZZ_rec;
	//int _NUM_rec;

	int  _RUNTIMES_eachstep;         // 100000
    int  _RUNTIMES_totalnum;         // 10000 (10000 * 100000)
	int  _RUNTIMES_output;           //50
	int  _RUNTIMES_remgap;
	//bool _MPI_OR_NOT;

	int _PROC_size;
	int _PROC_ID;
	//int  _RUNTIMES_iteration;
	//bool _START_fromzero;
	//int  _ITER;

	int _I_eachstep;
	int _I_totalnum;
	
 	////////////////  rand  ////////////////
	inline void   rand_init(const int procid);
	inline double rand_seed(int &iseed);
	//inline void   rand_update(const int procid);
	int *ir;
	int iy;
	int jj;

	double ran;

	double _onethird1;
	double _onethird2;
	double _movelength;
	double* _movelength_ea;
	//double* _movelength_pv;
	double _movelength_temp;
	int iseed_zero;
	int iseed_len1;
	int iseed_len2;
	int iseed_len3;
	int iseed_angle;
	int iseed_index;
	int iseed_rand;	

	double _MC_NUM_TOT; //total;
	double _MC_NUM_SUC; //succeed; ==1
	double* _MC_NUM_TOT_stat; // total;
	double* _MC_NUM_SUC_stat; // succeed; ==1
	double* _MC_NUM_TOT_all; // total;
	double* _MC_NUM_SUC_all; // succeed; ==1
	double _MC_NUM_FIL; //fail;    -1
	double* _MC_NUM_FIL_stat; //fail;    -1
	double* _MC_NUM_FIL_all; //fail;    -1
	double _MC_NUM_STT; //statistic;
	//void fout_ratio();
	//void cout_parameters();


};

////////////////////////////////////
#define mm 714025
#define ia 1366
#define iccc 150889
#define rm 1.0/mm
void cmc::rand_init(const int procid) {
	//different process has different seeds
	iseed_zero=(-2)*(procid+1)-1;
	iseed_len1=(-14)*(procid+1)-1;
	iseed_len2=(-26)*(procid+1)-1;
	iseed_len3=(-34)*(procid+1)-1;
	iseed_angle=(-6)*(procid+1)-1;
	iseed_index=(-10)*(procid+1)-1;
	iseed_rand=(-22)*(procid+1)-1;
}
void cmc::rand_update(const int procid) {
    iseed_zero=(-2)*(procid+1)-1;
	iseed_len1=(-14)*(procid+1)-1;
	iseed_len2=(-26)*(procid+1)-1;
	iseed_len3=(-34)*(procid+1)-1;
	iseed_angle=(-6)*(procid+1)-1;
	iseed_index=(-10)*(procid+1)-1;
	iseed_rand=(-22)*(procid+1)-1;
	MEMOSETZERO(ir, sizeof(int)*98);
	iy=0;
	jj=0;
}
double cmc::rand_seed(int &iseed) {
	if(iseed<0) {
		iseed=(int)fmod(double(iccc-iseed), double(mm));
		for(int j=1; j<=97; j++) {
			iseed=(int)fmod(double(ia*iseed+iccc), double(mm));
			ir[j]=iseed;
		}
		iseed=(int)fmod(double(ia*iseed+iccc), double(mm));
		iy=iseed;
	}
	jj=1+(97*iy)/mm;
	iy=ir[jj];
	ran=iy*rm;
	iseed=(int)fmod(double(ia*iseed+iccc), double(mm));
	ir[jj]=iseed;
	return ran;
}

/////////////////////////////////
#endif

