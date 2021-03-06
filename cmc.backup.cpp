//cmc.cpp

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#define HAVE_NO_VARIABLE_RETURN_TYPE_SUPPORT 1

//#include <mpi.h>

#include "cmc.h"

///////////////////////////xxxxxx
cmc::cmc() {		
	ee=0.0; ff=0.0; dd=0.0;
	dd_sqrt=0.0;
	tc=NULL; ec=NULL; //fc=NULL; 
	dc=NULL;
    vecDC=NULL; coorRES=NULL;
    len_vecDC=0.0;
	factor=0.0;
	tl=0.0;

	_PI_single=acos(-1);
	_PI_double=acos(-1.0)*2.0; _PI_half=acos(-1.0)/2.0;

	_onethird1=1.0/3.0;
	_onethird2=2.0/3.0;

	_movelength=0.1;
	_movelength_ea=NULL;
	//_movelength_pv=NULL;

	sin_PHI=0.0; cos_PHI=0.0;
	sin_THETA=0.0; cos_THETA=0.0;
	Len_CB=0.0; Len_CB_xy=0.0; Len_CB_xy_2=0.0;
	sin_Omega=0.0; cos_Omega=0.0;
	
	Omega=0.0; THETA=0.0;
	uc=0.0; vc=0.0; wc=0.0;
	TempBX=0.0; TempBY=0.0; TempBZ=0.0;
	_NUM_atoms=0; _NUM_chains=0; _NUM_chains_real=0;

	_ENER_LJ_chosenatom_backup=NULL;
	_ENER_total=0.0;
	_ENER_delta=0.0;
	_ENER_delta_LJ=0.0;
	_ENER_delta_LJ_intra=0.0;
	_ENER_delta_LJ_inter=NULL;
	_ENER_delta_BF_l=0.0;
	_ENER_delta_BF_r=0.0;
	_ENER_delta_AG_l=0.0;
	_ENER_delta_AG_c=0.0;
	_ENER_delta_AG_r=0.0;
	_ENER_delta_DH_1=0.0;
	_ENER_delta_DH_2=0.0;
	_ENER_delta_DH_3=0.0;
	_ENER_delta_DH_4=0.0;

	_Enery_initialized=false;
	_Memory_freed=true;

	_NEIGHBOR_lj=1;

	_TEMPERATURE=1.0;
	_TEMPERATURE_REP=1.0;

	_INDEX_chn_ind=0;        // for the chosen one;
	_INDEX_chn_or_not=false; // for the chosen one;
	_INDEX_atm_in_chn_ind=0;     // for the chosen one;
	_INDEX_atm_in_chn_ind_real=0;     // for the chosen one;
	_TYPE_atom_ind=0;        // for the chosen one;
	_TYPE_atom_ind_real=0;        // for the chosen one;
	
	_SIZE_of_chosen_chn=0;      // for the chosen one;
	_SIZE_of_chosen_chn_real=0;      // for the chosen one;
	_INDEX_chnhead=0;	_INDEX_chntail=0;

	_INDEX_chn_atom=NULL;
	_SIZE_of_chn=NULL;
	_SIZE_of_chn_real=NULL;
	_INDEX_CHN_HEAD=NULL;	_INDEX_CHN_TAIL=NULL;
	_INDEX_CHN_HEAD_real=NULL;	_INDEX_CHN_TAIL_real=NULL;
	_INDEX_chn_or_not_atom=NULL; 	_INDEX_chn_or_not_chain=NULL;
	_INDEX_atm_in_chn=NULL;   
	_INDEX_atm_in_chn_real=NULL;   
	_TYPE_atom=NULL; 
	_TYPE_atom_real=NULL; 
	_INDEX_LNEIGHBOR=NULL;	_INDEX_RNEIGHBOR=NULL;
	_Index_lneighbor_ind=0;	_Index_rneighbor_ind=0;
	_INDEX_LNEIGHBOR_real=NULL;	_INDEX_RNEIGHBOR_real=NULL;
	_Index_lneighbor_ind_real=0;	_Index_rneighbor_ind_real=0;
	_INDEX_RES_ATM=NULL;       

	_BF_flag=false; 
	_BF_coeff=1.0; 
	_BF_scale=1.0; 
	_AG_flag=false;
	_DH_flag=false;
	_OPlist=NULL;
	_len_OP=0;
	_EINTlist=NULL;
	_len_EINT=0;
	_FLAG_ener=false;
	_FLAG_pivot=1.0;
	_FLAG_rg2=false;
	_FLAG_dis=false;
	stat_i=0; stat_j=0; stat_k=0; stat_size=0; stat_head=0; stat_tail=0; 
	stat_com_x=NULL; stat_com_y=NULL; stat_com_z=NULL; tempnum=0.0;
	_SIZE_memo=0;
	_IFVERBOSE=false;
	_ACCINDEX=1000;

	_NUM_replicas=0;
	_T_rep_eachrep=NULL;      // each replica
	_T_eachrep=NULL;      // each replica
	_Index_T_rep=NULL;        // each replica, only on node0, index of which rep the T is in;
	_Index_rep_T=NULL;        // each replica, every node,    index of which T the rep is;
	_E_rep=NULL;      // each replica
	_NUM_rem=NULL;
	_NUM_rem_acc=NULL;
	_RUNTIMES_eachstep=0;        // 100000
    _RUNTIMES_totalnum=0;        // 10000
    _RUNTIMES_recording=0;     
	_RUNTIMES_output=0;          // 40 
	_B_delta=0.0;
	_INDEX_TEMPERATURE=0;

	_E_lowest=0.0;    //       -20
	_E_interval=0.0;  //        0.04
	_E_highest=0.0;
	_E_totalnum=0;  //        5000

	_EINT_lowest=0.0;    //       -20
	_EINT_interval=0.0;  //        0.04
	_EINT_highest=0.0;
	_EINT_totalnum=0;  //        5000
	
	_RG_lowest=0.0;    //       -20
	_RG_interval=0.0;  //        0.04
	_RG_highest=0.0;
	_RG_totalnum=0.0;  //        5000

	_RUNTIMES_remgap=0;
	//_MPI_OR_NOT=false;

	_I_eachstep=0; _I_totalnum=0;
	
	iy=0; jj=0; ran=0.0; ir=NULL;
	iseed_zero=-1; 
	iseed_len1=-1; iseed_len2=-1; iseed_len3=-1; 
	iseed_angle=-1; iseed_index=-1; iseed_rand=-1;

	//_RUNTIMES_iteration=0;
	//_START_fromzero=false;
	//_ITER=0;

	//_ddelta_l=0.0; _rmin_l=0.0; _rmax_l=0.0; 
	//_ddelta_r=0.0; _rmin_r=0.0; _rmax_r=0.0;
	r2min=0.0; r2max=0.0;
	
	//_DIS=0.0;
	_DIS2=0.0;_DIS3=0.0;_DIS6=0.0;_DIS9=0.0;_DIS12=0.0;
	_DIS_x=0.0; _DIS_y=0.0; _DIS_z=0.0;
	tempind=0; tempind_x=0; tempind_y=0;

	tempindex_judge=0; //tempindex_judge_bak=0;
	tempnumer_ret=-10; TempEner=0.0;

	_XX=NULL; _YY=NULL; _ZZ=NULL;
	_XX_rec=NULL; _YY_rec=NULL; _ZZ_rec=NULL;

	_MC_NUM_TOT=0;
	_MC_NUM_SUC=0;
	_MC_NUM_TOT_stat=NULL;
	_MC_NUM_SUC_stat=NULL;
	_MC_NUM_TOT_all=NULL;
	_MC_NUM_SUC_all=NULL;
	_MC_NUM_FIL=0;
	_MC_NUM_FIL_stat=NULL;
	_MC_NUM_FIL_all=NULL;
	_Statistic_over=true;
	//_MC_NUM_STT=0;
}
///////////////////////////
cmc::~cmc() {
	this->memo_free();
}
///////////////////////////
void cmc::tell_procid() {
	//printf(" @proc[%03d]", _PROC_ID);
	cout<<" @proc["<<setw(3)<<_PROC_ID<<"]";
}
void cmc::memo_allocation() {	
	//int i=0;
	string TempStr;
	string TempInterval;
	int i=0;
	tell_procid(); cout<<" memo_allocation start..."<<endl;
	///////////////   code start  //////////////////	
	ec=new double[3];
	//fc=new double[3];
	dc=new double[3];
	tc=new double[3];

	vecDC=new double[3];
	coorRES=new double[3];
	
	ir=new int[98];
	
	if(_NUM_atoms==0) {
		ErrorMSG(string("_NUM_atoms==0;::cmc::memo_allocation"));
		exit(SIZEERROR);
	}
	_T_rep_eachrep=new double[_NUM_replicas];      // each replica
	_T_eachrep=new double[_NUM_replicas];      // each replica
	_Index_T_rep=new int[_NUM_replicas];        // each replica, only on node0, index of which rep the T is in;
	_Index_rep_T=new int[_NUM_replicas];        // each replica, every node,    index of which T the rep is;
	_E_rep=new double[_NUM_replicas];      // each replica
	_NUM_rem=new int[_NUM_replicas];
	_NUM_rem_acc=new int[_NUM_replicas];

	_MC_NUM_SUC_stat=new double[_NUM_replicas];
	_MC_NUM_TOT_stat=new double[_NUM_replicas];
	_MC_NUM_SUC_all=new double[_NUM_replicas];
	_MC_NUM_TOT_all=new double[_NUM_replicas];

	_MC_NUM_FIL_stat=new double[_NUM_replicas];
	_MC_NUM_FIL_all=new double[_NUM_replicas];

	_movelength_ea=new double[_NUM_replicas];
	//_movelength_pv=new double[_NUM_replicas];

	_SIZE_memo=MEMOSIZE(_NUM_atoms);
	_ENER_LJ_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_ENER_LJ_chosenatom_backup=new double[_SIZE_memo];

	_ENER_BF_eachatom=new double[_SIZE_memo];
	_ENER_BF_chosenatom_backup_l=0.0;
	_ENER_BF_chosenatom_backup_r=0.0;

	_ENER_AG_eachatom=new double[_SIZE_memo+1];//important! for calculation convenience; the 0 and N+1 no use but for easy prog.
	_ENER_AG_chosenatom_backup_l=0.0;
	_ENER_AG_chosenatom_backup_c=0.0;
	_ENER_AG_chosenatom_backup_r=0.0;

	_ENER_DH_eachatom=new double[_SIZE_memo+2];//important! for calculation convenience; the 0 and N+1 no use but for easy prog.
	_ENER_DH_chosenatom_backup_1=0.0;
	_ENER_DH_chosenatom_backup_2=0.0;
	_ENER_DH_chosenatom_backup_3=0.0;
	_ENER_DH_chosenatom_backup_4=0.0;

	

	//_Series_eachain=new int*[_NUM_chains];

	_DIS2_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS2_chosenatom_backup=new double[_SIZE_memo];
	_DIS_x_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_x_chosenatom_backup=new double[_SIZE_memo];
	_DIS_y_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_y_chosenatom_backup=new double[_SIZE_memo];
	_DIS_z_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_z_chosenatom_backup=new double[_SIZE_memo];
	_COM_x=new double[_NUM_chains];
	_COM_y=new double[_NUM_chains];
	_COM_z=new double[_NUM_chains];
	_RG2_x=new double[_NUM_chains];
	_RG2_y=new double[_NUM_chains];
	_RG2_z=new double[_NUM_chains];
	_RG2_ec=new double[_NUM_chains];
	_RG2_actual_x.Build(_NUM_chains,_E_totalnum);
	_RG2_actual_xtot.Build(_NUM_chains,_E_totalnum);
	_RG2_actual_y.Build(_NUM_chains,_E_totalnum);
	_RG2_actual_ytot.Build(_NUM_chains,_E_totalnum);
	_RG2_actual_z.Build(_NUM_chains,_E_totalnum);
	_RG2_actual_ztot.Build(_NUM_chains,_E_totalnum);
	if( _NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) { 
		_DISSTAT=new double[(_NUM_chains-1)*_NUM_chains/2];
		_DIS_stat.Build((_NUM_chains-1)*_NUM_chains/2, _E_totalnum, _RG_totalnum);
		_DIS_stat_tot.Build((_NUM_chains-1)*_NUM_chains/2, _E_totalnum, _RG_totalnum);
		_DIS_stat_actual.Build((_NUM_chains-1)*_NUM_chains/2, _E_totalnum);
		_DIS_stat_actual_tot.Build((_NUM_chains-1)*_NUM_chains/2, _E_totalnum);
	}
	_CINDEXMAP=new int*[_NUM_chains];
	for(i=0; i<_NUM_chains; i++){
		_CINDEXMAP[i]=new int[_NUM_chains];
	}
	_indexRGSTAT=new int[_NUM_chains];
	stat_com_x=new double[_NUM_chains];
	stat_com_y=new double[_NUM_chains];
	stat_com_z=new double[_NUM_chains];
	
	_COM_x_stat.Build(_NUM_chains, _E_totalnum);
	_COM_y_stat.Build(_NUM_chains, _E_totalnum);
	_COM_z_stat.Build(_NUM_chains, _E_totalnum);
	
	/*_RG2_x_stat.Build(_NUM_chains, _E_totalnum, _RG_totalnum);
	_RG2_y_stat.Build(_NUM_chains, _E_totalnum, _RG_totalnum);
	_RG2_z_stat.Build(_NUM_chains, _E_totalnum, _RG_totalnum);*/
	_RG2_stat.Build(_NUM_chains, _E_totalnum, _RG_totalnum);

	if(_len_EINT>0 ) { 
		_CN_stat=new double[_E_totalnum];
		_CN_stat_tot=new double[_E_totalnum];
		if(_EINT_totalnum>1) {
			_EINT_stat.Build(_E_totalnum, _EINT_totalnum);
			_EINT_stat_tot.Build(_E_totalnum, _EINT_totalnum);
		}
	}

	_ENERLJ_stat.Build((_NUM_chains+1)*_NUM_chains/2, _E_totalnum);
	_ENERBF_stat.Build(_NUM_chains, _E_totalnum);
	_ENERAG_stat.Build(_NUM_chains, _E_totalnum);
	_ENERDH_stat.Build(_NUM_chains, _E_totalnum);
	_ENERLJ_stat_tot.Build((_NUM_chains+1)*_NUM_chains/2, _E_totalnum);
	_ENERBF_stat_tot.Build(_NUM_chains, _E_totalnum);
	_ENERAG_stat_tot.Build(_NUM_chains, _E_totalnum);
	_ENERDH_stat_tot.Build(_NUM_chains, _E_totalnum);
	_COM_x_stat_tot.Build(_NUM_chains, _E_totalnum);
	_COM_y_stat_tot.Build(_NUM_chains, _E_totalnum);
	_COM_z_stat_tot.Build(_NUM_chains, _E_totalnum);
	/*_DIS_x_stat_tot.Build(_NUM_chains, _E_totalnum, 50);
	_DIS_y_stat_tot.Build(_NUM_chains, _E_totalnum, 50);
	_DIS_z_stat_tot.Build(_NUM_chains, _E_totalnum, 50);*/
	
	/*_RG2_x_stat_tot.Build(_NUM_chains, _E_totalnum, _RG_totalnum);
	_RG2_y_stat_tot.Build(_NUM_chains, _E_totalnum, _RG_totalnum);
	_RG2_z_stat_tot.Build(_NUM_chains, _E_totalnum, _RG_totalnum);*/
	_RG2_stat_tot.Build(_NUM_chains, _E_totalnum, _RG_totalnum);
	
	//_ENERAB_stat.Build((_NUM_chains-1)*_NUM_chains/2, _E_totalnum);
	_ENER_delta_LJ_inter=new double[_NUM_chains];
	tempElj=new double[(_NUM_chains+1)*_NUM_chains/2];
	tempEbf=new double[_NUM_chains];
	tempEag=new double[_NUM_chains];
	tempEdh=new double[_NUM_chains];
	//_rgyration2=new double[_NUM_chains];
	//_rhead2=new double[_NUM_chains];

	_INDEX_chn_atom=new int[_SIZE_memo];
	_SIZE_of_chn=new int[_NUM_chains];
	_SIZE_of_chn_real=new int[_NUM_chains];
	_INDEX_CHN_HEAD=new int[_NUM_chains];
	_INDEX_CHN_TAIL=new int[_NUM_chains];
	_INDEX_CHN_HEAD_real=new int[_NUM_chains];
	_INDEX_CHN_TAIL_real=new int[_NUM_chains];
	_INDEX_atm_in_chn=new int[_SIZE_memo];
	_INDEX_atm_in_chn_real=new int[_SIZE_memo];
	_TYPE_atom=new int[_SIZE_memo];
	_TYPE_atom_real=new int[_SIZE_memo];
	_INDEX_LNEIGHBOR=new int[_SIZE_memo];
	_INDEX_RNEIGHBOR=new int[_SIZE_memo];
	_INDEX_LNEIGHBOR_real=new int[_SIZE_memo];
	_INDEX_RNEIGHBOR_real=new int[_SIZE_memo];
	_INDEX_RES_ATM=new int[_SIZE_memo];
    _XX=new double[_SIZE_memo];
	_YY=new double[_SIZE_memo];
	_ZZ=new double[_SIZE_memo];
	_XX_allrep.Build(_NUM_replicas, _SIZE_memo);
	_YY_allrep.Build(_NUM_replicas, _SIZE_memo);
	_ZZ_allrep.Build(_NUM_replicas, _SIZE_memo);
	/*_NUM_rec=_SIZE_memo*_RUNTIMES_eachstep;*/
	_XX_rec=new double[_SIZE_memo];
	_YY_rec=new double[_SIZE_memo];
	_ZZ_rec=new double[_SIZE_memo];
	_INDEX_chn_or_not_atom=new int[_SIZE_memo];
	_INDEX_chn_or_not_chain=new int[_NUM_chains];

	//tell_procid(); cout<<" I did it!!!"<<endl; 
	_PARA_KB.Build(_NUM_residues, _NUM_residues);
	_BOND_length.Build(_NUM_residues, _NUM_residues);
	_BOND_length2.Build(_NUM_residues, _NUM_residues);
	_BOND_delta.Build(_NUM_residues, _NUM_residues);
	_RMIN.Build(_NUM_residues, _NUM_residues);
	_RMAX.Build(_NUM_residues, _NUM_residues);
	_RMIN2.Build(_NUM_residues, _NUM_residues);
	_RMAX2.Build(_NUM_residues, _NUM_residues);
	_DDELTA.Build(_NUM_residues, _NUM_residues);
	_DDELTA_ea.Build(_NUM_residues, _NUM_residues);
	_PARA_KA.Build(_NUM_residues, _NUM_residues);
	_THETAZERO.Build(_NUM_residues, _NUM_residues);
	_THETAZERO_cos.Build(_NUM_residues, _NUM_residues);
	_THETAZERO_sin.Build(_NUM_residues, _NUM_residues);
	_PARA_KD.Build(_NUM_residues, _NUM_residues);
	_PHIZERO.Build(_NUM_residues, _NUM_residues);
	_PHIZERO_cos.Build(_NUM_residues, _NUM_residues);
	_PHIZERO_sin.Build(_NUM_residues, _NUM_residues);

	_EPSILON_eachres.Build(_NUM_residues, _NUM_residues);
	_LAMBDA_eachres.Build(_NUM_residues, _NUM_residues);
	_PPTYPE_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA2_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMADIS_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMAWELL_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMAWELL2_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA3_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA6_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA9_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA12_eachres.Build(_NUM_residues, _NUM_residues);
	_SIGMA24_eachres.Build(_NUM_residues, _NUM_residues);
	_E_cut_RR_eachres.Build(_NUM_residues, _NUM_residues);
	_R_cut_RR_eachres.Build(_NUM_residues, _NUM_residues);

	_Probability.Build(_NUM_replicas, _E_totalnum);
	_Probability_tot.Build(_NUM_replicas, _E_totalnum);
	_Probability_all=new double[_E_totalnum];

    ///////////////   code start  //////////////////
    tell_procid(); cout<<" memo_allocation done."<<endl;
	_Memory_freed=false;
	MPI_Barrier(MPI_COMM_WORLD);
}
////////////////////////////
void cmc::memo_setzero() {
	
	string TempStr;
	string TempInterval;
	int i=0;
	int j=0;
	tell_procid(); cout<<" memo_setzero start..."<<endl;
	///////////////   code start  //////////////////
	MEMOSETZERO(ec, sizeof(double)*3);
	//MEMOSETZERO(fc, sizeof(double)*3);
	MEMOSETZERO(dc, sizeof(double)*3);

	MEMOSETZERO(coorRES, sizeof(double)*3);
	MEMOSETZERO(vecDC, sizeof(double)*3);

	MEMOSETZERO(tc, sizeof(double)*3);

	MEMOSETZERO(ir, sizeof(int)*98);

	MEMOSETZERO(_T_rep_eachrep, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_T_eachrep, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_Index_T_rep, sizeof(int)*_NUM_replicas);
	MEMOSETZERO(_Index_rep_T, sizeof(int)*_NUM_replicas);
	MEMOSETZERO(_E_rep, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_NUM_rem, sizeof(int)*_NUM_replicas);
	MEMOSETZERO(_NUM_rem_acc, sizeof(int)*_NUM_replicas);

	MEMOSETZERO(_MC_NUM_SUC_stat, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_TOT_stat, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_SUC_all, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_TOT_all, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_FIL_all, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_FIL_stat, sizeof(double)*_NUM_replicas);

	MEMOSETZERO(_movelength_ea, sizeof(double)*_NUM_replicas);
	//MEMOSETZERO(_movelength_pv, sizeof(double)*_NUM_replicas);

	_ENER_LJ_eachatom.SetZero();
	MEMOSETZERO(_ENER_LJ_chosenatom_backup, sizeof(double)*_SIZE_memo);

	MEMOSETZERO(_ENER_BF_eachatom, sizeof(double)*_SIZE_memo);//i for bond i<->i+1
	MEMOSETZERO(_ENER_AG_eachatom, sizeof(double)*(_SIZE_memo+1));
	MEMOSETZERO(_ENER_DH_eachatom, sizeof(double)*(_SIZE_memo+2));
	//MEMOSETZERO(_Series_eachain, sizeof(int*)*_NUM_chains);

	_DIS2_eachatom.SetZero();
	MEMOSETZERO(_DIS2_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_x_eachatom.SetZero();
	MEMOSETZERO(_DIS_x_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_y_eachatom.SetZero();
	MEMOSETZERO(_DIS_y_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_z_eachatom.SetZero();
	MEMOSETZERO(_DIS_z_chosenatom_backup, sizeof(double)*_SIZE_memo);

	MEMOSETZERO(_COM_x, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_COM_y, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_COM_z, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_x, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_y, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_z, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_ec, sizeof(double)*_NUM_chains);
	_RG2_actual_x.SetZero();
	_RG2_actual_xtot.SetZero();
	_RG2_actual_y.SetZero();
	_RG2_actual_ytot.SetZero();
	_RG2_actual_z.SetZero();
	_RG2_actual_ztot.SetZero();

	if( _NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) {
		MEMOSETZERO(_DISSTAT, sizeof(double)*(_NUM_chains-1)*_NUM_chains/2);
		_DIS_stat.SetZero();
		_DIS_stat_tot.SetZero();
		_DIS_stat_actual.SetZero();
		_DIS_stat_actual_tot.SetZero();
	}
	for(i=0; i<_NUM_chains; i++){
		for(j=0; j<_NUM_chains; j++){
			_CINDEXMAP[i][j]=0.0;
		}
	}
	MEMOSETZERO(stat_com_x, sizeof(double)*_NUM_chains);
	MEMOSETZERO(stat_com_y, sizeof(double)*_NUM_chains);
	MEMOSETZERO(stat_com_z, sizeof(double)*_NUM_chains);

	_ENERLJ_stat_tot.SetZero();
	_ENERBF_stat_tot.SetZero();
	_ENERAG_stat_tot.SetZero();
	_ENERDH_stat_tot.SetZero();
	_COM_x_stat_tot.SetZero();
	_COM_y_stat_tot.SetZero();
	_COM_z_stat_tot.SetZero();
	/*_DIS_x_stat_tot.SetZero();
	_DIS_y_stat_tot.SetZero();
	_DIS_z_stat_tot.SetZero();*/

	/*_RG2_x_stat_tot.SetZero();
	_RG2_y_stat_tot.SetZero();
	_RG2_z_stat_tot.SetZero();*/
	_RG2_stat_tot.SetZero();

	if(_len_EINT>0) {
		MEMOSETZERO(_CN_stat, sizeof(double)*_E_totalnum);
		MEMOSETZERO(_CN_stat_tot, sizeof(double)*_E_totalnum);
		if(_EINT_totalnum>1) {
			_EINT_stat.SetZero();
			_EINT_stat_tot.SetZero();
		}
	}
		
	MEMOSETZERO(_ENER_delta_LJ_inter, sizeof(double)*_NUM_chains);
	MEMOSETZERO(tempElj, sizeof(double)*(_NUM_chains+1)*_NUM_chains/2);
	MEMOSETZERO(tempEag, sizeof(double)*_NUM_chains);
	MEMOSETZERO(tempEbf, sizeof(double)*_NUM_chains);
	MEMOSETZERO(tempEdh, sizeof(double)*_NUM_chains);
	_COM_x_stat.SetZero();
	_COM_y_stat.SetZero();
	_COM_z_stat.SetZero();
	/*_DIS_x_stat.SetZero();
	_DIS_y_stat.SetZero();
	_DIS_z_stat.SetZero();*/
	
	/*_RG2_x_stat.SetZero();
	_RG2_y_stat.SetZero();
	_RG2_z_stat.SetZero();*/
	_RG2_stat.SetZero();
	_ENERLJ_stat.SetZero();
	_ENERBF_stat.SetZero();
	_ENERAG_stat.SetZero();
	_ENERDH_stat.SetZero();
	//MEMOSETZERO(_rgyration2, sizeof(double)*_NUM_chains);
	//MEMOSETZERO(_rhead2, sizeof(double)*_NUM_chains);

	MEMOSETZERO(_INDEX_chn_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_SIZE_of_chn, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_SIZE_of_chn_real, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_HEAD, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_TAIL, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_HEAD_real, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_TAIL_real, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_atm_in_chn, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_atm_in_chn_real, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_TYPE_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_TYPE_atom_real, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_LNEIGHBOR, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_RNEIGHBOR, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_LNEIGHBOR_real, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_RNEIGHBOR_real, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_RES_ATM, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_XX, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_YY, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_ZZ, sizeof(double)*_SIZE_memo);
	_XX_allrep.SetZero();
	_YY_allrep.SetZero();
	_ZZ_allrep.SetZero();
	MEMOSETZERO(_XX_rec, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_YY_rec, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_ZZ_rec, sizeof(double)*_SIZE_memo);

	MEMOSETZERO(_INDEX_chn_or_not_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_chn_or_not_chain, sizeof(int)*_NUM_chains);

	_PARA_KB.SetZero();
	_BOND_length.SetZero();
	_BOND_length2.SetZero();
	_BOND_delta.SetZero();
	_RMIN.SetZero();
	_RMAX.SetZero();
	_RMIN2.SetZero();
	_RMAX2.SetZero();
	_DDELTA.SetZero();
	_DDELTA_ea.SetZero();
	_PARA_KA.SetZero();
	_THETAZERO.SetZero();
	_THETAZERO_cos.SetZero();
	_THETAZERO_sin.SetZero();
	_PARA_KD.SetZero();
	_PHIZERO.SetZero();
	_PHIZERO_cos.SetZero();
	_PHIZERO_sin.SetZero();

	_EPSILON_eachres.SetZero();
	_LAMBDA_eachres.SetZero();
	_PPTYPE_eachres.SetZero();
	_SIGMA_eachres.SetZero();
	_SIGMA2_eachres.SetZero();
	_SIGMADIS_eachres.SetZero();
	_SIGMAWELL_eachres.SetZero();
	_SIGMAWELL2_eachres.SetZero();
	_SIGMA3_eachres.SetZero();
	_SIGMA6_eachres.SetZero();
	_SIGMA9_eachres.SetZero();
	_SIGMA12_eachres.SetZero();
	_SIGMA24_eachres.SetZero();
	_E_cut_RR_eachres.SetZero();
	_R_cut_RR_eachres.SetZero();

	_Probability.SetZero();
	_Probability_tot.SetZero();

	MEMOSETZERO(_Probability_all, sizeof(double)*_E_totalnum);

    ///////////////   code start  //////////////////
    tell_procid(); cout<<" memo_setzero done."<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
}
///////////////////////////
void cmc::memo_evaluation_fnode() {
	if( _PROC_ID!=0 ) {
		return;
	}

	tell_procid(); cout<<" memo_evaluation_fnode start..."<<endl;

	int Size_RES=0;
	int Size_ATM=0;
	int index_atm=0;
	int index_res=0;

	int i=0;
	int j=0;
	int k=0;
	int l=0;
	int numerator=0;
	/*for(i=1; i<=_NUM_atoms; i++) {
		_INDEX_chn_atom[i]=_system_._Index_CHN[i];
		_INDEX_atm_in_chn[i]=_system_._Index_ATM_in_CHN[i];
		_TYPE_atom[i]=_system_._Type_ATOM[i];
		_XX[i]=_system_._XX[i];
		_YY[i]=_system_._YY[i];
		_ZZ[i]=_system_._ZZ[i];
	}*/
	int _CHCK_UPPER_POINT=20;
	_NUM_chains_real=1;
	for(i=0; i<_NUM_atoms-1; i++) {
		if(_system_.allatoms[i].atomname_special!=_system_.allatoms[i+1].atomname_special) {
			_NUM_chains_real+=1;
		}
	}
	int index_atm_real=0;
	cout<<" there are actually "<<_NUM_chains_real<<" chains moving around in the system."<<endl;
	for(i=0; i<_NUM_chains; i++) {
		//cout<<i<<endl;
		Size_RES=_system_.chains[i].nresidues;
		if( _system_.chains[i].index_chn_real<_CHCK_UPPER_POINT ) {
			_INDEX_chn_or_not_chain[i]=true;
		} else {
			_INDEX_chn_or_not_chain[i]=false;
		}
		for(l=numerator-1;l>=0;l--) {
			if(_system_.allatoms[numerator].atomname_special!=_system_.allatoms[l].atomname_special) {
				break;
			}
		}
		_INDEX_CHN_HEAD_real[i]=l+2;
		for(l=numerator+1;l<_NUM_atoms;l++) {
			if(_system_.allatoms[numerator].atomname_special!=_system_.allatoms[l].atomname_special) {
				break;
			}
		}
		_INDEX_CHN_TAIL_real[i]=l;
		_SIZE_of_chn_real[i]=_INDEX_CHN_TAIL_real[i]-_INDEX_CHN_HEAD_real[i]+1;
		cout<<" chn["<<i<<"]s("<<_SIZE_of_chn_real[i]<<")<->["<<_INDEX_CHN_HEAD_real[i]<<":"<<_INDEX_CHN_TAIL_real[i]<<"]"<<endl;
		index_atm=0;
		for(j=0; j<Size_RES; j++) {
			//cout<<j<<endl;
			Size_ATM=_system_.chains[i].residues[j].natoms;
			//cout<<Size_ATM<<" "<<Size_RES<<endl;
			for(k=0; k<Size_ATM; k++) {
				//cout<<"k="<<k<<endl;
				/*_XX[numerator+1]=_system_.chains[i].residues[j].atoms[k].x_coordinate;
				_YY[numerator+1]=_system_.chains[i].residues[j].atoms[k].y_coordinate;
				_ZZ[numerator+1]=_system_.chains[i].residues[j].atoms[k].z_coordinate;*/
				_XX[numerator+1]=_system_.allatoms[numerator].x_coordinate;
				_YY[numerator+1]=_system_.allatoms[numerator].y_coordinate;
				_ZZ[numerator+1]=_system_.allatoms[numerator].z_coordinate;
				//cout<<_XX[numerator+1]<<" "<<_YY[numerator+1]<<" "<<_ZZ[numerator+1]<<endl; 
				if( j==0 && k==0 ) { // head of this chain;
					_TYPE_atom[numerator+1]=-1; //head
					_INDEX_CHN_HEAD[i]=numerator+1;
					if(numerator==0) {//1st must be 1st;
						_TYPE_atom_real[numerator+1]=-1;
						index_atm_real=0;
					} else {//not 1st of the molecule;
						if(_system_.allatoms[numerator].atomname_special==_system_.allatoms[numerator-1].atomname_special) {
							_TYPE_atom_real[numerator+1]=0; //not head
							/*for(l=numerator-1;l>=0;l--) {
								if(_system_.allatoms[numerator].atomname_special==_system_.allatoms[l].atomname_special) {
									break;
								}
							}*/
						} else {
							_TYPE_atom_real[numerator+1]=-1; //head
							index_atm_real=0;
						}
					}
				} else { // not head of this chain;
					_TYPE_atom[numerator+1]=0;
					_TYPE_atom_real[numerator+1]=0;
				} 
				if( j==(Size_RES-1) && k==(Size_ATM-1) ) {
					_TYPE_atom[numerator+1]=1; //tail
					_INDEX_CHN_TAIL[i]=numerator+1;
					if(numerator==_NUM_atoms-1) {//last must be last;
						_TYPE_atom_real[numerator+1]=1;
					} else {//not 1st of the molecule;
						if(_system_.allatoms[numerator].atomname_special==_system_.allatoms[numerator+1].atomname_special) {
							if(numerator!=0) {
								_TYPE_atom_real[numerator+1]=0; //not tail
							} else {
								_TYPE_atom_real[numerator+1]=-1; //head
								index_atm_real=0;
							}
						} else {
							_TYPE_atom_real[numerator+1]=1; //tail
						}
					}
				}
				_INDEX_chn_atom[numerator+1]=i;
				_INDEX_chn_or_not_atom[numerator+1]=_INDEX_chn_or_not_chain[i];
				_INDEX_RES_ATM[numerator+1]=_system_.indexatmres[numerator]; //right
				//_Index_ATM_in_RES[numerator+1]=k;
				_INDEX_atm_in_chn[numerator+1]=index_atm+k; //right
				_INDEX_atm_in_chn_real[numerator+1]=index_atm_real; //right
				index_atm_real++;
				numerator++;
			    //cout<<_Index_CHN_real[numerator+1]<<endl;
			}
			index_atm+=Size_ATM;
		}
		_SIZE_of_chn[i]=index_atm;
		
		index_res+=Size_RES;
	}
	/*numerator=0;
	for(i=0; i<_NUM_chains; i++) {
		_Series_eachain[i]=new int[_SIZE_of_chn[i]];
		for(j=0; j<_NUM_atoms; j++) {
			if( _system_.atoms[numerator].atomindex==_system_.chains[i].atoms )
		if(numerator==0) {
			_TYPE_atom[numerator+1]=-1; //head
			_INDEX_CHN_HEAD[i]=numerator+1;
		}
		numerator++;		
	}*/
	//_sigma_surf_1_6=pow(0.4, 1.0/6.0)*_sigma_surf; //new_edition for 9, 3 LJ of semi-finite substrate;
	tell_procid(); cout<<" memo_evaluation_fnode done."<<endl;
}
///////////////////////////
void cmc::memo_evaluation_bcast() {
	MPI_Bcast(_XX, _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_YY, _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_ZZ, _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_TYPE_atom, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_TYPE_atom_real, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_CHN_HEAD, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_CHN_TAIL, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_CHN_HEAD_real, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_CHN_TAIL_real, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_chn_atom, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_chn_or_not_atom, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_RES_ATM, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_atm_in_chn, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_atm_in_chn_real, _SIZE_memo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_INDEX_chn_or_not_chain, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIZE_of_chn, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIZE_of_chn_real, _NUM_chains, MPI_INT, 0, MPI_COMM_WORLD);
}
///////////////////////////
inline void cmc::loadrestart(double ogboxlx, double ogboxly, double ogboxlz) {
	if(_RESTARTFlag) {
		cout<<" @Proc"<<setw(3)<<_PROC_ID<<" :: restart being loading."<<endl;
		char* tFILENAME_conf;
		tFILENAME_conf=new char[10];
		chck_bond_len();//after each make mapping!
		sprintf(tFILENAME_conf, "confT%03d.pdb", _PROC_ID+1);
		_system_.memo_free();
		cout<<" @Proc"<<setw(3)<<_PROC_ID<<" loading [ "<<tFILENAME_conf<<" ] ... "<<endl;
		_system_.readpdbinfo(tFILENAME_conf, 1.0, false, true);
		int i=0;
		for(i=0; i<_NUM_atoms; i++) {
			_XX[i+1]=_system_.allatoms[i].x_coordinate;
			_YY[i+1]=_system_.allatoms[i].y_coordinate;
			_ZZ[i+1]=_system_.allatoms[i].z_coordinate;
			//cout<<" "<<_XX[i]<<" "<<_YY[i]
		}
		cout<<" @Proc"<<setw(3)<<_PROC_ID<<" rebuilding the system..."<<endl;
		for(i=1; i<_NUM_atoms; i++) {
			while( (_XX[i+1]-_XX[i])>(ogboxlx/2.0) ) {
				_XX[i+1]=_XX[i+1]-ogboxlx;
			}
			while( (_XX[i+1]-_XX[i])<(-ogboxlx/2.0) ) {
				_XX[i+1]=_XX[i+1]+ogboxlx;
			}
			while( (_YY[i+1]-_YY[i])>(ogboxly/2.0) ) {
				_YY[i+1]=_YY[i+1]-ogboxly;
			}
			while( (_YY[i+1]-_YY[i])<(-ogboxly/2.0) ) {
				_YY[i+1]=_YY[i+1]+ogboxly;
			}
			while( (_ZZ[i+1]-_ZZ[i])>(ogboxlz/2.0) ) {
				_ZZ[i+1]=_ZZ[i+1]-ogboxlz;
			}
			while( (_ZZ[i+1]-_ZZ[i])<(-ogboxlz/2.0) ) {
				_ZZ[i+1]=_ZZ[i+1]+ogboxlz;
			}
		}
		/*cout<<" @Proc"<<setw(3)<<_PROC_ID<<" rebuilt."<<endl;
		for(i=0; i<_NUM_atoms; i++) {
			cout<<" id"<<setw(3)<<i+1<<" "<<setw(10)<<_XX[i+1]<<" "<<setw(10)<<_YY[i+1]<<" "<<setw(10)<<_ZZ[i+1]<<endl;
			//cout<<" "<<_XX[i]<<" "<<_YY[i]
		}*/
	}
}
///////////////////////////
void cmc::memo_evaluation() {
	init_epsilonsigma(); //can not change the order; done
	broadcast_epsilonsigma();
	load_temperatures();
	scatter_temperatures();
	memo_evaluation_fnode();
	memo_evaluation_bcast();
	loadrestart(_ogboxlx, _ogboxly, _ogboxlz);

	int i=0;
	int j=0;
	int num=0;
	for(i=0; i<_NUM_chains; i++){
		for(j=0; j<i; j++){
			_CINDEXMAP[i][j]=num++;
			_CINDEXMAP[j][i]=_CINDEXMAP[i][j];
		}
	}
	for(i=0; i<_NUM_chains; i++){
		_CINDEXMAP[i][i]=num++;
	}
	if(_PROC_ID==0) {
		for(i=0; i<_NUM_chains; i++){
			for(j=0; j<i; j++){
				if(i!=j) {
					cout<<" "<<i<<":"<<j<<":"<<_CINDEXMAP[i][j];
				}
			}
			cout<<endl;
		}
		for(i=0; i<_NUM_chains; i++){
			cout<<" "<<i<<":"<<i<<":"<<_CINDEXMAP[i][i];
		}
		cout<<endl;
	}

	if(_len_EINT>0) {
		_EINT_index=_CINDEXMAP[_EINTlist[0]][_EINTlist[1]];
		tell_procid(); cout<<" eint index: "<<_EINTlist[0]<<":"<<_EINTlist[1]<<":"<<_EINT_index<<endl;
		_CN_backup.Build(_SIZE_of_chn[_EINTlist[0]], _SIZE_of_chn[_EINTlist[1]]);
		_CN_backup.SetZero();
		/*for(int i=0; i<_SIZE_of_chn[_EINTlist[0]]; i++) {
			for(int j=0; j<_SIZE_of_chn[_EINTlist[1]]; j++) {
				cout<<" "<<_CN_backup.pArray[i][j];
			}
			cout<<endl;
		}*/
	}

	_RecordingNAME=new char[50];
	
	for(i=1;i<_SIZE_memo;i++) {
		_INDEX_LNEIGHBOR[i]=(i-_NEIGHBOR_lj)<=_INDEX_CHN_HEAD[_INDEX_chn_atom[i]]?_INDEX_CHN_HEAD[_INDEX_chn_atom[i]]:(i-_NEIGHBOR_lj); //[
		_INDEX_RNEIGHBOR[i]=(i+_NEIGHBOR_lj)>=_INDEX_CHN_TAIL[_INDEX_chn_atom[i]]?_INDEX_CHN_TAIL[_INDEX_chn_atom[i]]:(i+_NEIGHBOR_lj); //]
	}
	for(i=1;i<_SIZE_memo;i++) {
		_INDEX_LNEIGHBOR_real[i]=(i-_NEIGHBOR_lj)<=_INDEX_CHN_HEAD_real[_INDEX_chn_atom[i]]?_INDEX_CHN_HEAD_real[_INDEX_chn_atom[i]]:(i-_NEIGHBOR_lj); //[
		_INDEX_RNEIGHBOR_real[i]=(i+_NEIGHBOR_lj)>=_INDEX_CHN_TAIL_real[_INDEX_chn_atom[i]]?_INDEX_CHN_TAIL_real[_INDEX_chn_atom[i]]:(i+_NEIGHBOR_lj); //]
	}
	for(i=1;i<_SIZE_memo-1;i++) { // important, you have to do this before calculation of eachatom: <angle calculation>
		for(j=i+1;j<_SIZE_memo;j++) {
			_DIS_x_eachatom.pArray[i][j]=DIS_PBC_X(_XX[i],_XX[j]);
			_DIS_y_eachatom.pArray[i][j]=DIS_PBC_Y(_YY[i],_YY[j]);
			_DIS_z_eachatom.pArray[i][j]=DIS_PBC_Z(_ZZ[i],_ZZ[j]);
			_DIS2_eachatom.pArray[i][j]=_DIS_x_eachatom.pArray[i][j]*_DIS_x_eachatom.pArray[i][j]
					                   +_DIS_y_eachatom.pArray[i][j]*_DIS_y_eachatom.pArray[i][j]
					                   +_DIS_z_eachatom.pArray[i][j]*_DIS_z_eachatom.pArray[i][j];	
			/*if( (i==6) && (j==7) ) {
				cout<<" x"<<j<<"-x"<<i<<"="<<_XX[j]-_XX[i]<<endl
				    <<" y"<<j<<"-y"<<i<<"="<<_YY[j]-_YY[i]<<endl
				    <<" z"<<j<<"-z"<<i<<"="<<_ZZ[j]-_ZZ[i]<<endl;
				cout<<" x"<<j<<"-x"<<i<<"="<<_DIS_x_eachatom.pArray[i][j]<<endl
				    <<" y"<<j<<"-y"<<i<<"="<<_DIS_y_eachatom.pArray[i][j]<<endl
				    <<" z"<<j<<"-z"<<i<<"="<<_DIS_z_eachatom.pArray[i][j]<<endl;
			}*/
			_DIS_x_eachatom.pArray[j][i]=-_DIS_x_eachatom.pArray[i][j];
			_DIS_y_eachatom.pArray[j][i]=-_DIS_y_eachatom.pArray[i][j];
			_DIS_z_eachatom.pArray[j][i]=-_DIS_z_eachatom.pArray[i][j];
			_DIS2_eachatom.pArray[j][i]=_DIS2_eachatom.pArray[i][j];
		}
	}
	/*for(i=0; i<_NUM_chains; i++) {
		if(_SIZE_of_chn[i]<3) {
			_RG2_x_stat.SetZero();
			RG
		}
	}*/

	chck_bond_len();
	if( _NUM_atoms==1 && _INDEX_chn_or_not_atom[1] ) {
		cout<<" _NUM_atoms="<<_NUM_atoms<<endl;
		cout<<" _INDEX_chn_or_not_atom[1]="<<_INDEX_chn_or_not_atom[1]<<endl;
		cout<<" chain with only one atom can not move ... "<<endl;
		cout<<" try again and make sure if this atom belongs to a chain or not (idx_chn > 20 ?)"<<endl;
		//exit(-1);
		exit(LOGICERROR);
	}
	tell_procid(); cout<<" memo_evaluation done."<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	check_memo();
	///////////////   code start  //////////////////
	///////////////   code start  //////////////////
}
//////////////////////////
inline void cmc::check_memo() {
	//tell_procid();
	tell_procid(); cout<<" check the molecule information .... "<<endl;
	int i;
	for(i=0; i<_NUM_chains; i++) {
		if(_INDEX_chn_or_not_chain[i]) {
			tell_procid(); cout<<" the szchn["<<i<<"] is: "<<_SIZE_of_chn[i]<<"r("<<_SIZE_of_chn_real[i]<<") ("<<_INDEX_CHN_HEAD[i]<<","<<_INDEX_CHN_TAIL[i]<<")"<<"("
				<<_INDEX_CHN_HEAD_real[i]<<","<<_INDEX_CHN_TAIL_real[i]<<")"<<endl;
		} else {
			tell_procid(); cout<<" the szncn["<<i<<"] is: "<<_SIZE_of_chn[i]<<"r("<<_SIZE_of_chn_real[i]<<") ("<<_INDEX_CHN_HEAD[i]<<","<<_INDEX_CHN_TAIL[i]<<")"<<"("
				<<_INDEX_CHN_HEAD_real[i]<<","<<_INDEX_CHN_TAIL_real[i]<<")"<<endl;
		}
	}
	if(!_IFVERBOSE) { 
		return;
	}
	int tnum_chn=int( log10(double(_NUM_chains)));
	int tidx_atm=int(log10(_NUM_atoms/10))+1;
	//if(_PROC_ID==0){ return;}
	for(i=1; i<_SIZE_memo; i++) {
		//tidx_atm=int( log10(double(_SIZE_of_chn[_INDEX_chn_atom[i]]) ))+2;
		if(_INDEX_chn_or_not_atom[i]) {
			cout<<" ccatom [";
		} else {
			cout<<" ncatom [";
		}
		cout<<setw(tidx_atm)<<i<<"] on chn:["<<setw(tnum_chn)<<_INDEX_chn_atom[i]<<"] the ["
		    <<setw(tidx_atm)<<_INDEX_atm_in_chn[i]+1<<"]inchn ["
		    <<setw(tidx_atm)<<_INDEX_atm_in_chn_real[i]+1<<"]inrcn ["
		    <<setw(tidx_atm)<<_INDEX_LNEIGHBOR[i]<<","<<setw(tidx_atm)<<_INDEX_RNEIGHBOR[i]<<"]th r["
		    <<setw(tidx_atm)<<_INDEX_LNEIGHBOR_real[i]<<","<<setw(tidx_atm)<<_INDEX_RNEIGHBOR_real[i]<<"], tp["
		    <<setw(2)<<_TYPE_atom[i]<<"], rtp["
		    <<setw(2)<<_TYPE_atom_real[i]
		    <<"] idx_r_in_m: ["<<_INDEX_RES_ATM[i]
		    <<"] c:["<<_XX[i]<<","<<_YY[i]<<","<<_ZZ[i]<<"];"<<endl;
	}
	cout<<" check the molecule information done! "<<endl<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
}
///////////////////////////
void cmc::memo_free() {	
	tell_procid(); cout<<" memory deallocation ... "<<endl;
	if(_Memory_freed) { 
		tell_procid(); cout<<" memory freed already~ "<<endl;
		return; 
	}
	delete[] _RecordingNAME;
	delete[] ec;
	ec=NULL;
	//delete[] fc;
	delete[] dc;
	dc=NULL;
	delete[] coorRES;
	coorRES=NULL;
	delete[] vecDC;
	vecDC=NULL;
	delete[] tc;
	tc=NULL;
	delete[] ir;
	ir=NULL;

	delete[] _T_rep_eachrep; _T_rep_eachrep=NULL;
	delete[] _T_eachrep; _T_eachrep=NULL;
	delete[] _Index_T_rep; _Index_T_rep=NULL;
	delete[] _Index_rep_T; _Index_rep_T=NULL;
	delete[] _E_rep; _E_rep=NULL;
	delete[] _NUM_rem; _NUM_rem=NULL;
	delete[] _NUM_rem_acc; _NUM_rem_acc=NULL;

	delete[] _MC_NUM_SUC_stat; 
	delete[] _MC_NUM_TOT_stat; 
	delete[] _MC_NUM_SUC_all; 
	delete[] _MC_NUM_TOT_all; 

	delete[] _MC_NUM_FIL_all; 
	delete[] _MC_NUM_FIL_stat; 

	delete[] _movelength_ea; 
	//delete[] _movelength_pv; 

	_ENER_LJ_eachatom.Release();
	delete[] _ENER_LJ_chosenatom_backup;
	_ENER_LJ_chosenatom_backup=NULL;

	delete[] _ENER_BF_eachatom;
	_ENER_BF_eachatom=NULL;
	delete[] _ENER_AG_eachatom;
	_ENER_AG_eachatom=NULL;
	delete[] _ENER_DH_eachatom;
	_ENER_DH_eachatom=NULL;

	//delete[] _Series_eachain;
	//_Series_eachain=NULL;

	_DIS2_eachatom.Release();
	delete[] _DIS2_chosenatom_backup;
	_DIS2_chosenatom_backup=NULL;
	_DIS_x_eachatom.Release();
	delete[] _DIS_x_chosenatom_backup;
	_DIS_x_chosenatom_backup=NULL;
	_DIS_y_eachatom.Release();
	delete[] _DIS_y_chosenatom_backup;
	_DIS_y_chosenatom_backup=NULL;
	_DIS_z_eachatom.Release();
	delete[] _DIS_z_chosenatom_backup;
	_DIS_z_chosenatom_backup=NULL;
	if(_len_OP!=0) {
		delete[] _OPlist;
		_OPlist=NULL;
	}
	if(_len_EINT>0 ) {
		_CN_backup.Release();
		delete[] _CN_stat;
		delete[] _CN_stat_tot;
		delete[] _EINTlist;
		_EINTlist=NULL;
		if(_EINT_totalnum>1) {
			_EINT_stat.Release();
			_EINT_stat_tot.Release();
		}
	}
	delete[] _COM_x;
	delete[] _COM_y;
	delete[] _COM_z;
	_COM_x=NULL;
	_COM_y=NULL;
	_COM_z=NULL;
	delete[] _RG2_x;
	delete[] _RG2_y;
	delete[] _RG2_z;
	delete[] _RG2_ec;
	_RG2_x=NULL;
	_RG2_y=NULL;
	_RG2_z=NULL;
	_RG2_ec=NULL;
	if( _NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) {
		delete[] _DISSTAT; _DISSTAT=NULL;
		_DIS_stat.Release();
		_DIS_stat_tot.Release();
		_DIS_stat_actual.Release();
		_DIS_stat_actual_tot.Release();
	}

	delete[] stat_com_x;
	stat_com_x=NULL;
	delete[] stat_com_y;
	stat_com_y=NULL;
	delete[] stat_com_z;
	stat_com_z=NULL;
	
	if(_ACCINDEX==0) { 
		//_EBF_ACC.Release();
		_ELJ_ACC.Release();
		//_BF_interval.Release();
		_LJ_interval.Release();
	}
	for(int i=0; i<_NUM_chains; i++){
		delete[] _CINDEXMAP[i];
	}
	delete[] _CINDEXMAP;
	delete[] _indexRGSTAT;

	delete[] _ENER_delta_LJ_inter;
	delete[] tempElj;
	delete[] tempEag;
	delete[] tempEbf;
	delete[] tempEdh;
	_RG2_actual_x.Release();
	_RG2_actual_xtot.Release();
	_RG2_actual_y.Release();
	_RG2_actual_ytot.Release();
	_RG2_actual_z.Release();
	_RG2_actual_ztot.Release();
	_COM_x_stat.Release();
	_COM_y_stat.Release();
	_COM_z_stat.Release();
	/*_DIS_x_stat.Release();
	_DIS_y_stat.Release();
	_DIS_z_stat.Release();*/
	
	/*_RG2_x_stat.Release();
	_RG2_y_stat.Release();
	_RG2_z_stat.Release();*/
	_RG2_stat.Release();
	_ENERLJ_stat.Release();
	_ENERBF_stat.Release();
	_ENERAG_stat.Release();
	_ENERDH_stat.Release();


	_ENERLJ_stat_tot.Release();
	_ENERAG_stat_tot.Release();
	_ENERBF_stat_tot.Release();
	_ENERDH_stat_tot.Release();
	_COM_x_stat_tot.Release();
	_COM_y_stat_tot.Release();
	_COM_z_stat_tot.Release();
	/*_DIS_x_stat_tot.Release();
	_DIS_y_stat_tot.Release();
	_DIS_z_stat_tot.Release();*/

	/*_RG2_x_stat_tot.Release();
	_RG2_y_stat_tot.Release();
	_RG2_z_stat_tot.Release(); */ 
	_RG2_stat_tot.Release();  

	//delete[] _rgyration2;
	//_rgyration2=NULL;
	//delete[] _rhead2;
	//_rhead2=NULL;

	delete[] _INDEX_chn_atom;
	_INDEX_chn_atom=NULL;
	delete[] _SIZE_of_chn;
	_SIZE_of_chn=NULL;
	delete[] _SIZE_of_chn_real;
	_SIZE_of_chn_real=NULL;
	delete[] _INDEX_CHN_HEAD;
	_INDEX_CHN_HEAD=NULL;
	delete[] _INDEX_CHN_TAIL;
	_INDEX_CHN_TAIL=NULL;
	delete[] _INDEX_CHN_HEAD_real;
	_INDEX_CHN_HEAD_real=NULL;
	delete[] _INDEX_CHN_TAIL_real;
	_INDEX_CHN_TAIL_real=NULL;
	
	delete[] _INDEX_atm_in_chn;
	_INDEX_atm_in_chn=NULL;  
	delete[] _INDEX_atm_in_chn_real;
	_INDEX_atm_in_chn_real=NULL;  
	
	delete[] _TYPE_atom;
	_TYPE_atom=NULL;
	delete[] _TYPE_atom_real;
	_TYPE_atom_real=NULL;
	
	
	delete[] _INDEX_LNEIGHBOR;
	_INDEX_LNEIGHBOR=NULL;
	delete[] _INDEX_RNEIGHBOR;
	_INDEX_RNEIGHBOR=NULL;
	delete[] _INDEX_LNEIGHBOR_real;
	_INDEX_LNEIGHBOR_real=NULL;
	delete[] _INDEX_RNEIGHBOR_real;
	_INDEX_RNEIGHBOR_real=NULL;
	delete[] _INDEX_RES_ATM;
	_INDEX_RES_ATM=NULL;
    delete[] _XX;
	delete[] _YY;
	delete[] _ZZ;
    _XX=NULL;
    _YY=NULL;
    _ZZ=NULL;
    _XX_allrep.Release();
    _YY_allrep.Release();
    _ZZ_allrep.Release();
	delete[] _XX_rec;
	delete[] _YY_rec;
	delete[] _ZZ_rec;

	delete[] _INDEX_chn_or_not_atom;
	_INDEX_chn_or_not_atom=NULL;
	delete[] _INDEX_chn_or_not_chain;
	_INDEX_chn_or_not_chain=NULL;
	

	_PARA_KB.Release();
	_BOND_length.Release();
	_BOND_length2.Release();
	_BOND_delta.Release();
	_RMIN.Release();
	_RMAX.Release();
	_RMIN2.Release();
	_RMAX2.Release();
	_DDELTA.Release();
	_DDELTA_ea.Release();
	_PARA_KA.Release();
	_THETAZERO.Release();
	_THETAZERO_cos.Release();
	_THETAZERO_sin.Release();
	_PARA_KD.Release();
	_PHIZERO.Release();
	_PHIZERO_cos.Release();
	_PHIZERO_sin.Release();

	_EPSILON_eachres.Release();
	_LAMBDA_eachres.Release();
	_PPTYPE_eachres.Release();
	_SIGMA_eachres.Release();
	_SIGMA2_eachres.Release();
	_SIGMADIS_eachres.Release();
	_SIGMAWELL_eachres.Release();
	_SIGMAWELL2_eachres.Release();
	_SIGMA3_eachres.Release();
	_SIGMA6_eachres.Release();
	_SIGMA9_eachres.Release();
	_SIGMA12_eachres.Release();
	_SIGMA24_eachres.Release();
	_E_cut_RR_eachres.Release();
	_R_cut_RR_eachres.Release();

	_Probability.Release();
	_Probability_tot.Release();
	delete[] _Probability_all;
	//cout<<"mamamam11"<<endl;

	if(!_Statistic_over) { close_statistic(); }
    //if( _PROC_ID==0 ) { tell_procid(); cout<<" _ZZm freed successfully ! "<<endl; }
	tell_procid(); cout<<" memory deallocation done!"<<endl;
	_Memory_freed=true;
	//cout<<"mamamam11"<<endl;
}
///////////////////////////
void cmc::load_parameters(const string FILENAME_para) { //only for fnode
	if(_PROC_ID!=0) {
		return;
	}
	ifstream para_ifstream(FILENAME_para.c_str());
	if(para_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
		exit(IOERROR);
	}
	cout<<endl<<" Loading parameters from file [ "<<FILENAME_para<<" ] ..."<<endl;
	string tempstr;
	vector<string> tempvec;
	while( getline(para_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		if(tempvec.size()!=0) {
			readparameter(tempstr, string("_IFVERBOSE"), _IFVERBOSE);
			readparameter(tempstr, string("_FILENAME_conf"), _FILENAME_conf);
			
			readparameter(tempstr, string("_CoorX1"), _CoorX1);
			readparameter(tempstr, string("_CoorX2"), _CoorX2);
			readparameter(tempstr, string("_CoorY1"), _CoorY1);
			readparameter(tempstr, string("_CoorY2"), _CoorY2);
			readparameter(tempstr, string("_CoorZ1"), _CoorZ1);
			readparameter(tempstr, string("_CoorZ2"), _CoorZ2);
			
			readparameter(tempstr, string("_PBC_Dim"), _PBC_Dim);
			readparameter(tempstr, string("_BF_flag"), _BF_flag);
			readparameter(tempstr, string("_BF_coeff"), _BF_coeff);
			readparameter(tempstr, string("_BF_scale"), _BF_scale);

			readparameter(tempstr, string("_AG_flag"), _AG_flag);
			readparameter(tempstr, string("_DH_flag"), _DH_flag);
			readparameter(tempstr, string("_FFPARA_FN"), _FFPARA_FN);
			readparameter(tempstr, string("_NEIGHBOR_lj"), _NEIGHBOR_lj);
			readparameter(tempstr, string("_FLAG_ener"), _FLAG_ener);
			readparameter(tempstr, string("_FLAG_pivot"), _FLAG_pivot);
			readparameter(tempstr, string("_FLAG_rg2"), _FLAG_rg2);
			readparameter(tempstr, string("_FLAG_dis"), _FLAG_dis);

			readparameter(tempstr, string("_ACCINDEX"), _ACCINDEX);

			readparameter(tempstr, string("_NUM_replicas"), _NUM_replicas);
			readparameter(tempstr, string("_RUNTIMES_eachstep"), _RUNTIMES_eachstep);
			readparameter(tempstr, string("_RUNTIMES_totalnum"), _RUNTIMES_totalnum);
			readparameter(tempstr, string("_RUNTIMES_recording"), _RUNTIMES_recording);
			//readparameter(tempstr, string("_RESTARTFlag"), _RESTARTFlag);
			readparameter(tempstr, string("_RUNTIMES_output"), _RUNTIMES_output);
			readparameter(tempstr, string("_RUNTIMES_remgap"), _RUNTIMES_remgap);

			readparameter(tempstr, string("_E_lowest"), _E_lowest);
			readparameter(tempstr, string("_E_interval"), _E_interval);
			readparameter(tempstr, string("_E_totalnum"), _E_totalnum);
			readparameter(tempstr, string("_RG_lowest"), _RG_lowest);
			readparameter(tempstr, string("_RG_interval"), _RG_interval);
			readparameter(tempstr, string("_RG_totalnum"), _RG_totalnum);
			readparameter(tempstr, string("_EINT_lowest"), _EINT_lowest);
			readparameter(tempstr, string("_EINT_interval"), _EINT_interval);
			readparameter(tempstr, string("_EINT_totalnum"), _EINT_totalnum);
			//readparameter(tempstr, string("_MPI_OR_NOT"), _MPI_OR_NOT);
			if( tempvec[0] == string("_OPlist") ) {
				cout<<" oplist: ";
				_len_OP=tempvec.size()-1;
				if(_len_OP!=0) {
					_OPlist=new int[_len_OP];
				}
				for(int i=0; i<_len_OP; i++) {
					_OPlist[i]=atoi(tempvec[i+1].c_str());
					cout<<_OPlist[i]<<" ";
				}
				cout<<endl;
			}
			if( tempvec[0] == string("_EINTlist") ) {
				cout<<" _EINTlist: ";
				_len_EINT=tempvec.size()-1;
				if(_len_EINT>2) {
					cout<<"_EINTlist now supports only Eint between 2 chains;"<<endl;
					exit(LOGICERROR);
				}
				if(_len_EINT!=0) {
					_EINTlist=new int[_len_EINT];
				}
				for(int i=0; i<_len_EINT; i++) {
					_EINTlist[i]=atoi(tempvec[i+1].c_str())-1;
					cout<<_EINTlist[i]<<" ";
				}
				cout<<endl;
			}
			if( tempvec[0] == string("_RESTARTFlag") ) {
				_ogboxlx=0.0;
				_ogboxly=0.0;
				_ogboxlz=0.0;
				if( tempvec.size()<2 ) {
					cout<<" you have to specify the value of _RESTARTFlag; exit() load_para;"<<endl;
					exit(IOERROR);
				} 
				readparameter(tempstr, string("_RESTARTFlag"), _RESTARTFlag);
				if(_RESTARTFlag>0 && tempvec.size()<5) {
					cout<<" you have to specify the oglx ogly oglz of the og box; exit() load_para;"<<endl;
					exit(IOERROR);
				} 
				if(_RESTARTFlag>0) {
					cout<<" _RESTARTFlag: TRUE"<<endl;
					_ogboxlx=atof(tempvec[2].c_str());
					_ogboxly=atof(tempvec[3].c_str());
					_ogboxlz=atof(tempvec[4].c_str());
				} else {
					cout<<" _RESTARTFlag: FALSE"<<endl;
				}
				cout<<" _RESTARTboxlen: "<<_ogboxlx<<" "<<_ogboxly<<" "<<_ogboxlz<<endl;
			}
		}
		
		cout<<" now the string is: "<<tempstr<<endl;
	}
	cout<<" what happend?"<<endl;
	if(_RUNTIMES_output>_RUNTIMES_totalnum) {
		ErrorMSG("_RUNTIMES_output>_RUNTIMES_totalnum::load_parameters");
		exit(LOGICERROR);	
	}
	if(_NEIGHBOR_lj<0) {
		ErrorMSG("_NEIGHBOR_lj<0.0::load_parameters");
		exit(LOGICERROR);
	}
	if(_CoorX1>=_CoorX2 || _CoorY1>=_CoorY2 || _CoorZ1>=_CoorZ2) {
		ErrorMSG("_CoorX1>=_CoorX2 or _CoorY1>=_CoorY2 or _CoorZ1>=_CoorZ2::load_parameters");
		exit(LOGICERROR);
	}
	if(_PBC_Dim<1) {
		ErrorMSG("_PBC_Dim<1::load_parameters");
		exit(LOGICERROR);
	}
	
	if(_IFVERBOSE) {
		cout<<setiosflags(ios::left);

		cout<<setw(30)<<" ::_IFVERBOSE: "<<_IFVERBOSE<<endl;
		cout<<setw(30)<<" ::_FILENAME_conf: "<<_FILENAME_conf.c_str()<<endl;
		cout<<setw(30)<<" ::_CoorX1: "<<_CoorX1<<endl;
		cout<<setw(30)<<" ::_CoorX2: "<<_CoorX2<<endl;
		cout<<setw(30)<<" ::_CoorY1: "<<_CoorY1<<endl;
		cout<<setw(30)<<" ::_CoorY2: "<<_CoorY2<<endl;
		cout<<setw(30)<<" ::_CoorZ1: "<<_CoorZ1<<endl;
		cout<<setw(30)<<" ::_CoorZ2: "<<_CoorZ2<<endl;
		cout<<setw(30)<<" ::_PBC_Dim: "<<_PBC_Dim<<endl;
	}
	_PBL_X=_CoorX2-_CoorX1;
	_PBL_Y=_CoorY2-_CoorY1;
	_PBL_Z=_CoorZ2-_CoorZ1;
	_PBL_X_2=_PBL_X/2.0;
	_PBL_Y_2=_PBL_Y/2.0;
	_PBL_Z_2=_PBL_Z/2.0;
	if(_IFVERBOSE) {
		cout<<setw(30)<<" ::_PBL_X: "<<_PBL_X<<endl;	
		cout<<setw(30)<<" ::_PBL_Y: "<<_PBL_Y<<endl;
		cout<<setw(30)<<" ::_PBL_Z: "<<_PBL_Z<<endl;

		cout<<setw(30)<<" ::_BF_flag: "<<_BF_flag<<endl;
		cout<<setw(30)<<" ::_BF_coeff: "<<_BF_coeff<<endl;
		cout<<setw(30)<<" ::_BF_scale: "<<_BF_scale<<endl;
		cout<<setw(30)<<" ::_AG_flag: "<<_AG_flag<<endl;
		cout<<setw(30)<<" ::_DH_flag: "<<_DH_flag<<endl;
	}
	_E_highest=_E_lowest+_E_interval*_E_totalnum;
	_RG_highest=_RG_lowest+_RG_interval*_RG_totalnum;
	if(_len_EINT>0) {
		_EINT_highest=_EINT_lowest+_EINT_interval*_EINT_totalnum;
	}
	if(_IFVERBOSE) {

		cout<<setw(30)<<" ::_FFPARA_FN: "<<_FFPARA_FN<<endl;
		cout<<setw(30)<<" ::_NEIGHBOR_lj: "<<_NEIGHBOR_lj<<endl;
		cout<<setw(30)<<" ::_FLAG_ener: "<<_FLAG_ener<<endl;
		cout<<setw(30)<<" ::_FLAG_pivot: "<<_FLAG_pivot<<endl;
		cout<<setw(30)<<" ::_FLAG_rg2: "<<_FLAG_rg2<<endl;
		cout<<setw(30)<<" ::_FLAG_dis: "<<_FLAG_dis<<endl;
		//_ACCINDEX_bf=sqrt(_ACCINDEX);
		//_ACCINDEX_ag=_ACCINDEX/2.0;
		cout<<setw(30)<<" ::_ACCINDEX: "<<_ACCINDEX<<endl;
		//cout<<setw(30)<<" ::_ACCINDEX_ag: "<<_ACCINDEX_ag<<endl;
		//cout<<setw(30)<<" ::_ACCINDEX_bf: "<<_ACCINDEX_bf<<endl;

		cout<<setw(30)<<" ::_NUM_replicas: "<<_NUM_replicas<<endl;
		cout<<setw(30)<<" ::_RUNTIMES_eachstep: "<<_RUNTIMES_eachstep<<endl;
		cout<<setw(30)<<" ::_RUNTIMES_recording: "<<_RUNTIMES_recording<<endl;
		cout<<setw(30)<<" ::_RESTARTFlag: "<<_RESTARTFlag<<endl;
		if(_RESTARTFlag>0) {
			cout<<setw(30)<<" ::_ogboxlx: "<<_ogboxlx<<endl;
			cout<<setw(30)<<" ::_ogboxly: "<<_ogboxly<<endl;
			cout<<setw(30)<<" ::_ogboxlz: "<<_ogboxlz<<endl;
		}
		cout<<setw(30)<<" ::_RUNTIMES_totalnum: "<<_RUNTIMES_totalnum<<endl;
		cout<<setw(30)<<" ::_RUNTIMES_output: "<<_RUNTIMES_output<<endl;
		cout<<setw(30)<<" ::_RUNTIMES_remgap: "<<_RUNTIMES_remgap<<endl;

		cout<<setw(30)<<" ::_E_lowest: "<<_E_lowest<<endl;
		cout<<setw(30)<<" ::_E_interval: "<<_E_interval<<endl;
		cout<<setw(30)<<" ::_E_highest: "<<_E_highest<<endl;
		cout<<setw(30)<<" ::_E_totalnum: "<<_E_totalnum<<endl;
		cout<<setw(30)<<" ::_EINT_lowest: "<<_EINT_lowest<<endl;
		cout<<setw(30)<<" ::_EINT_interval: "<<_EINT_interval<<endl;
		cout<<setw(30)<<" ::_EINT_highest: "<<_EINT_highest<<endl;
		cout<<setw(30)<<" ::_EINT_totalnum: "<<_EINT_totalnum<<endl;
		cout<<setw(30)<<" ::_RG_lowest: "<<_RG_lowest<<endl;
		cout<<setw(30)<<" ::_RG_interval: "<<_RG_interval<<endl;
		cout<<setw(30)<<" ::_RG_highest: "<<_RG_highest<<endl;
		cout<<setw(30)<<" ::_RG_totalnum: "<<_RG_totalnum<<endl;

		cout<<resetiosflags(ios::left);
	}

	para_ifstream.close();
	cout<<" Parameters loaded!"<<endl;
}
void cmc::write_parameters(const string FILENAME_para) {
	if(_PROC_ID!=0) {
		return;
	}
	string OFName=FILENAME_para;
	int OFNameLen=FILENAME_para.size();
	OFName=OFName.substr(0, OFNameLen-3)+string(".cmp")+OFName.substr(OFNameLen-3, OFNameLen);
	ofstream para_ofstream( OFName.c_str() );
	if(para_ofstream==NULL) {
		cout<<" Error: Can not open file: "<<OFName<<endl;
		exit(IOERROR);
	}
	cout<<setiosflags(ios::left);

	writeparameter(para_ofstream, string("_IFVERBOSE"), _IFVERBOSE);
	writeparameter(para_ofstream, string("_FILENAME_conf"), _FILENAME_conf.c_str());
	writeparameter(para_ofstream, string("_CoorX1"), _CoorX1);
	writeparameter(para_ofstream, string("_CoorX2"), _CoorX2);
	writeparameter(para_ofstream, string("_CoorY1"), _CoorY1);
	writeparameter(para_ofstream, string("_CoorY2"), _CoorY2);
	writeparameter(para_ofstream, string("_CoorZ1"), _CoorZ1);
	writeparameter(para_ofstream, string("_CoorZ2"), _CoorZ2);
	writeparameter(para_ofstream, string("_PBC_Dim"), _PBC_Dim);
	writeparameter(para_ofstream, string("_PBL_X"), _PBL_X);
	writeparameter(para_ofstream, string("_PBL_Y"), _PBL_Y);
	writeparameter(para_ofstream, string("_PBL_Z"), _PBL_Z);
	writeparameter(para_ofstream, string("_PBL_X_2"), _PBL_X_2);
	writeparameter(para_ofstream, string("_PBL_Y_2"), _PBL_Y_2);
	writeparameter(para_ofstream, string("_PBL_Z_2"), _PBL_Z_2);
	writeparameter(para_ofstream, string("_BF_flag"), _BF_flag);
	writeparameter(para_ofstream, string("_BF_coeff"), _BF_coeff);
	writeparameter(para_ofstream, string("_BF_scale"), _BF_scale);
	writeparameter(para_ofstream, string("_AG_flag"), _AG_flag);
	writeparameter(para_ofstream, string("_DH_flag"), _DH_flag);

	if(_IFVERBOSE) {
		cout<<setw(30)<<" ::_NUM_atoms: "<<_NUM_atoms<<endl;
    }
    writeparameter(para_ofstream, string("_NUM_atoms"), _NUM_atoms);
    if(_IFVERBOSE) {
    	cout<<setw(30)<<" ::_NUM_residues: "<<_NUM_residues<<endl;
    }
    writeparameter(para_ofstream, string("_NUM_residues"), _NUM_residues);
    if(_IFVERBOSE) {
    	cout<<setw(30)<<" ::_NUM_chains: "<<_NUM_chains<<endl;
    }
    writeparameter(para_ofstream, string("_NUM_chains"), _NUM_chains);

	writeparameter(para_ofstream, string("_FFPARA_FN"), _FFPARA_FN);
	writeparameter(para_ofstream, string("_NEIGHBOR_lj"), _NEIGHBOR_lj);
	writeparameter(para_ofstream, string("_FLAG_ener"), _FLAG_ener);
	writeparameter(para_ofstream, string("_FLAG_pivot"), _FLAG_pivot);

	writeparameter(para_ofstream, string("_ACCINDEX"), _ACCINDEX);
	//writeparameter(para_ofstream, string("_ACCINDEX_ag"), _ACCINDEX_ag);
	//writeparameter(para_ofstream, string("_ACCINDEX_bf"), _ACCINDEX_bf);

	writeparameter(para_ofstream, string("_FLAG_rg2"), _FLAG_rg2);
	writeparameter(para_ofstream, string("_FLAG_dis"), _FLAG_dis);
	para_ofstream<<" _OPlist ";
	int i=0;
	for(i=0; i<_len_OP; i++) {
		para_ofstream<<_OPlist[i]<<" ";
	}
	para_ofstream<<endl;
	for(i=0; i<(_len_OP/2); i++) {
		para_ofstream<<" Group["<<setw(3)<<i+1<<"] : "<<_OPlist[2*i]<<"-"<<_OPlist[2*i+1]<<endl;
	}
	int tempN=_len_OP/2;
	//int tempnum=(tempN+1)*tempN/2;
	int j=0;
	int tempa=1;
	for(i=1; i<=tempN; i++) {
		for(j=i; j<=tempN; j++) {
			para_ofstream<<" GroupPair["<<setw(3)<<tempa++<<"] : G"<<i<<"-G"<<j<<endl;
		}
	}
	if(_len_EINT>0) {
		para_ofstream<<" _EINTlist ";
		i=0;
		for(i=0; i<_len_EINT; i++) {
			para_ofstream<<_EINTlist[i]<<" ";
		}
		para_ofstream<<endl;
		for(i=0; i<(_len_EINT/2); i++) {
			para_ofstream<<" EINT["<<setw(3)<<i+1<<"] : "<<_EINTlist[2*i]<<"-"<<_EINTlist[2*i+1]<<endl;
		}
		tempN=_len_EINT/2;
		//int tempnum=(tempN+1)*tempN/2;
		j=0;
		tempa=1;
		for(i=1; i<=tempN; i++) {
			for(j=i; j<=tempN; j++) {
				para_ofstream<<" EINTPair["<<setw(3)<<tempa++<<"] : G"<<i<<"-G"<<j<<endl;
			}
		}
	}
	//writeparameter(para_ofstream, string("_NUM_replica"), _NUM_replica);
	_RUNTIMES_eachstep*=_NUM_atoms;
	if(_IFVERBOSE) {
		cout<<setw(30)<<" ::_RUNTIMES_eachstep: "<<_RUNTIMES_eachstep<<endl;
	}
	_RUNTIMES_recording*=_NUM_atoms;
	if(_IFVERBOSE) {
		cout<<setw(30)<<" ::_RUNTIMES_recording: "<<_RUNTIMES_recording<<endl;
	}
	writeparameter(para_ofstream, string("_NUM_replicas"), _NUM_replicas);
	writeparameter(para_ofstream, string("_RUNTIMES_eachstep"), _RUNTIMES_eachstep);
	if(system("rm -fr pdbfiles")) {};
	if(system("mkdir pdbfiles")) {};
	writeparameter(para_ofstream, string("_RUNTIMES_totalnum"), _RUNTIMES_totalnum);
	writeparameter(para_ofstream, string("_RUNTIMES_recording"), _RUNTIMES_recording);
	writeparameter(para_ofstream, string("_RESTARTFlag"), _RESTARTFlag);
	if(_RESTARTFlag>0) {
		writeparameter(para_ofstream, string("_ogboxlx"), _ogboxlx);
		writeparameter(para_ofstream, string("_ogboxly"), _ogboxly);
		writeparameter(para_ofstream, string("_ogboxlz"), _ogboxlz);
	}
	writeparameter(para_ofstream, string("_RUNTIMES_output"), _RUNTIMES_output);
	writeparameter(para_ofstream, string("_RUNTIMES_remgap"), _RUNTIMES_remgap);
	//writeparameter(para_ofstream, string("_MPI_OR_NOT"), _MPI_OR_NOT);

	writeparameter(para_ofstream, string("_E_lowest"), _E_lowest);
	writeparameter(para_ofstream, string("_E_interval"), _E_interval);
	writeparameter(para_ofstream, string("_E_highest"), _E_highest);
	writeparameter(para_ofstream, string("_E_totalnum"), _E_totalnum);

	writeparameter(para_ofstream, string("_EINT_lowest"), _EINT_lowest);
	writeparameter(para_ofstream, string("_EINT_interval"), _EINT_interval);
	writeparameter(para_ofstream, string("_EINT_highest"), _EINT_highest);
	writeparameter(para_ofstream, string("_EINT_totalnum"), _EINT_totalnum);
	
	//_RG_totalnum=int(double(_RG_totalnum)/double(_NUM_chains)/double(_NUM_chains));
	_RG_interval=(_RG_highest-_RG_lowest)/_RG_totalnum;
	writeparameter(para_ofstream, string("_RG_lowest"), _RG_lowest);
	writeparameter(para_ofstream, string("_RG_interval"), _RG_interval);
	writeparameter(para_ofstream, string("_RG_highest"), _RG_highest);
	_DISSTAT_lowest=0.0;
	_DISSTAT_highest=sqrt(_PBL_X_2*_PBL_X_2+_PBL_Y_2*_PBL_Y_2+_PBL_Z_2*_PBL_Z_2);
	_DISSTAT_highest=_DISSTAT_highest*_DISSTAT_highest;
	_DISSTAT_interval=(_DISSTAT_highest-_DISSTAT_lowest)/_RG_totalnum;
	//_DISSTAT_lowest+=_DISSTAT_interval;
	//_DISSTAT_highest+=_DISSTAT_interval;
	writeparameter(para_ofstream, string("_DISSTAT_lowest"), _DISSTAT_lowest);
	writeparameter(para_ofstream, string("_DISSTAT_interval"), _DISSTAT_interval);
	writeparameter(para_ofstream, string("_DISSTAT_highest"), _DISSTAT_highest);
	tell_procid(); cout<<" rg_totnum="<<_RG_totalnum<<endl;
	tell_procid(); cout<<" _dis_low="<<_DISSTAT_lowest<<endl;
	tell_procid(); cout<<" _dis_high="<<_DISSTAT_highest<<endl;
	tell_procid(); cout<<" _dis_interval="<<_DISSTAT_interval<<endl;
	writeparameter(para_ofstream, string("_RG_totalnum"), _RG_totalnum);

	cout<<resetiosflags(ios::left);
	para_ofstream.close();
	cout<<" check [ "<<OFName<<" ] to c if paras loaded right or not."<<endl;
}
void cmc::broadcast_parameters() {
	if(_PROC_ID==0) {
		tell_procid(); cout<<" broadcasting parameters ... ";
	}
	MPI_Bcast(&_IFVERBOSE, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorX1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorX2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorY1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorY2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorZ1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_CoorZ2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBC_Dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_X, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_Y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_Z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_X_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_Y_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_PBL_Z_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_BF_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_BF_coeff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_BF_scale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_AG_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_DH_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&_NUM_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&_NUM_residues, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&_NUM_chains, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&_NEIGHBOR_lj, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_FLAG_ener, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_FLAG_pivot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_ACCINDEX, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&_ACCINDEX_ag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&_ACCINDEX_bf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_FLAG_rg2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_FLAG_dis, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_len_OP, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_len_EINT, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if( !_BF_flag && (_FLAG_pivot>1e-12) ) {
		cout<<" error: "<<endl;
		cout<<" _BF_flag="<<_BF_flag<<endl;
		cout<<" _FLAG_pivot="<<_FLAG_pivot<<endl;
		exit(LOGICERROR);
	}

	if(_PROC_ID!=0) {
		if(_len_OP!=0) {
			_OPlist=new int[_len_OP];
		}
		if(_len_EINT!=0) {
			_EINTlist=new int[_len_EINT];
		}
	}
	
	MPI_Bcast(_OPlist, _len_OP, MPI_INT, 0, MPI_COMM_WORLD);
	_len_OP=_len_OP/2;
	if(_len_EINT>0) {
		MPI_Bcast(_EINTlist, _len_EINT, MPI_INT, 0, MPI_COMM_WORLD);
		_len_EINT=_len_EINT/2;
	}
	int i=0; 
	int j=0;
	bool tempflag=false;
	cout<<" >>--------------"<<endl<<" Real Running list is: ";
	for(i=1; i<=_NUM_atoms; i++) {
		tempflag=false;
		//cout<<" "<<i<<" len_OP:"<<_len_OP<<endl;
		for(j=0; j<_len_OP; j++) {
			if( i>=_OPlist[2*j] && i<=_OPlist[2*j+1] ) {
				cout<<"["<<i<<"] "<<endl;
				//cout<<"happened"<<endl;
				tempflag=true;
				break; // fixed group, return false;
			}
		}
		if(!tempflag) {
			_Runninglist.push_back(i);
			cout<<i<<" ";
		}

	}
	cout<<endl;
	_realSize=_Runninglist.size();
	cout<<" @Proc["<<_PROC_ID<<"]: there are "<<_realSize<<" atoms moving around after kicking out fixed groups."<<endl;
	cout<<" >>--------------"<<endl;
	

	MPI_Bcast(&_NUM_replicas, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RUNTIMES_eachstep, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RUNTIMES_totalnum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RUNTIMES_recording, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RESTARTFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(_RESTARTFlag>0) {
		MPI_Bcast(&_ogboxlx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&_ogboxly, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&_ogboxlz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	MPI_Bcast(&_RUNTIMES_output, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RUNTIMES_remgap, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&_E_lowest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_E_interval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_E_highest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_E_totalnum, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&_EINT_lowest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_EINT_interval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_EINT_highest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_EINT_totalnum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&_RG_lowest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RG_interval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RG_highest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_RG_totalnum, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&_DISSTAT_lowest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_DISSTAT_interval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_DISSTAT_highest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(_PROC_ID==0) {
		cout<<" done!"<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
////////////////////////////////
void cmc::init_conformation_a(const int NUM_atoms, const int NUM_chains, const double LEN_bond, const string FILENAME_conf,
	const double xorigin, const double yorigin, const double zorigin) {
	//const double xoffset, const double yoffset, const double zoffset) {
	int i=0;
	int j=0;
	string AType=string("ATOM");
	string AChem_Name=string("C");
	int ASpec_Name=1;
	string AResi_Name=string("CHN");
	string AChn_Name=string("A");
	char AChn_Name_temp[2];//for chain index use!
	char ASpec_Name_temp[2];//for chain index use!
	int AResi_Index=0;
	double sqrt_3_d_2_BOND=sqrt(3)/2*LEN_bond;
	double sqrt_1_d_2_BOND=0.5*LEN_bond;
	int dim=int( pow(NUM_atoms*NUM_chains,1.0/3.0) )+1;
	cout<<" dim is: "<<dim<<endl;
	int temp_atoms=NUM_atoms;
	int dim_xx=0;
	int dim_yy=0;
	int dim_zz=0;
	int direction_x=1;
	int direction_y=1;
	int direction_z=1;
	for(i=0; i<NUM_chains; i++) {
		//xorigin+=xoffset;
		AChn_Name_temp[0]=AChn_Name.c_str()[0];
		AChn_Name_temp[1]='\0';
		AChn_Name=string(AChn_Name_temp);
		sprintf(ASpec_Name_temp, "%d", ASpec_Name++);
		AResi_Index++;
		AChn_Name_temp[0]=AChn_Name.c_str()[0]+1;
		AChn_Name_temp[1]='\0';
		for(j=0; j<temp_atoms; j++) {
			_system_.add_arbitrary_info(AType, j+i*temp_atoms+1, AChem_Name, string(ASpec_Name_temp), AResi_Name, AChn_Name, AResi_Index,
										sqrt_3_d_2_BOND*dim_xx+xorigin, //x
										sqrt_1_d_2_BOND*(dim_xx%2+2*dim_yy)+yorigin, //y 
										LEN_bond*dim_zz+zorigin, 1.0, false); //z
			cout<<i+1<<":"<<j+1<<" dim_x: "<<dim_xx<<" dim_y: "<<dim_yy<<" dim_z: "<<dim_zz<<endl;
			dim_xx+=direction_x;
			if( dim_xx==dim || dim_xx==-1 ) {
				direction_x=-direction_x;
				cout<<" direct_x: "<<direction_x<<endl;
				dim_xx+=direction_x;
				dim_yy+=direction_y;
				if( dim_yy==dim || dim_yy==-1 ) {
					direction_y=-direction_y;
					cout<<" direct_y: "<<direction_y<<endl;
					dim_yy+=direction_y;
					dim_zz+=direction_z;
					if( dim_zz==dim || dim_zz==-1 ) {
						direction_z=-direction_z;
						cout<<" direct_z: "<<direction_z<<endl;
						dim_zz+=direction_z;
					}
				}
			} 
		}
	}
	_system_.writepdbinfo(FILENAME_conf.c_str(), true);
	_system_.memo_free();
	cout<<endl<<" "<<NUM_atoms*NUM_chains<<" atoms was constructed on "<<NUM_chains<<" chains!"<<endl;
}
//////see loadffpara.cpp for the parameters-loading code;
///////////////////////////
void cmc::load_conformation() {//initialization conformation; init.pdb (xyz);
	if(_PROC_ID!=0) {
		return;
	}
	cout<<" loading conformation from: [ "<<_FILENAME_conf<<" ];"<<endl;
	_system_.readpdbinfo(_FILENAME_conf.c_str(), 1.0, false, _IFVERBOSE);
	_NUM_chains=_system_.nchains;
	_NUM_atoms=0;
	int j;
	int tnresidues;
	_NUM_residues=_system_.nresidues;
	for(int i=0; i<_NUM_chains; i++) {
		tnresidues=_system_.chains[i].nresidues;
		for(j=0; j<tnresidues; j++) {
			_NUM_atoms+=_system_.chains[i].residues[j].natoms;
		}
		//_NUM_residues+=tnresidues;
	}
}
///////////////////////////
void cmc::load_temperatures() {//initialization temperatures;
	if(_PROC_ID!=0) {
		return;
	}
	cout<<" loading temperatures from: [ _temperaturelist.pls ];"<<endl;
	string TempStr;
	ifstream IFSTREAM_temperatruelist("_temperaturelist.pls");
	if(IFSTREAM_temperatruelist==NULL) {
		ErrorMSG("can not open file: _temperaturelist.pls :: cmc::load_temperatures");
		exit(IOERROR);
	}
	for(int i=0; i<_NUM_replicas; i++) {
		if(!getline(IFSTREAM_temperatruelist, TempStr)) {
			ErrorMSG("no enough input : _temperaturelist.pls :: cmc::load_temperatures");
			exit(IOERROR);
		} else {
			_T_eachrep[i]=atof( Split( BFilter( TempStr ) )[0].c_str() );
			_T_rep_eachrep[i]=1.0/_T_eachrep[i];
			_Index_T_rep[i]=i;
			_Index_rep_T[i]=i;
			//cout<<setw(3)<<setiosflags(ios::left)<<"t:"<<_T_rep_eachrep[i]<<resetiosflags(ios::left)<<endl;
		}
	}
	IFSTREAM_temperatruelist.close();
}
///////////////////////////
void cmc::scatter_temperatures() {
	MPI_Bcast(_T_rep_eachrep, _NUM_replicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_T_eachrep, _NUM_replicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(_Index_rep_T, 1, MPI_INT, &_INDEX_TEMPERATURE, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Scatter(_Index_T_rep, 1, MPI_INT, &_INDEX_TEMPERATURE, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	/*if(_PROC_ID) {
		for(int i=0; i<_NUM_replicas; i++) {
			tell_procid(); cout<<setw(3)<<setiosflags(ios::left)<<"t:"
			<<_T_rep_eachrep[i]<<resetiosflags(ios::left)<<endl;
		}
	}*/
	_TEMPERATURE_REP=_T_rep_eachrep[_INDEX_TEMPERATURE];
	_TEMPERATURE=_T_eachrep[_INDEX_TEMPERATURE];
	if(_TEMPERATURE_REP<=0.0) {
		tell_procid(); cout<<"_TEMPERATURE_REP="<<_TEMPERATURE_REP<<"<=0.0"<<endl;
		cout<<"exit!"<<endl;
		exit(LOGICERROR);
	} 
	tell_procid(); cout<<" T["<<_INDEX_TEMPERATURE<<"]="<<_TEMPERATURE<<"; T_rep="<<_TEMPERATURE_REP<<";"<<endl; 
	MPI_Barrier(MPI_COMM_WORLD);
}
///////////////////////////
void cmc::broadcast_epsilonsigma() {
	MPI_Bcast(_PARA_KB.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_BOND_length.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_BOND_length2.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(_PROC_ID==0) {
		_movelength=_movelength*_BF_coeff;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&_movelength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(_PROC_ID==0 && _BF_flag) {
		int i=0;
		double tangentK=(-_BF_coeff+_BF_scale)/double(_NUM_replicas-1);
		for(i=0; i<_NUM_replicas; i++) {
			//_movelength_ea[i]=_movelength*sqrt(3.0)/3.0*(_BF_coeff+tangentK*double(i));
			_movelength_ea[i]=_movelength*(_BF_coeff+tangentK*double(i));

			//_movelength_ea[i]=_movelength*sqrt(3.0)/3.0*pow(_BF_scale, double(i)/double(_NUM_replicas-1));
			
			//_movelength_pv[i]=_PI_double*pow(_BF_scale, double(i)/double(_NUM_replicas-1));
			cout<<" @T"<<i<<": _movelen random pivot is: "<<_movelength_ea[i]*sqrt(3.0)
				       <<": movelen random=_movelen/sqrt(3)="<<_movelength_ea[i]<<endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(_movelength_ea, _NUM_replicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(_movelength_pv, _NUM_replicas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_BOND_delta.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_RMIN.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_RMAX.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_RMIN2.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_RMAX2.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_DDELTA.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_DDELTA_ea.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PARA_KA.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_THETAZERO.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_THETAZERO_cos.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_THETAZERO_sin.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PARA_KD.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PHIZERO.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PHIZERO_cos.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PHIZERO_sin.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(_EPSILON_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_LAMBDA_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_PPTYPE_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA2_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMADIS_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMAWELL_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMAWELL2_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA3_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA6_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA9_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA12_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_SIGMA24_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_E_cut_RR_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_R_cut_RR_eachres.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(_ACCINDEX==0) {
		return;
	}
	if(_PROC_ID) {
		_ELJ_ACC.Build(_NUM_residues, _NUM_residues, _ACCINDEX);
		_LJ_interval.Build(_NUM_residues, _NUM_residues);
		_ELJ_ACC.SetZero();
		_LJ_interval.SetZero();
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(_ELJ_ACC.pArray[0][0], _NUM_residues*_NUM_residues*_ACCINDEX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(_EBF_ACC.pArray[0][0], _NUM_residues*_NUM_residues*_ACCINDEX_bf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//MPI_Bcast(_BF_interval.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(_LJ_interval.pArray[0], _NUM_residues*_NUM_residues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Barrier(MPI_COMM_WORLD);
	/*if(_PROC_ID!=0) {
		tell_procid(); cout<<" _NUM_residues="<<_NUM_residues<<endl;
		for(int i=0; i<_NUM_residues; i++) {
			for(int j=0; j<_NUM_residues; j++) {
				cout<<_PARA_KB.pArray[i][j]<<" ";
			}
			cout<<endl;
		}
	}*/
}
///////////////////////////
void cmc::init_epsilonsigma() {//temporary arbitrary!
	if(_PROC_ID!=0) {
		return;
	}
	ifstream ffpara_ifstream(_FFPARA_FN.c_str());
	if(ffpara_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<_FFPARA_FN<<endl;
		exit(IOERROR);
	}
	
	string tempstr;
	vector<string> tempvec;
	
	int EP_SIG_X=_NUM_residues;
	int EP_SIG_Y=_NUM_residues;
	int xx=0;
	int yy=0;
	while( getline(ffpara_ifstream, tempstr) ) {
		if( BFilter(tempstr)==string("[PPTYPE]") ) {//start reading PPTYPE info;
			cout<<" Loading pptype from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;PPTYPE]")) {
					cout<<" PPTYPE loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" PPTYPE size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_PPTYPE_eachres.pArray[xx][yy]=atoi(tempvec[yy].c_str());
								if(_IFVERBOSE) { 
									cout<<setw(6)<<_PPTYPE_eachres.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" PPTYPE size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if( BFilter(tempstr)==string("[LAMBDA]") ) {//start reading LAMBDA info;
			cout<<" Loading LAMBDA from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;LAMBDA]")) {
					cout<<" LAMBDA loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" LAMBDA size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_LAMBDA_eachres.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) { 
									cout<<setw(6)<<_LAMBDA_eachres.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" LAMBDA size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if( BFilter(tempstr)==string("[EPSILON]") ) {//start reading epsilon info;
			cout<<" Loading epsilon from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;EPSILON]")) {
					cout<<" epsilon loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" epsilon size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_EPSILON_eachres.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) { 
									cout<<setw(6)<<_EPSILON_eachres.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" epsilon size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[SIGMA]")) {
			cout<<" Loading sigma from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;SIGMA]")) {
					cout<<" sigma loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" sigma size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_SIGMA_eachres.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) {
									cout<<setw(6)<<_SIGMA_eachres.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" sigma size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[RCUT]")) {
			cout<<" Loading r_cut from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;RCUT]")) {
					cout<<" r_cut loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" r_cut size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								if(tempvec[yy]==string("#")) {
									if(_PPTYPE_eachres.pArray[xx][yy]==2) {
										_R_cut_RR_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy];
										//_R_cut_RR_eachres.pArray[xx][yy]=pow(2.0,1.0/6.0)*_SIGMA_eachres.pArray[xx][yy];
									} else if(_PPTYPE_eachres.pArray[xx][yy]==4) {
										_R_cut_RR_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy];
										//_R_cut_RR_eachres.pArray[xx][yy]=pow(2.0,1.0/12.0)*_SIGMA_eachres.pArray[xx][yy];
									} else {
										_R_cut_RR_eachres.pArray[xx][yy]=pow(2.0,1.0/6.0)*_SIGMA_eachres.pArray[xx][yy];
									}
								} else if(tempvec[yy]==string("@")) {
									if(_PPTYPE_eachres.pArray[xx][yy]==2) {
										_R_cut_RR_eachres.pArray[xx][yy]=pow(_MAX_DOUBLE, 1.0/12.0);
									} else if(_PPTYPE_eachres.pArray[xx][yy]==4) {
										_R_cut_RR_eachres.pArray[xx][yy]=pow(_MAX_DOUBLE, 1.0/24.0);
									} else {
										_R_cut_RR_eachres.pArray[xx][yy]=pow(_MAX_DOUBLE, 1.0/12.0);
									}
								} else {
									_R_cut_RR_eachres.pArray[xx][yy]=atof(tempvec[yy].c_str());
								}
								if(_IFVERBOSE) {
									cout<<setw(6)<<_R_cut_RR_eachres.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" r_cut size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[KANGLE]")) {
			cout<<" Loading k_angle from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;KANGLE]")) {
					cout<<" k_angle loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" k_angle size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_PARA_KA.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) {
									cout<<setw(6)<<_PARA_KA.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" k_angle size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[THETA0]")) {
			cout<<" Loading theta0 from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;THETA0]")) {
					cout<<" theta0 loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" theta0 size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_THETAZERO.pArray[xx][yy]=atof(tempvec[yy].c_str())*_PI_single;
								_THETAZERO_cos.pArray[xx][yy]=cos(_THETAZERO.pArray[xx][yy]);
								_THETAZERO_sin.pArray[xx][yy]=sin(_THETAZERO.pArray[xx][yy]);
								if(_IFVERBOSE) {
									cout<<setw(6)<<_THETAZERO.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" theta0 size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[KDIHEDRAL]")) {
			cout<<" Loading k_dihedral from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;KDIHEDRAL]")) {
					cout<<" k_dihedral loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info insufficient...
							cout<<" k_dihedral size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_PARA_KD.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) {
									cout<<setw(6)<<_PARA_KD.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" k_dihedral size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[PHI0]")) {
			cout<<" Loading PHI0 from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;PHI0]")) {
					cout<<" PHI0 loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" PHI0 size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_PHIZERO.pArray[xx][yy]=atof(tempvec[yy].c_str())*_PI_single;
								_PHIZERO_cos.pArray[xx][yy]=cos(_PHIZERO.pArray[xx][yy]);
								_PHIZERO_sin.pArray[xx][yy]=sin(_PHIZERO.pArray[xx][yy]);
								if(_IFVERBOSE) {
									cout<<setw(6)<<_PHIZERO.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" PHI0 size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[KBOND]")) {
			cout<<" Loading k_bond from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;KBOND]")) {
					cout<<" k_bond loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" k_bond size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_PARA_KB.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) {
									cout<<setw(6)<<_PARA_KB.pArray[xx][yy]<<" ";
								}
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" k_bond size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[BLENGTH]")) {
			cout<<" Loading bond_length from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;BLENGTH]")) {
					cout<<" bond_length loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" bond_length size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_BOND_length.pArray[xx][yy]=atof(tempvec[yy].c_str());
								_BOND_length2.pArray[xx][yy]=_BOND_length.pArray[xx][yy]*_BOND_length.pArray[xx][yy];
								if(_movelength<_BOND_length.pArray[xx][yy]) {
									_movelength=_BOND_length.pArray[xx][yy];
								}
								if(_IFVERBOSE) {
									cout<<setw(6)<<_BOND_length.pArray[xx][yy]<<" ";
								}
							}	
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" bond_length size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else if (BFilter(tempstr)==string("[BDELTA]")) {
			cout<<" Loading bond_delta from file [ "<<_FFPARA_FN<<" ] ..."<<endl;
			while( getline(ffpara_ifstream, tempstr) ) {
				if(BFilter(tempstr)==string("[;BDELTA]")) {
					cout<<" bond_delta loaded."<<endl;
					break;
				} else {
					if( BFilter(tempstr).size()!=0 ) {//if a valid line;
						tempvec=Split(BFilter(tempstr));
						if(int(tempvec.size())<EP_SIG_Y) {// if info unsufficient...
							cout<<" bond_delta size wrong... "<<endl;
							exit(LOGICERROR);
						} else { // read info;
							if(xx>=EP_SIG_X) { continue; }
							for(yy=0; yy<EP_SIG_Y; yy++ ) {
								_BOND_delta.pArray[xx][yy]=atof(tempvec[yy].c_str());
								if(_IFVERBOSE) {
									cout<<setw(6)<<_BOND_delta.pArray[xx][yy]<<" ";
								}	
							}
							xx++;
							if(_IFVERBOSE) cout<<endl;
						}
					} else {
						continue; //read another line;
					}
				}
			}
			if(xx!=EP_SIG_X) {
				cout<<" bond_delta size wrong... "<<endl;
				exit(LOGICERROR);
			} else {
				xx=0;
			}
		} else {
			continue;
		}
	}
	
	ffpara_ifstream.close();

	for(xx=0; xx<EP_SIG_X; xx++) {
		for(yy=0; yy<EP_SIG_Y; yy++) {
			_RMIN.pArray[xx][yy]=_BOND_length.pArray[xx][yy]-_BOND_delta.pArray[xx][yy];
			_RMAX.pArray[xx][yy]=_BOND_length.pArray[xx][yy]+_BOND_delta.pArray[xx][yy];
			_RMIN2.pArray[xx][yy]=_RMIN.pArray[xx][yy]*_RMIN.pArray[xx][yy];
			_RMAX2.pArray[xx][yy]=_RMAX.pArray[xx][yy]*_RMAX.pArray[xx][yy];
			_DDELTA.pArray[xx][yy]=2.0*_BOND_delta.pArray[xx][yy];
			_DDELTA_ea.pArray[xx][yy]=_DDELTA.pArray[xx][yy]*sqrt(3.0)/3.0;
			_SIGMA2_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy]*_SIGMA_eachres.pArray[xx][yy];
			_SIGMADIS_eachres.pArray[xx][yy]=_SIGMA2_eachres.pArray[xx][yy]*_CN_coeff*_CN_coeff;
			if(_PPTYPE_eachres.pArray[xx][yy]==2) {
				_SIGMAWELL_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy];
				//_SIGMAWELL_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy]*pow(2.0,1.0/6.0);
			} else if(_PPTYPE_eachres.pArray[xx][yy]==4) {
				_SIGMAWELL_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy];
				//_SIGMAWELL_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy]*pow(2.0,1.0/12.0);
			} else {
				_SIGMAWELL_eachres.pArray[xx][yy]=_SIGMA_eachres.pArray[xx][yy]*pow(2.0,1.0/6.0);
			}
			_SIGMAWELL2_eachres.pArray[xx][yy]=_SIGMAWELL_eachres.pArray[xx][yy]*_SIGMAWELL_eachres.pArray[xx][yy];
			_SIGMA3_eachres.pArray[xx][yy]=_SIGMA2_eachres.pArray[xx][yy]*_SIGMA_eachres.pArray[xx][yy];
			_SIGMA6_eachres.pArray[xx][yy]=_SIGMA3_eachres.pArray[xx][yy]*_SIGMA3_eachres.pArray[xx][yy];
			_SIGMA9_eachres.pArray[xx][yy]=_SIGMA3_eachres.pArray[xx][yy]*_SIGMA6_eachres.pArray[xx][yy];
			_SIGMA12_eachres.pArray[xx][yy]=_SIGMA6_eachres.pArray[xx][yy]*_SIGMA6_eachres.pArray[xx][yy];
			_SIGMA24_eachres.pArray[xx][yy]=_SIGMA12_eachres.pArray[xx][yy]*_SIGMA12_eachres.pArray[xx][yy];
			_R_cut_RR_eachres.pArray[xx][yy]=_R_cut_RR_eachres.pArray[xx][yy]*_R_cut_RR_eachres.pArray[xx][yy];
			_DIS2=_R_cut_RR_eachres.pArray[xx][yy];
			//_DIS6=_DIS2*_DIS2*_DIS2;
			//_DIS12=_DIS6*_DIS6;
			_E_cut_RR_eachres.pArray[xx][yy]=Energy_LJ(_PPTYPE_eachres.pArray[xx][yy], _LAMBDA_eachres.pArray[xx][yy], _EPSILON_eachres.pArray[xx][yy],
				_SIGMA12_eachres.pArray[xx][yy], _SIGMA6_eachres.pArray[xx][yy],
				_SIGMA2_eachres.pArray[xx][yy], _SIGMA_eachres.pArray[xx][yy], _DIS2);
			//cout<<_E_cut_RR_eachres.pArray[xx][yy]<<endl;
			
			    //cout<<i<<" "<<j<<" succ!"<<endl;
			/*if( _system_.atoms[i-1].residuename==string("ROD") && _system_.atoms[j-1].residuename==string("ROD") )
			{
				_IF_RR_eachatom.pArray[i][j]=1;//rod is a bond, not impossible a single atom;
			}
			else
			{
				_IF_RR_eachatom.pArray[i][j]=0;//rod is a bond, not impossible a single atom;
			}*/
			cout<<setw(4)<<" epsilon="<<_EPSILON_eachres.pArray[xx][yy]
			             <<" sigma="<<_SIGMA_eachres.pArray[xx][yy]
			             <<" rcut="<<sqrt(_R_cut_RR_eachres.pArray[xx][yy])
			             <<" Ecut="<<_E_cut_RR_eachres.pArray[xx][yy]
			             <<" Ecut_x_LJ_og="<<4.0*_EPSILON_eachres.pArray[xx][yy]*(
			             	pow(_SIGMA_eachres.pArray[xx][yy],12.0)/pow(_R_cut_RR_eachres.pArray[xx][yy],6.0)
			               -pow(_SIGMA_eachres.pArray[xx][yy],6.0)/pow(_R_cut_RR_eachres.pArray[xx][yy],3.0))<<endl;
			//	<<" "<<_IF_RR_eachatom.pArray[xx][yy]<<endl;
		}
		//cout<<endl;
	}
	int ep_sig_x=0;
	int ep_sig_y=0;

	double sigma_max;
	int num_inteval=1000;
	double dis_inteval;
	double t_dis;
	double t_dis_2;
	//double t_dis_6;
	//double t_dis_12;
	double t_ener;
	double temp_sig6,temp_sig12,temp_sig2;;
	//double temp_r6,temp_r12;
	//int t_res=0;
	string i_name;
	string j_name;
	//int i=0;
	//int j=0;
	char tempiname[5];
	cout<<" [ ---------- gnuploting the potential profile ---------- ]"<<endl;
	ofstream pararecord("pararecord.par");
	for(ep_sig_x=0; ep_sig_x<EP_SIG_X; ep_sig_x++) {
		for(ep_sig_y=ep_sig_x; ep_sig_y<EP_SIG_Y; ep_sig_y++) {
			sprintf(tempiname,"%d", _system_.residues[ep_sig_x].residueiname);
			i_name=tempiname+_system_.residues[ep_sig_x].residuename;
			sprintf(tempiname,"%d", _system_.residues[ep_sig_y].residueiname);
			j_name=tempiname+_system_.residues[ep_sig_y].residuename;
			i_name=i_name+string("_")+j_name;
			ofstream para4gnuplot("para4gnuplot.par");
			para4gnuplot<<"fn=\'"<<i_name<<"\'\n"
			            <<"epsilon="<<showpoint<<setprecision(22)<<_EPSILON_eachres.pArray[ep_sig_x][ep_sig_y]<<"\n"
			            <<"sigma="<<showpoint<<setprecision(22)<<_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y]<<"\n"
			            <<"cut="<<showpoint<<setprecision(22)<<sqrt(_R_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y])<<endl;
			pararecord<<"fn=\'"<<i_name<<"\'\n"
			          <<"epsilon="<<showpoint<<setprecision(22)<<_EPSILON_eachres.pArray[ep_sig_x][ep_sig_y]<<"\n"
			          <<"sigma="<<showpoint<<setprecision(22)<<_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y]<<"\n"
			          <<"cut="<<showpoint<<setprecision(22)<<sqrt(_R_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y])<<"\n"<<"\n";
			para4gnuplot.close();
			i_name=i_name+string(".nrg");
			ofstream outener( i_name.c_str() );
			sigma_max=4.0*_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y];
			dis_inteval=sigma_max/num_inteval;
			temp_sig2=_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y]*_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y];
			temp_sig6=temp_sig2*temp_sig2*temp_sig2;
			temp_sig12=temp_sig6*temp_sig6;
			//temp_r6=_R_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y];
			//temp_r12=temp_r6*temp_r6;
			cout<<" blen2="<<_BOND_length2.pArray[ep_sig_x][ep_sig_y]<<endl;
			cout<<" sigw2="<<_SIGMAWELL2_eachres.pArray[ep_sig_x][ep_sig_y]<<endl;
			for(xx=1; xx<=num_inteval; xx++) {
				t_dis=double(xx)*dis_inteval;
				t_dis_2=t_dis*t_dis;
				//t_dis_6=t_dis_2*t_dis_2*t_dis_2;
				//t_dis_12=t_dis_6*t_dis_6;
				t_ener=0.0;
				if(t_dis_2<_R_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y]) {
					t_ener=Energy_LJ(_PPTYPE_eachres.pArray[ep_sig_x][ep_sig_y],_LAMBDA_eachres.pArray[ep_sig_x][ep_sig_y],
						_EPSILON_eachres.pArray[ep_sig_x][ep_sig_y],
						temp_sig12,temp_sig6,temp_sig2,_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y],t_dis_2)
				    	  -_E_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y];	
				}
				outener<<t_dis<<"\t"<<t_ener<<"\t"<<Energy_BF(_PARA_KB.pArray[ep_sig_x][ep_sig_y], t_dis_2, 
					_BOND_length.pArray[ep_sig_x][ep_sig_y],
					_BOND_length2.pArray[ep_sig_x][ep_sig_y], 
					_EPSILON_eachres.pArray[ep_sig_x][ep_sig_y],
					_SIGMAWELL2_eachres.pArray[ep_sig_x][ep_sig_y], temp_sig6)<<"\n";
			}
			outener.close();
			if(system("gnuplot < draw.gpl")) {};
		}
	}
	pararecord.close();
	cout<<" [ ---------- gnuploting the potential profile ---------- ]"<<endl;
	cout<<" [ NOTE ]: you can check your energy profile now... those .eps files. "<<endl;
	cout<<" "<<endl;

	if(_ACCINDEX==0) { return; }
	
	//_EAG_ACC.Build(EP_SIG_X, EP_SIG_Y, _ACCINDEX_ag);
	//_EBF_ACC.Build(EP_SIG_X, EP_SIG_Y, _ACCINDEX_bf);
	_ELJ_ACC.Build(EP_SIG_X, EP_SIG_Y, _ACCINDEX);
	//_EDH_ACC.Build(EP_SIG_X, EP_SIG_Y, _ACCINDEX);
	//_BF_interval.Build(EP_SIG_X, EP_SIG_Y);
	_LJ_interval.Build(EP_SIG_X, EP_SIG_Y);
	//_EAG_ACC.SetZero();
	//_EBF_ACC.SetZero();
	_ELJ_ACC.SetZero();
	//_EDH_ACC.SetZero();
	//_BF_interval.SetZero();
	_LJ_interval.SetZero();
	//_AG_interval=_PI_single/_ACCINDEX_ag;
	//_DH_interval=_PI_double/_ACCINDEX;
	
	for(ep_sig_x=0; ep_sig_x<EP_SIG_X; ep_sig_x++) {
		for(ep_sig_y=0; ep_sig_y<EP_SIG_Y; ep_sig_y++) {
			//_BF_interval.pArray[ep_sig_x][ep_sig_y]=_DDELTA.pArray[ep_sig_x][ep_sig_y]/_ACCINDEX_bf; 
			_LJ_interval.pArray[ep_sig_x][ep_sig_y]=_R_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y]/_ACCINDEX;
			
			//1.LJ
			for(xx=0; xx<_ACCINDEX; xx++) {
				t_dis_2=(double(xx)+0.5)*_LJ_interval.pArray[ep_sig_x][ep_sig_y];				
				_ELJ_ACC.pArray[ep_sig_x][ep_sig_y][xx]=Energy_LJ(
					_PPTYPE_eachres.pArray[ep_sig_x][ep_sig_y],
					_LAMBDA_eachres.pArray[ep_sig_x][ep_sig_y],
					_EPSILON_eachres.pArray[ep_sig_x][ep_sig_y],
					_SIGMA12_eachres.pArray[ep_sig_x][ep_sig_y],
					_SIGMA6_eachres.pArray[ep_sig_x][ep_sig_y],
					_SIGMA2_eachres.pArray[ep_sig_x][ep_sig_y],
					_SIGMA_eachres.pArray[ep_sig_x][ep_sig_y],
					t_dis_2)-_E_cut_RR_eachres.pArray[ep_sig_x][ep_sig_y];	
				
			}
			//.BF
			/*for(xx=0; xx<_ACCINDEX_bf; xx++) {
				t_dis=_RMIN.pArray[ep_sig_x][ep_sig_y]+(double(xx)+0.5)*_BF_interval.pArray[ep_sig_x][ep_sig_y];				
				_EBF_ACC.pArray[ep_sig_x][ep_sig_y][xx]=Energy_BF(_PARA_KB.pArray[ep_sig_x][ep_sig_y], 
					t_dis, _BOND_length.pArray[ep_sig_x][ep_sig_y]);
				
			}*/
			//.AG
			//for(xx=0; xx<_ACCINDEX_ag; xx++) {
			//	t_dis=(double(xx)+0.5)*_AG_interval.pA;				
				/*_EAG_ACC.pArray[ep_sig_x][ep_sig_y][xx]=Energy_AG(_PARA_KB.pArray[ep_sig_x][ep_sig_y], 
					t_dis, _BOND_length.pArray[ep_sig_x][ep_sig_y]);*/
			//	_EAG_ACC.pArray[ep_sig_x][ep_sig_y][xx]=cos(xx);				
			//}
			//.DH
			//for(xx=0; xx<_ACCINDEX; xx++) {
			//	t_dis=(double(xx)+0.5)*_PI_;				
				/*_EAG_ACC.pArray[ep_sig_x][ep_sig_y][xx]=Energy_AG(_PARA_KB.pArray[ep_sig_x][ep_sig_y], 
					t_dis, _BOND_length.pArray[ep_sig_x][ep_sig_y]);*/
			//	_EAG_ACC.pArray[ep_sig_x][ep_sig_y][xx]=cos(xx);				
			//}
		}
	}
}
///////////////////////////
///////////////////////////
bool cmc::make_change() {
	if( _INDEX_chn_or_not ) {// if atom belongs to a chain;
		if(_FLAG_pivot>_TEMPERATURE) { // not pivot algorithm;
		//if( (_INDEX_TEMPERATURE%2) != 0 ) {
			if( _TYPE_atom_ind_real == 0 ) {
				tempddelta=_movelength_ea[_INDEX_TEMPERATURE];
				//temprandom=int(rand_seed(iseed_len2)*3.00);
				temprandom=rand()%3;
				if(temprandom==0) {
					_XX[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_XX[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				} else if(temprandom==1) {
					_YY[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_YY[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				} else {
					_ZZ[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_ZZ[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				}
				//some rules here;
				dc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen];
				dc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen];
				dc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen];
				dd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
				if(    dd>_RMAX2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]]
				    || dd<_RMIN2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]] ) {
					if(temprandom==0) {
						_XX[_INDEX_chosen]=_XX[0];
					} else if(temprandom==1) {
						_YY[_INDEX_chosen]=_YY[0];
					} else {
						_ZZ[_INDEX_chosen]=_ZZ[0];
					}
					_MC_NUM_FIL_stat[_INDEX_TEMPERATURE]+=1.0;
					return false;
				}
				ec[0]=_XX[_INDEX_chosen]-_XX[_INDEX_chosen-1];
				ec[1]=_YY[_INDEX_chosen]-_YY[_INDEX_chosen-1];
				ec[2]=_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-1];
				ee=ec[0]*ec[0]+ec[1]*ec[1]+ec[2]*ec[2];
				if(    ee>_RMAX2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-1]]
				    || ee<_RMIN2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-1]] ) {
					if(temprandom==0) {
						_XX[_INDEX_chosen]=_XX[0];
					} else if(temprandom==1) {
						_YY[_INDEX_chosen]=_YY[0];
					} else {
						_ZZ[_INDEX_chosen]=_ZZ[0];
					}
					_MC_NUM_FIL_stat[_INDEX_TEMPERATURE]+=1.0;
					return false;
				}
				//return false; 
			} else {
				//tempddelta=_DDELTA_ea.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]];
				tempddelta=_movelength_ea[_INDEX_TEMPERATURE];
				//temprandom=int(rand_seed(iseed_len2)*3.00);
				temprandom=rand()%3;
				if(temprandom==0) {
					_XX[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_XX[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				} else if(temprandom==1) {
					_YY[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_YY[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				} else {
					_ZZ[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
					//_ZZ[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
				}
				dc[0]=_XX[_INDEX_chosen]-_XX[_INDEX_chosen-_TYPE_atom_ind_real];
				dc[1]=_YY[_INDEX_chosen]-_YY[_INDEX_chosen-_TYPE_atom_ind_real];
				dc[2]=_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-_TYPE_atom_ind_real];
				dd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
				if(    dd>_RMAX2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]]
				    || dd<_RMIN2.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]] ) {
					if(temprandom==0) {
						_XX[_INDEX_chosen]=_XX[0];
					} else if(temprandom==1) {
						_YY[_INDEX_chosen]=_YY[0];
					} else {
						_ZZ[_INDEX_chosen]=_ZZ[0];
					}
					_MC_NUM_FIL_stat[_INDEX_TEMPERATURE]+=1.0;
					return false;
				}
			}
		} else { // pivot algorithm;
			Omega=rand_seed(iseed_angle)*_PI_double;//*_movelength_pv[_INDEX_TEMPERATURE];
			sin_Omega=sin(Omega);
			cos_Omega=cos(Omega);	
			if( _TYPE_atom_ind_real == 0 ) {// if atom is in the middle of the chain;  	
				//vector AB -> ;  
				dc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen-1];
				dc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen-1];
				dc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen-1];
				//dd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
				dd=_DIS2_eachatom.pArray[_INDEX_chosen+1][_INDEX_chosen-1];
				dd_sqrt=sqrt(dd);

				ec[0]=_XX[_INDEX_chosen]-_XX[_INDEX_chosen-1];
				ec[1]=_YY[_INDEX_chosen]-_YY[_INDEX_chosen-1];
				ec[2]=_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-1];
				//ee=ec[0]*ec[0]+ec[1]*ec[1]+ec[2]*ec[2];
				ee=_DIS2_eachatom.pArray[_INDEX_chosen][_INDEX_chosen-1];

				/*fc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen];
				fc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen];
				fc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen];
				ff=fc[0]*fc[0]+fc[1]*fc[1]+fc[2]*fc[2];*/
				ff=_DIS2_eachatom.pArray[_INDEX_chosen+1][_INDEX_chosen];

				//just for safe, although this will be impossible in a real system.
				//printf(" dd=%32.30f\n", dd);
				if(dd<1e-12) {
					/*cout<<"TempEner="<<TempEner<<endl;
					cout<<"_ENER_total="<<_ENER_total<<endl;
					cout<<"tempindex_judge="<<tempindex_judge<<endl;
					cout<<"tempindex_judge_bak="<<tempindex_judge_bak<<endl;*/
					factor=0.5;
					cout<<" factor="<<factor<<" truly happened!"<<endl;
				} else {
					factor=(ee-ff)/2.0/dd+0.5;
				}
				//vector AD -> ;
				uc=dc[0]*factor;//+_XX[_INDEX_chosen-1];
				vc=dc[1]*factor;//+_YY[_INDEX_chosen-1];
				wc=dc[2]*factor;//+_ZZ[_INDEX_chosen-1];
				
				//vector DC -> ; DC = AC - AD ;  
				vecDC[0]=ec[0]-uc;
				vecDC[1]=ec[1]-vc;
				vecDC[2]=ec[2]-wc;
				//cout<<" e+f="<<sqrt(ee)+sqrt(ff)<<" d="<<sqrt(dd)<<endl;

				if(_BF_flag) { // without bond fluctuation
					len_vecDC=sqrt(vecDC[0]*vecDC[0]+vecDC[1]*vecDC[1]+vecDC[2]*vecDC[2]);
					
					//
					ee=max(dd_sqrt-_RMAX.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]],
						_RMIN.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-1]]
					  +rand_seed(iseed_len1)*_DDELTA.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-1]]);

					r2max=min((dd_sqrt+ee), _RMAX.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]]);
					r2min=max(fabs(dd_sqrt-ee), _RMIN.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]]);
					if(r2max>=r2min) {
						//cout<<" dd="<<dd_sqrt<<", ee="<<ee<<", r2max="<<r2max<<", r2min="<<r2min<<endl;
						ff=r2min+rand_seed(iseed_len1)*(r2max-r2min);
						//cout<<ff<<endl;
						////getchar();
						//if(ee>)
					} else { //ee too short. change it.
						//ff=r2max;
						//ee=dd_sqrt
						/*cout<<" i="<<_INDEX_chosen<<endl;
						cout<<" d="<<dd_sqrt<<endl;
						cout<<" e="<<ee<<endl;
						cout<<" rmax="<<_RMAX.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]]<<endl;
						cout<<" vs.. "<<(dd_sqrt+ee)<<endl;
						cout<<" rmin="<<_RMIN.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen+1]]<<endl;
						cout<<" vs.. "<<fabs(dd_sqrt-ee)<<endl;
						cout<<" error: r2max="<<r2max<<"; r2min="<<r2min<<endl;*/
						_MC_NUM_FIL_stat[_INDEX_TEMPERATURE]+=1.0;
						return false; 
					}

					//if(dd_sqrt>ee+ff) {
					//	cout<<"d>e+f truly happened"<<endl;
					//}
					ee=ee*ee;
					ff=ff*ff;

					if(dd<1e-12) {
						factor=0.5;
						cout<<" factor="<<factor<<" truly happened!"<<endl;
					//} else if( (ee+dd) < ff ) { //e side obtuse angle;
						//factor=(ee+dd-ff)/2.0/dd;
					//} else if( (ff+dd) < ee ) { //f side obtuse angle;
						//factor=(ee-ff+dd)/2.0/dd;
					} else {
						//factor=(ee-ff+dd)/2.0/dd;
						factor=(ee-ff)/2.0/dd+0.5;
					}
					
					tl=sqrt(ee-dd*factor*factor);
					/*cout<<" ee="<<ee<<" e="<<sqrt(ee)<<endl;
					cout<<" ff="<<ff<<" f="<<sqrt(ff)<<endl;
					cout<<" dd="<<dd<<" d="<<dd_sqrt<<endl;
					cout<<" factor="<<factor<<endl;
					cout<<" g= "<<factor*dd_sqrt<<endl;
					cout<<" tl= "<<tl<<endl;*/
					if(tl!=tl) {//when e+f=d; tl=nan; this is important
						tl=0.0;
					}
					//cout<<" veclen="<<len_vecDC<<" :: tl="<<tl<<endl;
					//vector AD -> ;
					uc=dc[0]*factor;//+_XX[_INDEX_chosen-1];
					vc=dc[1]*factor;//+_YY[_INDEX_chosen-1];
					wc=dc[2]*factor;//+_ZZ[_INDEX_chosen-1];
					//cout<<"uc="<<uc<<endl;
					//cout<<"vc="<<vc<<endl;
					//cout<<"wc="<<wc<<endl;

					//vector DC -> ; DC = AC - AD ;  
					if(len_vecDC>1e-12) {
						vecDC[0]=vecDC[0]*tl/len_vecDC;
						vecDC[1]=vecDC[1]*tl/len_vecDC;
						vecDC[2]=vecDC[2]*tl/len_vecDC;
					} 
				}
				dd_sqrt=sin_Omega/dd_sqrt; //not sqrt(d) but for the convenience for next step calc;

				//coorRES=Add_p(C_p(cos_Omega,vecDC), C_p(sin_Omega/sqrt(dd), p_X_p(dc, vecDC)));
				if(len_vecDC>1e-12) {
					_XX[_INDEX_chosen]=_XX[_INDEX_chosen-1]+uc+vecDC[0]*cos_Omega+(dc[1]*vecDC[2]-dc[2]*vecDC[1])*dd_sqrt;
					_YY[_INDEX_chosen]=_YY[_INDEX_chosen-1]+vc+vecDC[1]*cos_Omega+(dc[2]*vecDC[0]-dc[0]*vecDC[2])*dd_sqrt;
					_ZZ[_INDEX_chosen]=_ZZ[_INDEX_chosen-1]+wc+vecDC[2]*cos_Omega+(dc[0]*vecDC[1]-dc[1]*vecDC[0])*dd_sqrt;
				} else {
					_MC_NUM_FIL_stat[_INDEX_TEMPERATURE]+=1.0;
					return false; 
				}
				
				/**/
				double tempee=(_XX[_INDEX_chosen]-_XX[_INDEX_chosen-1])*(_XX[_INDEX_chosen]-_XX[_INDEX_chosen-1])
								+(_YY[_INDEX_chosen]-_YY[_INDEX_chosen-1])*(_YY[_INDEX_chosen]-_YY[_INDEX_chosen-1])
								+(_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-1])*(_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-1]);
				double tempff=(_XX[_INDEX_chosen]-_XX[_INDEX_chosen+1])*(_XX[_INDEX_chosen]-_XX[_INDEX_chosen+1])
								+(_YY[_INDEX_chosen]-_YY[_INDEX_chosen+1])*(_YY[_INDEX_chosen]-_YY[_INDEX_chosen+1])
								+(_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen+1])*(_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen+1]);
				double tempdd=(_XX[_INDEX_chosen-1]-_XX[_INDEX_chosen+1])*(_XX[_INDEX_chosen-1]-_XX[_INDEX_chosen+1])
								+(_YY[_INDEX_chosen-1]-_YY[_INDEX_chosen+1])*(_YY[_INDEX_chosen-1]-_YY[_INDEX_chosen+1])
								+(_ZZ[_INDEX_chosen-1]-_ZZ[_INDEX_chosen+1])*(_ZZ[_INDEX_chosen-1]-_ZZ[_INDEX_chosen+1]);
				if( (fabs(tempee-ee)>1e-6) || (fabs(tempff-ff)>1e-6) ) {
					tell_procid();
					cout<<endl;
					cout<<_I_totalnum<<" :: "<<_I_eachstep<<endl;
					cout<<_INDEX_chosen<<endl;
					cout<<" veclen="<<sqrt(vecDC[0]*vecDC[0]+vecDC[1]*vecDC[1]+vecDC[2]*vecDC[2])<<endl;
					cout<<" after "<<endl;
					cout<<" tl_calc="<<sqrt(
						(_XX[_INDEX_chosen-1]+uc-_XX[_INDEX_chosen])*(_XX[_INDEX_chosen-1]+uc-_XX[_INDEX_chosen])
						+(_YY[_INDEX_chosen-1]+vc-_YY[_INDEX_chosen])*(_YY[_INDEX_chosen-1]+vc-_YY[_INDEX_chosen])
						+(_ZZ[_INDEX_chosen-1]+wc-_ZZ[_INDEX_chosen])*(_ZZ[_INDEX_chosen-1]+wc-_ZZ[_INDEX_chosen]))<<endl;
					cout<<" veclen="<<sqrt(vecDC[0]*vecDC[0]+vecDC[1]*vecDC[1]+vecDC[2]*vecDC[2])<<endl;
					cout<<" g="<<sqrt(uc*uc+vc*vc+wc*wc)<<endl;
					double realfactor=(tempee-tempff)/2.0/tempdd+0.5;
					cout<<" realfactor="<<realfactor<<" factor="<<factor<<endl;
					cout<<" real e="<<sqrt(tempee)<<" e="<<sqrt(ee)<<endl;
					cout<<" real f="<<sqrt(tempff)<<" f="<<sqrt(ff)<<endl;
					cout<<" real t="<<sqrt(tempee-dd*realfactor*realfactor)<<" t="<<tl<<endl;
					cout<<" temp t="<<sqrt(tempee-dd*factor*factor)<<endl;
					cout<<" e+f="<<sqrt(ee)+sqrt(ff)<<" d="<<sqrt(dd)<<endl;
					cout<<" real e+f="<<sqrt(tempee)+sqrt(tempff)<<endl;
					
					//cout<<"uc="<<uc<<endl;
					//cout<<"vc="<<vc<<endl;
					//cout<<"wc="<<wc<<endl;
					exit(LOGICERROR);
				}
				
			} else {// if atom is at the front or at the end of the chain; 
				THETA=asin(rand_seed(iseed_angle)*2.0-1.0)+_PI_half;
				//cout<<" THE="<<THE<<endl;
				sin_THETA=sin(THETA);
				cos_THETA=cos(THETA);

				/*tc[0]=_XX[_INDEX_chosen-_TYPE_atom_ind]-_XX[_INDEX_chosen];
				tc[1]=_YY[_INDEX_chosen-_TYPE_atom_ind]-_YY[_INDEX_chosen];
				tc[2]=_ZZ[_INDEX_chosen-_TYPE_atom_ind]-_ZZ[_INDEX_chosen];

				tl=sqrt(tc[0]*tc[0]+tc[1]*tc[1]+tc[2]*tc[2]);*/
				if(_BF_flag) {
					tl=_RMIN.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]]
					  +rand_seed(iseed_len1)*_DDELTA.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]];
				} else  {
					tl=_BOND_length.pArray[_INDEX_RES_ATM[_INDEX_chosen]][_INDEX_RES_ATM[_INDEX_chosen-_TYPE_atom_ind_real]];
				} 
				//cout<<" tl="<<tl<<endl;
				sin_THETA=sin_THETA*tl;
				_XX[_INDEX_chosen]=_XX[_INDEX_chosen-_TYPE_atom_ind_real]+cos_Omega*sin_THETA;
				_YY[_INDEX_chosen]=_YY[_INDEX_chosen-_TYPE_atom_ind_real]+sin_Omega*sin_THETA;
				_ZZ[_INDEX_chosen]=_ZZ[_INDEX_chosen-_TYPE_atom_ind_real]+cos_THETA*tl;
			}
		}
	} else {// if atom does not belongs to a chain;	
		/*tempddelta=_movelength_ea[_INDEX_TEMPERATURE];
		_XX[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
		_ZZ[_INDEX_chosen]+=(rand_seed(iseed_len2)-0.5)*tempddelta;
		_YY[_INDEX_chosen]+=(rand_seed(iseed_len3)-0.5)*tempddelta;*/
		tempddelta=_movelength_ea[_INDEX_TEMPERATURE];
		//temprandom=int(rand_seed(iseed_len2)*3.00);
		temprandom=rand()%3;
		if(temprandom==0) {
			_XX[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
			//_XX[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
		} else if(temprandom==1) {
			_YY[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
			//_YY[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
		} else {
			_ZZ[_INDEX_chosen]+=(rand_seed(iseed_len1)-0.5)*tempddelta;
			//_ZZ[_INDEX_chosen]+=(double(rand()%1000000)/1e6-0.5)*tempddelta;
		}
	}
	//_MC_NUM_TOT+=1; /*************  for check program ************/
	_MC_NUM_TOT_stat[_INDEX_TEMPERATURE]+=1.0;
	return true;
}
///////////////////////////////////////////////////////////////
///////////////////////////
bool cmc::chck_bond_len() {//after each make mapping!
	if(!_Enery_initialized) {
		cout<<" energy not initialized, can not check bond length..."<<endl;
		return LOGICERROR;
	}
	int i=0;
	double neighbor_len=0.0;
	//tell_procid();
	for(i=1; i<_NUM_atoms; i++) {
		if( _INDEX_chn_or_not_atom[i+1] && _INDEX_chn_or_not_atom[i] && _INDEX_CHN_HEAD_real[_INDEX_chn_atom[i+1]]<=i) {
			neighbor_len=sqrt(_DIS2_eachatom.pArray[i][i+1]);
			if( !_BF_flag ) {
				if( fabs(neighbor_len-_BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]])>2e-2 ) {
					tell_procid();
					cout<<endl<<" distance between "<<i<<" and "<<i+1<<" atoms has some problem: "<<endl;
					cout<<"_INDEX_chn_or_not_atom[i+1]="<<_INDEX_chn_or_not_atom[i+1]<<endl;
					cout<<"_INDEX_chn_or_not_atom[i]="<<_INDEX_chn_or_not_atom[i]<<endl;
					cout<<"_INDEX_CHN_HEAD_real[i+1]="<<_INDEX_CHN_HEAD_real[_INDEX_chn_atom[i+1]]<<endl;
					printf(" neighbor_len=%2.18lf\n", neighbor_len);
					printf(" _DIS2=%4.18lf\n", _DIS2_eachatom.pArray[i][i+1]);
					printf(" _BOND_length=%2.18lf\n", _BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]]);
					printf(" neighbor_len-_BOND_length=%e\n", neighbor_len-_BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]]);
					cout<<" error: |neighbor_len-_BOND_length|>2e-2, check your program."<<endl;
					exit(-1);
				}
			} else {
				if( fabs(neighbor_len-_BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]])-_BOND_delta.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]]>1e-2 ) {
					tell_procid();
					cout<<_XX[i]<<endl;
					cout<<_YY[i]<<endl;
					cout<<_ZZ[i]<<endl;
					cout<<_XX[i+1]<<endl;
					cout<<_YY[i+1]<<endl;
					cout<<_ZZ[i+1]<<endl;
					cout<<_DIS2_eachatom.pArray[i][i+1]<<endl;

					cout<<endl<<" distance between "<<i<<" and "<<i+1<<" atoms has some problem: "<<endl;
					cout<<"_INDEX_chn_or_not_atom[i+1]="<<_INDEX_chn_or_not_atom[i+1]<<endl;
					cout<<"_INDEX_chn_or_not_atom[i]="<<_INDEX_chn_or_not_atom[i]<<endl;
					cout<<"_INDEX_CHN_HEAD_real[i+1]="<<_INDEX_CHN_HEAD_real[_INDEX_chn_atom[i+1]]<<endl;
					printf(" neighbor_len=%e\n", neighbor_len);
					printf(" _BOND_length=%e\n", _BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]]);
					printf(" neighbor_len-_BOND_length=%e\n", neighbor_len-_BOND_length.pArray[_INDEX_RES_ATM[i]][_INDEX_RES_ATM[i+1]]);
					cout<<" error: |neighbor_len-_BOND_length|>_BOND_delta, check your program."<<endl;
					exit(-1);
				}
			}
		}
	}
	return true;
}
///////////////////////////
void cmc::make_choice(const int INDEX_coor) {	
	if( (INDEX_coor>_NUM_atoms) || (INDEX_coor<=0) ) {
		cout<<"Error: INDEX_coor="<<INDEX_coor<<", should be in [1, "<<_NUM_atoms<<"]"<<endl;
		exit(LOGICERROR);
	}
	_INDEX_chosen           = INDEX_coor;
	_INDEX_chn_ind          =_INDEX_chn_atom[_INDEX_chosen];        // for the chosen one;
	_SIZE_of_chosen_chn     =_SIZE_of_chn[_INDEX_chn_ind];          // for the chosen one; * different!
	_SIZE_of_chosen_chn_real     =_SIZE_of_chn_real[_INDEX_chn_ind];          // for the chosen one; * different!
	_INDEX_chnhead          =_INDEX_CHN_HEAD[_INDEX_chn_ind];
	_INDEX_chntail          =_INDEX_CHN_TAIL[_INDEX_chn_ind]+1;
	_INDEX_chnhead_real     =_INDEX_CHN_HEAD_real[_INDEX_chn_ind];
	_INDEX_chntail_real     =_INDEX_CHN_TAIL_real[_INDEX_chn_ind]+1;
	_INDEX_chn_or_not       =_INDEX_chn_or_not_atom[_INDEX_chosen]; // for the chosen one;
	_INDEX_atm_in_chn_ind   =_INDEX_atm_in_chn[_INDEX_chosen]; // for the chosen one;
	_INDEX_atm_in_chn_ind_real   =_INDEX_atm_in_chn_real[_INDEX_chosen]; // for the chosen one;
	_TYPE_atom_ind          =_TYPE_atom[_INDEX_chosen];        // for the chosen one;
	_TYPE_atom_ind_real     =_TYPE_atom_real[_INDEX_chosen];        // for the chosen one;
	_Index_lneighbor_ind    =_INDEX_LNEIGHBOR[_INDEX_chosen];
	_Index_rneighbor_ind    =_INDEX_RNEIGHBOR[_INDEX_chosen]+1;
	_Index_lneighbor_ind_real    =_INDEX_LNEIGHBOR_real[_INDEX_chosen];
	_Index_rneighbor_ind_real    =_INDEX_RNEIGHBOR_real[_INDEX_chosen]+1;
	/*if(_INDEX_chosen==21) {
		cout<<_Index_lneighbor_ind<<"  "<<_Index_lneighbor_ind<<endl;
	}*/
	_Res_chosen             =_INDEX_RES_ATM[_INDEX_chosen];
	//// [preparation for change...]
	_XX[0]=_XX[_INDEX_chosen];
	_YY[0]=_YY[_INDEX_chosen];
	_ZZ[0]=_ZZ[_INDEX_chosen];
	//calc_energy_chosen_atom();// not necessary any more in this detailed energy edition !!!!!!!!!!!!!!!!!
	int i;
	//#pragma ompparallel for
	for(i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_chosenatom_backup[i]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_chosenatom_backup[i]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_chosenatom_backup[i]=_DIS_x_eachatom.pArray[_INDEX_chosen][i];
		_DIS_y_chosenatom_backup[i]=_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_chosenatom_backup[i]=_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}	
	/*for(i=_INDEX_chosen+1; i<; i++) {
		_ENER_LJ_chosenatom_backup[i]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_chosenatom_backup[i]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_chosenatom_backup[i]=_DIS_x_eachatom.pArray[_INDEX_chosen][i];
		_DIS_y_chosenatom_backup[i]=_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_chosenatom_backup[i]=_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}*/	
	if(_BF_flag) {
		_ENER_BF_chosenatom_backup_l=_ENER_BF_eachatom[_INDEX_chosen-1];
		_ENER_BF_chosenatom_backup_r=_ENER_BF_eachatom[_INDEX_chosen];
	}
	if(_AG_flag) {
		_ENER_AG_chosenatom_backup_l=_ENER_AG_eachatom[_INDEX_chosen-1];
		_ENER_AG_chosenatom_backup_c=_ENER_AG_eachatom[_INDEX_chosen];
		_ENER_AG_chosenatom_backup_r=_ENER_AG_eachatom[_INDEX_chosen+1];
	}
	if(_DH_flag) {
		_ENER_DH_chosenatom_backup_1=_ENER_DH_eachatom[_INDEX_chosen-1];
		_ENER_DH_chosenatom_backup_2=_ENER_DH_eachatom[_INDEX_chosen];
		_ENER_DH_chosenatom_backup_3=_ENER_DH_eachatom[_INDEX_chosen+1];
		_ENER_DH_chosenatom_backup_4=_ENER_DH_eachatom[_INDEX_chosen+2];
	}
}
///////////////////////////
void cmc::make_judge() {	
    calc_energy_chosen_atom();
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    _ENER_delta_LJ=0.0;
	_ENER_delta_LJ_intra=0.0;
	for(i=0; i<_NUM_chains; i++) {
		_ENER_delta_LJ_inter[i]=0.0;
	}
	/*_ENER_delta_BF_l=0.0;
	_ENER_delta_BF_r=0.0;
	_ENER_delta_AG_l=0.0;
	_ENER_delta_AG_c=0.0;
	_ENER_delta_AG_r=0.0;
	_ENER_delta_DH_1=0.0;
	_ENER_delta_DH_2=0.0;
	_ENER_delta_DH_3=0.0;
	_ENER_delta_DH_4=0.0;*/
	double temp=0.0;
	//cout<<endl;
	/*for(i=0;i<_NUM_chains;i++) {
		if(i!=_INDEX_chn_ind) {
			k=_INDEX_CHN_HEAD[i];
			l=_INDEX_CHN_TAIL[i];	
			for(j=k; j<=l; j++) {
				_ENER_delta_LJ_inter[i]+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][j]-_ENER_LJ_chosenatom_backup[j];
			}
			//cout<<" Einter["<<i<<"]="<<_ENER_delta_LJ_inter[i]<<" @ make_judge() "<<endl;
			temp+=_ENER_delta_LJ_inter[i];
		} else {
			for(j=_INDEX_chnhead; j<_INDEX_chntail; j++) {
				_ENER_delta_LJ_intra+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][j]-_ENER_LJ_chosenatom_backup[j];
			}
			temp+=_ENER_delta_LJ_intra;
			//cout<<" Eintra["<<_INDEX_chn_ind<<"]="<<_ENER_delta_LJ_intra<<" @ make_judge() "<<endl;
		}
	}*/
	//#pragma ompparallel for
	for(j=_INDEX_chnhead; j<_INDEX_chntail; j++) {
		_ENER_delta_LJ_intra+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][j]-_ENER_LJ_chosenatom_backup[j];
	}
	temp+=_ENER_delta_LJ_intra;
	for(i=0;i<_INDEX_chn_ind;i++) {
		k=_INDEX_CHN_HEAD[i];
		l=_INDEX_CHN_TAIL[i];	
		//#pragma ompparallel for
		for(j=k; j<=l; j++) {
			_ENER_delta_LJ_inter[i]+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][j]-_ENER_LJ_chosenatom_backup[j];
		}
		//cout<<" Einter["<<i<<"]="<<_ENER_delta_LJ_inter[i]<<" @ make_judge() "<<endl;
		temp+=_ENER_delta_LJ_inter[i];
	} 
	for(i=_INDEX_chn_ind+1;i<_NUM_chains;i++) {
		k=_INDEX_CHN_HEAD[i];
		l=_INDEX_CHN_TAIL[i];	
		//#pragma ompparallel for
		for(j=k; j<=l; j++) {
			_ENER_delta_LJ_inter[i]+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][j]-_ENER_LJ_chosenatom_backup[j];
		}
		//cout<<" Einter["<<i<<"]="<<_ENER_delta_LJ_inter[i]<<" @ make_judge() "<<endl;
		temp+=_ENER_delta_LJ_inter[i];
	}
	/*for(i=1; i<_SIZE_memo; i++) {
		_ENER_delta_LJ+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i]-_ENER_LJ_chosenatom_backup[i];
	}
	if( fabs(temp-_ENER_delta_LJ)>1e-12 ) {
		cout<<" temp="<<temp<<endl;
		cout<<" _ENER_delta_LJ="<<_ENER_delta_LJ<<endl;
		cout<<" _ENER_delta_LJ-temp="<<_ENER_delta_LJ-temp<<endl;
		cout<<" something wrong ... "<<endl;
		cout<<" fabs(temp-_ENER_delta_LJ)>1e-12 :: make_judge() "<<endl;
		exit(-1);
	}*/
	if(_BF_flag) {
		_ENER_delta_BF_r=_ENER_BF_eachatom[_INDEX_chosen]  -_ENER_BF_chosenatom_backup_r;
	    _ENER_delta_BF_l=_ENER_BF_eachatom[_INDEX_chosen-1]-_ENER_BF_chosenatom_backup_l;
	}
	if(_AG_flag) {
		_ENER_delta_AG_r=_ENER_AG_eachatom[_INDEX_chosen+1]-_ENER_AG_chosenatom_backup_r;
	    _ENER_delta_AG_c=_ENER_AG_eachatom[_INDEX_chosen]  -_ENER_AG_chosenatom_backup_c;
	    _ENER_delta_AG_l=_ENER_AG_eachatom[_INDEX_chosen-1]-_ENER_AG_chosenatom_backup_l;
	}
	if(_DH_flag) {
		//add your code here;
		_ENER_delta_DH_4=_ENER_DH_eachatom[_INDEX_chosen+2]-_ENER_DH_chosenatom_backup_4;
		_ENER_delta_DH_3=_ENER_DH_eachatom[_INDEX_chosen+1]-_ENER_DH_chosenatom_backup_3;
	    _ENER_delta_DH_2=_ENER_DH_eachatom[_INDEX_chosen]  -_ENER_DH_chosenatom_backup_2;
	    _ENER_delta_DH_1=_ENER_DH_eachatom[_INDEX_chosen-1]-_ENER_DH_chosenatom_backup_1;
	}
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_delta+=_ENER_ele_new-_ENER_ele_old;
	}*/
	_ENER_delta=temp+_ENER_delta_BF_l+_ENER_delta_BF_r
	                +_ENER_delta_AG_l+_ENER_delta_AG_c+_ENER_delta_AG_r
	                +_ENER_delta_DH_1+_ENER_delta_DH_2+_ENER_delta_DH_3+_ENER_delta_DH_4;
	if(_ENER_delta<=0.0) {
		tempnumer_ret=2;//succ!
	} else {
		if( exp(-_ENER_delta*_TEMPERATURE_REP) > rand_seed(iseed_rand) ) {
			tempnumer_ret=2;// succ;
		} else {
			tempnumer_ret=1;// fail;
		}
	}
	//}
}
///////////////////////////
void cmc::make_accept() {
	_ENER_total+=_ENER_delta;
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_ele_old=_ENER_ele_new;
		//cout<<_ENER_ele_old/_ENER_total<<endl;//for test!
	}*/
	int i;
	//#pragma ompparallel for
	for(i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_eachatom.pArray[i][_INDEX_chosen]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_eachatom.pArray[i][_INDEX_chosen]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_eachatom.pArray[i][_INDEX_chosen]=-_DIS_x_eachatom.pArray[_INDEX_chosen][i];
		_DIS_y_eachatom.pArray[i][_INDEX_chosen]=-_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_eachatom.pArray[i][_INDEX_chosen]=-_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}
	/*for(i=_INDEX_chosen+1; i<_SIZE_memo; i++) {
		_ENER_LJ_eachatom.pArray[i][_INDEX_chosen]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_eachatom.pArray[i][_INDEX_chosen]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_eachatom.pArray[i][_INDEX_chosen]=-_DIS_x_eachatom.pArray[_INDEX_chosen][i];//minus;
		_DIS_y_eachatom.pArray[i][_INDEX_chosen]=-_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_eachatom.pArray[i][_INDEX_chosen]=-_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}*/
	if(_FLAG_ener) {
		tempElj[_CINDEXMAP[_INDEX_chn_ind][_INDEX_chn_ind]]+=_ENER_delta_LJ_intra;
		for(i=0;i<_INDEX_chn_ind;i++) {
			tempElj[_CINDEXMAP[_INDEX_chn_ind][i]]+=_ENER_delta_LJ_inter[i];
		}
		for(i=_INDEX_chn_ind+1;i<_NUM_chains;i++) {
			tempElj[_CINDEXMAP[i][_INDEX_chn_ind]]+=_ENER_delta_LJ_inter[i];
		}
		if(_BF_flag) {
			if( _INDEX_chosen > 1 ) {
				tempEbf[_INDEX_chn_atom[_INDEX_chosen-1]]+=_ENER_delta_BF_l;
			}
			tempEbf[_INDEX_chn_ind]+=_ENER_delta_BF_r;
		}
		if(_AG_flag) {
			if( _INDEX_chosen > 1 ) {
				tempEag[_INDEX_chn_atom[_INDEX_chosen-1]]+=_ENER_delta_AG_l;
			}
			tempEag[_INDEX_chn_ind]+=_ENER_delta_AG_c;
			if( _INDEX_chosen < _NUM_atoms ) {
				tempEag[_INDEX_chn_atom[_INDEX_chosen+1]]+=_ENER_delta_AG_r;
			}
		}
		if(_DH_flag) {
			if( _INDEX_chosen > 1 ) {
				tempEdh[_INDEX_chn_atom[_INDEX_chosen-1]]+=_ENER_delta_DH_1;
			}
			tempEdh[_INDEX_chn_ind]+=_ENER_delta_DH_2;
			if( _INDEX_chosen < _NUM_atoms ) {
				tempEdh[_INDEX_chn_atom[_INDEX_chosen+1]]+=_ENER_delta_DH_3;
			}
			if( _INDEX_chosen < ( _NUM_atoms-1 ) ) {
				tempEdh[_INDEX_chn_atom[_INDEX_chosen+2]]+=_ENER_delta_DH_4;
			}
		}
	}
	if(_FLAG_rg2) {
		_COM_x[_INDEX_chn_ind]+=_XX[_INDEX_chosen]-_XX[0];
		_COM_y[_INDEX_chn_ind]+=_YY[_INDEX_chosen]-_YY[0];
		_COM_z[_INDEX_chn_ind]+=_ZZ[_INDEX_chosen]-_ZZ[0];
	}
	//make_mapping_individual();//if just _XX when calc_ener, useless!
	//_MC_NUM_SUC+=1;/*************  for check program ************/
	_MC_NUM_SUC_stat[_INDEX_TEMPERATURE]+=1.0;/*************  for check program ************/
}
///////////////////////////
void cmc::make_reject_all()
{
	int i=0;
	if(_FLAG_ener) {
		_ENER_delta_LJ=0.0;
		_ENER_delta_LJ_intra=0.0;
		for(i=0; i<_NUM_chains; i++) {
			_ENER_delta_LJ_inter[i]=0.0;
		}
		/*_ENER_delta_BF=0.0;
		_ENER_delta_AG=0.0;
		_ENER_delta_DH=0.0;*/
	}
	
	_XX[_INDEX_chosen]=_XX[0];
	_YY[_INDEX_chosen]=_YY[0];
	_ZZ[_INDEX_chosen]=_ZZ[0];
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_ele_new=_ENER_ele_old;
	}*/
	//#pragma ompparallel for
	for(i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][i]=_ENER_LJ_chosenatom_backup[i];
		_DIS2_eachatom.pArray[_INDEX_chosen][i]=_DIS2_chosenatom_backup[i];
		_DIS_x_eachatom.pArray[_INDEX_chosen][i]=_DIS_x_chosenatom_backup[i];
		_DIS_y_eachatom.pArray[_INDEX_chosen][i]=_DIS_y_chosenatom_backup[i];
		_DIS_z_eachatom.pArray[_INDEX_chosen][i]=_DIS_z_chosenatom_backup[i];
	}
	/*for(i=_INDEX_chosen+1; i<; i++) {
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][i]=_ENER_LJ_chosenatom_backup[i];
		_DIS2_eachatom.pArray[_INDEX_chosen][i]=_DIS2_chosenatom_backup[i];
		_DIS_x_eachatom.pArray[_INDEX_chosen][i]=_DIS_x_chosenatom_backup[i];
		_DIS_y_eachatom.pArray[_INDEX_chosen][i]=_DIS_y_chosenatom_backup[i];
		_DIS_z_eachatom.pArray[_INDEX_chosen][i]=_DIS_z_chosenatom_backup[i];
	}*/
	if(_BF_flag) {
		_ENER_BF_eachatom[_INDEX_chosen-1]=_ENER_BF_chosenatom_backup_l;
		_ENER_BF_eachatom[_INDEX_chosen]=  _ENER_BF_chosenatom_backup_r;
	}
	if(_AG_flag) {
		_ENER_AG_eachatom[_INDEX_chosen-1]=_ENER_AG_chosenatom_backup_l;
		_ENER_AG_eachatom[_INDEX_chosen]=  _ENER_AG_chosenatom_backup_c;
		_ENER_AG_eachatom[_INDEX_chosen+1]=_ENER_AG_chosenatom_backup_r;
	}
	if(_DH_flag) {
		_ENER_DH_eachatom[_INDEX_chosen-1]=_ENER_DH_chosenatom_backup_1;
		_ENER_DH_eachatom[_INDEX_chosen]=  _ENER_DH_chosenatom_backup_2;
		_ENER_DH_eachatom[_INDEX_chosen+1]=_ENER_DH_chosenatom_backup_3;
		_ENER_DH_eachatom[_INDEX_chosen+2]=_ENER_DH_chosenatom_backup_4;
	}
	//_MC_NUM_FIL+=1;/*************  for check program ************/
}
///////////////////////////
void cmc::make_reject_only_move() {
	_XX[_INDEX_chosen]=_XX[0];
	_YY[_INDEX_chosen]=_YY[0];
	_ZZ[_INDEX_chosen]=_ZZ[0];
	cout<<" how could this happen?"<<endl;
	//_MC_NUM_FIL+=1;/*************  for check program ************///no need!
}
///////////////////////////
int cmc::make_mcmove(const int INDEX_coor) {
	make_choice(_Runninglist[INDEX_coor]);
	if( make_change() ) {
		make_judge();//in judge _ENER_total can not be changed, should be in accept! 
		if(tempnumer_ret==2) {//succ==2
			//cout<<" succ!! "<<_INDEX_chosen<<endl;
			make_accept();
		/*cout<<" "<<_INDEX_chosen-1<<": "<<_ENER_AG_eachatom[_INDEX_chosen-1]
		    <<" "<<_INDEX_chosen<<": "<<_ENER_AG_eachatom[_INDEX_chosen]
		    <<" "<<_INDEX_chosen+1<<": "<<_ENER_AG_eachatom[_INDEX_chosen+1]<<endl;*/
		} else if(tempnumer_ret==1) {//fail==1
			make_reject_all();
		} else if(tempnumer_ret==-1) {//false==-1, cos' no new energy calc!
			make_reject_only_move();
		}
		//cout<<" _ENER_total="<<_ENER_total<<endl;
		return tempnumer_ret;
	} else {
		//do make sure that there is no change in coordinates,
		//or you must add make_reject_only_move here!
		return -1; //not count, cuz this is the calculation problem, 
		           //actually there is a solution but we did not find it; 
	}
	
}
///////////////////////////
void cmc::make_mapping_wrap() {
	translate2box();
	int i=0;
	int j=0;
	int k=0;
	int numerator=0;
	int sz_res=0;
	int sz_atm=0;
	for(i=0; i<_NUM_chains; i++) {
		sz_res=_system_.chains[i].nresidues;
		for(j=0; j<sz_res; j++) {
			sz_atm=_system_.chains[i].residues[j].natoms;
			for(k=0; k<sz_atm; k++) {
				_system_.allatoms[numerator].x_coordinate=_XX_rec[numerator+1];
				_system_.allatoms[numerator].y_coordinate=_YY_rec[numerator+1];
				_system_.allatoms[numerator].z_coordinate=_ZZ_rec[numerator+1];
				numerator++;	
			}
		}
	}
	if(numerator!=_NUM_atoms) {
		cout<<" Total_Num="<<numerator<<endl;
		cout<<" _NUM_atoms="<<_NUM_atoms<<endl;
		ErrorMSG("Total_Num!=_NUM_atoms :: cmc::make_mapping_wrap");
		exit(SIZEERROR);
	}
}
///////////////////////////
void cmc::make_mapping() {
	int i=0;
	int j=0;
	int k=0;
	int numerator=0;
	int sz_res=0;
	int sz_atm=0;
	for(i=0; i<_NUM_chains; i++) {
		sz_res=_system_.chains[i].nresidues;
		for(j=0; j<sz_res; j++) {
			sz_atm=_system_.chains[i].residues[j].natoms;
			for(k=0; k<sz_atm; k++) {
				_system_.allatoms[numerator].x_coordinate=_XX_rec[numerator+1];
				_system_.allatoms[numerator].y_coordinate=_YY_rec[numerator+1];
				_system_.allatoms[numerator].z_coordinate=_ZZ_rec[numerator+1];
				numerator++;	
			}
		}
	}
	if(numerator!=_NUM_atoms) {
		cout<<" Total_Num="<<numerator<<endl;
		cout<<" _NUM_atoms="<<_NUM_atoms<<endl;
		ErrorMSG("Total_Num!=_NUM_atoms :: cmc::make_mapping");
		exit(SIZEERROR);
	}
}
///////////////////////////
void cmc::fout_conformation(int CHCK_BOND_LEN_OR_NOT) {
	MPI_Gather(_XX, _SIZE_memo, MPI_DOUBLE, _XX_allrep.pArray[0], _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(_YY, _SIZE_memo, MPI_DOUBLE, _YY_allrep.pArray[0], _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(_ZZ, _SIZE_memo, MPI_DOUBLE, _ZZ_allrep.pArray[0], _SIZE_memo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(_PROC_ID!=0) {
		return;
	}
	int i=0;
	int j=0;
	char* FILENAME_conf;
	FILENAME_conf=new char[10];
	for(j=0; j<_NUM_replicas; j++){
		for(i=0; i<_SIZE_memo; i++) {
			_XX_rec[i]=_XX_allrep.pArray[_Index_T_rep[j]][i];
			_YY_rec[i]=_YY_allrep.pArray[_Index_T_rep[j]][i];
			_ZZ_rec[i]=_ZZ_allrep.pArray[_Index_T_rep[j]][i];
			//cout<<" "<<_XX[i]<<" "<<_YY[i]
		}
		make_mapping_wrap();
		if(CHCK_BOND_LEN_OR_NOT==1) {
			chck_bond_len();//after each mapping made;
		}
		sprintf(FILENAME_conf, "confT%03d.pdb", j+1);
		_system_.writepdbinfo(FILENAME_conf, _IFVERBOSE);
	}
	delete[] FILENAME_conf;
}
/*//////////////////////////
inline void cmc::trajectory_rec() {
	
	_OUTfile=fopen( _FNoutput,"ab");    
    //fwrite(&_ENER_total,sizeof(double),1,_OUTfile);
    fwrite(&_XX_rec,sizeof(double),_NUM_rec,_OUTfile);
    fwrite(&_YY_rec,sizeof(double),_NUM_rec,_OUTfile);
    fwrite(&_ZZ_rec,sizeof(double),_NUM_rec,_OUTfile);
    fclose(_OUTfile);
}*/
////////////////////////////
inline void cmc::recording() {
	////// add your code here /////
	sprintf(_RecordingNAME, "pdbfiles/T%03dE%07d.xyz", _INDEX_TEMPERATURE, tempindex_judge);
	FILE* fq;
	fq = fopen(_RecordingNAME,"ab");
	if(fq==NULL) {
		cout<<" can not open file: "<<_RecordingNAME<<endl;
		exit(IOERROR);
	}
	fwrite(&_XX[1],sizeof(double),_SIZE_memo,fq);
	fwrite(&_YY[1],sizeof(double),_SIZE_memo,fq);
	fwrite(&_ZZ[1],sizeof(double),_SIZE_memo,fq);
	/*for(int i=1; i<_SIZE_memo; i++) {
		fprintf(fq, "ATOM  %5d                   %8.3f%8.3f%8.3f\n", i, _XX[i], _YY[i], _ZZ[i]);
	}*/
	/*targetstream<<setw(6)<<setiosflags(ios::left)<<type.c_str()<<resetiosflags(ios::left)
		<<setw(5)<<atomindex<<" "
		<<setw(2)<<atomname_chemical.c_str()
		<<setw(2)<<setiosflags(ios::left)<<atomname_special.c_str()<<resetiosflags(ios::left)<<" "
		<<setw(3)<<residuename.c_str()<<" "
		<<chainname.c_str()
		<<setw(4)<<residueiname<<setw(4)<<" "
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<x_coordinate-double(int(x_coordinate)/1000)*1000.0
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<y_coordinate-double(int(y_coordinate)/1000)*1000.0
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<z_coordinate-double(int(z_coordinate)/1000)*1000.0
		<<setw(26)<<" "<<endl;*/
	fclose(fq);
}
////////////////////////////
void cmc::calc_energy_chosen_atom() {
	int k=0;
	//int l=0;
	//int x=0;
	//int y=0;
	//int z=0;
	//int it_nb_start=_INDEX_atm_in_chn_ind>_NEIGHBOR_lj?(_INDEX_chosen-_NEIGHBOR_lj):(_INDEX_chosen-_INDEX_atm_in_chn_ind);
	//int it_nb_end=(_INDEX_atm_in_chn_ind+1+_NEIGHBOR_lj)<_SIZE_of_chosen_chn?(_INDEX_chosen+_NEIGHBOR_lj+1):(_INDEX_chosen+_SIZE_of_chosen_chn-_INDEX_atm_in_chn_ind);

	//for test	
	/*int it_nb_start_bak=0;
	int it_nb_end_bak=0;
	if( _INDEX_atm_in_chn_ind >= _NEIGHBOR_lj ) {
		it_nb_start_bak=_INDEX_chosen-_NEIGHBOR_lj;
	} else {
		it_nb_start_bak=_INDEX_chosen-_INDEX_atm_in_chn_ind;
	}
	if( _INDEX_atm_in_chn_ind+1+_NEIGHBOR_lj <= _SIZE_of_chosen_chn ) {
		it_nb_end_bak=_INDEX_chosen+_NEIGHBOR_lj+1;
	} else {
		it_nb_end_bak=_INDEX_chosen+_SIZE_of_chosen_chn-_INDEX_atm_in_chn_ind;
	}
	if( it_nb_start!=it_nb_start_bak || it_nb_end!=it_nb_end_bak ) {
		cout<<it_nb_start<<" : "<<it_nb_start_bak<<endl;
		cout<<it_nb_end<<" : "<<it_nb_end_bak<<endl;
	}*/
	//test ends here;
	/*for(k=1;k<_INDEX_chnhead;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[k]-_XX[_INDEX_chosen]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[k]-_YY[_INDEX_chosen]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[k]-_ZZ[_INDEX_chosen]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];				
	}
	for(k=_INDEX_chnhead;k<_INDEX_chosen;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=_XX[k]-_XX[_INDEX_chosen];
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=_YY[k]-_YY[_INDEX_chosen];
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=_ZZ[k]-_ZZ[_INDEX_chosen];
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];				
	}
	for(k=_INDEX_chosen+1;k<=_INDEX_chntail;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=_XX[k]-_XX[_INDEX_chosen];
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=_YY[k]-_YY[_INDEX_chosen];
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=_ZZ[k]-_ZZ[_INDEX_chosen];
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];				
	}
	for(k=_INDEX_chnhead+1;k<_SIZE_memo;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[k]-_XX[_INDEX_chosen]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[k]-_YY[_INDEX_chosen]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[k]-_ZZ[_INDEX_chosen]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];				
	}*/
	//#pragma ompparallel for                           
	for(k=1;k<_INDEX_chosen;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen],_XX[k]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen],_YY[k]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen],_ZZ[k]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];
		/*cout<<_INDEX_chosen<<" "<<k<<" "
		    <<_XX[_INDEX_chosen]<<"-"<<_XX[k]<<"="<<_DIS_x_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_YY[_INDEX_chosen]<<"-"<<_YY[k]<<"="<<_DIS_y_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_ZZ[_INDEX_chosen]<<"-"<<_ZZ[k]<<"="<<_DIS_z_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_DIS2_eachatom.pArray[_INDEX_chosen][k]<<endl;	*/	                               			
	}
	//#pragma ompparallel for
	for(k=_INDEX_chosen+1;k<_SIZE_memo;k++) {//good!! for every update, never miss one coordinate.
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen],_XX[k]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen],_YY[k]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen],_ZZ[k]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                               +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];				
		/*cout<<_INDEX_chosen<<" "<<k<<" "
		    <<_XX[_INDEX_chosen]<<"-"<<_XX[k]<<"="<<_DIS_x_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_YY[_INDEX_chosen]<<"-"<<_YY[k]<<"="<<_DIS_y_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_ZZ[_INDEX_chosen]<<"-"<<_ZZ[k]<<"="<<_DIS_z_eachatom.pArray[_INDEX_chosen][k]<<" "
		    <<_DIS2_eachatom.pArray[_INDEX_chosen][k]<<endl;*/
	}

	if( _INDEX_chn_or_not ) {
		if( ( _SIZE_of_chosen_chn_real > 1 ) && _BF_flag ) {//if chain && with Bond Fluctuation
			if( _TYPE_atom_ind_real == -1 ) { //the start(-1) or the end(1) of the chain
				tempind=_INDEX_chosen+1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_ENER_BF_eachatom[_INDEX_chosen]=Energy_BF(_PARA_KB.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]], 
					_DIS2_eachatom.pArray[_INDEX_chosen][tempind], 
					_BOND_length.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_BOND_length2.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_SIGMAWELL2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
				    _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]);
				//enumer++;
			} else if (_TYPE_atom_ind_real == 1) { 
				tempind=_INDEX_chosen-1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_ENER_BF_eachatom[tempind]=Energy_BF(_PARA_KB.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_DIS2_eachatom.pArray[_INDEX_chosen][tempind],
					_BOND_length.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_BOND_length2.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_SIGMAWELL2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
				    _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]);
		    } else {
				tempind=_INDEX_chosen+1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_ENER_BF_eachatom[_INDEX_chosen]=Energy_BF(_PARA_KB.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_DIS2_eachatom.pArray[_INDEX_chosen][tempind],
					_BOND_length.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_BOND_length2.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_SIGMAWELL2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
				    _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]);
				
				tempind=_INDEX_chosen-1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_ENER_BF_eachatom[tempind]=Energy_BF(_PARA_KB.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_DIS2_eachatom.pArray[_INDEX_chosen][tempind],
					_BOND_length.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_BOND_length2.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
					_SIGMAWELL2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]],
				    _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]);
			}
		}
		//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
		if ( ( _SIZE_of_chosen_chn_real > 2 ) && _AG_flag ) { //if chain and with angle calculation; >3 [3,..]
			if( _TYPE_atom_ind_real == -1 ) {
				tempind=_INDEX_chosen+1; // in ag calculation _dis is theta. and then used as theta-theta0;
				//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				//    fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])>1e-12) {
					_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind+1])
					            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind+1]) 
					            +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind+1]) )
					           /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind+1]) );
					_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0;
					//_DIS2=_DIS*_DIS;
					_ENER_AG_eachatom[tempind]=Energy_AG(
						    (_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+
					    	 _PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0,_DIS);
					//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
				//}
			} else if (_TYPE_atom_ind_real == 1) {
				tempind=_INDEX_chosen-1;
				//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				    //fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])>1e-12) {
					_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind-1])
					            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind-1]) 
					            +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind-1]) )
					           /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind-1]) );
					_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0;
					//_DIS2=_DIS*_DIS;
					_ENER_AG_eachatom[tempind]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+
						_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0,_DIS);
					//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
				//}
			} else { // in the middle of the chain;
			    //... , chosen , ...
				//when <chosen> chosen, important to note that _dis_x and -_dis_x;
				//but in this model this angle calculation is negligible; think!
				
				//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen-1]])>1e-12 &&
				    //fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen+1]])>1e-12) {
				    _DIS=acos( ( ( _DIS_x_eachatom.pArray[_INDEX_chosen][_INDEX_chosen-1]*(-_DIS_x_eachatom.pArray[_INDEX_chosen][_INDEX_chosen+1]) )
					            +( _DIS_y_eachatom.pArray[_INDEX_chosen][_INDEX_chosen-1]*(-_DIS_y_eachatom.pArray[_INDEX_chosen][_INDEX_chosen+1]) )
					            +( _DIS_z_eachatom.pArray[_INDEX_chosen][_INDEX_chosen-1]*(-_DIS_z_eachatom.pArray[_INDEX_chosen][_INDEX_chosen+1]) ) )
					           /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][_INDEX_chosen-1]*_DIS2_eachatom.pArray[_INDEX_chosen][_INDEX_chosen+1]) );
					_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen-1]]+_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen+1]])/2.0;
					//_DIS2=_DIS*_DIS;
					//if( fabs((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0)>1e-6) {
						//cout<<" (theta-theta0)^2="<<_DIS2;}
					//if( fabs((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0)>1e-6) {
					//	cout<<" ka="<<(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0; }
					_ENER_AG_eachatom[_INDEX_chosen]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen-1]]
						+_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[_INDEX_chosen+1]])/2.0,_DIS);
					//cout<<" index "<<_INDEX_chosen<<": "<<_ENER_AG_eachatom[_INDEX_chosen]<<endl;
				//}
				//if( fabs((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0)>1e-6) { 
						//cout<<" index "<<_INDEX_chosen<<": "<<_ENER_AG_eachatom[_INDEX_chosen];}
				if( _SIZE_of_chosen_chn_real>3 ) {
					if(_INDEX_atm_in_chn_ind_real==1) { //0 <1> 2 3 ... n-1;
						tempind=_INDEX_chosen+1;
						//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				   			//fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])>1e-12) {							
							_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind+1])
							            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind+1]) 
						    	        +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind+1]) )
						     	      /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind+1]) );
							_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0;
							//_DIS2=_DIS*_DIS;
							_ENER_AG_eachatom[tempind]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]
								+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0,_DIS);
							//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
						//}
					} else if ( _INDEX_atm_in_chn_ind_real==(_SIZE_of_chosen_chn_real-2) ) { //0 1 ... n-3 <n-2> n-1;
						tempind=_INDEX_chosen-1;
						//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				   			//fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])>1e-12) {	
							_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind-1])
							            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind-1]) 
							            +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind-1]) )
							           /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind-1]) );
							
							_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0;
							//_DIS2=_DIS*_DIS;
							_ENER_AG_eachatom[tempind]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]
								+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0,_DIS);
							//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
						//}
					} else {// not head not tail, not <1> not <n-2>, what else? when 4: 0[h] 1<1> 2<4-2> 3[t].
					  // the only result: size of chain >=5 and in the middle of a chain not tail not head not <1> not <n-2>
						tempind=_INDEX_chosen+1;
						//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				   			//fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])>1e-12) {	
							_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind+1])
							            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind+1]) 
						    	        +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind+1]) )
						     	      /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind+1]) );
							
							_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0;
							//_DIS2=_DIS*_DIS;
							_ENER_AG_eachatom[tempind]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]
								+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind+1]])/2.0,_DIS);
							//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
							//
						//}
						tempind=_INDEX_chosen-1;
						//if( fabs(_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]])>1e-12 &&
				   			//fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])>1e-12) {	
							_DIS=acos( ( (_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[tempind][tempind-1])
							            +(_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[tempind][tempind-1]) 
							            +(_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[tempind][tempind-1]) )
							           /sqrt(_DIS2_eachatom.pArray[_INDEX_chosen][tempind]*_DIS2_eachatom.pArray[tempind][tempind-1]) );
							_DIS=_DIS-(_THETAZERO.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0;
							//_DIS2=_DIS*_DIS;
							_ENER_AG_eachatom[tempind]=Energy_AG((_PARA_KA.pArray[_Res_chosen][_INDEX_RES_ATM[tempind]]
								+_PARA_KA.pArray[_INDEX_RES_ATM[tempind]][_INDEX_RES_ATM[tempind-1]])/2.0,_DIS);
							//cout<<" index "<<tempind<<": "<<_ENER_AG_eachatom[tempind]<<endl;
						//}
					}
				}
			}
			//cout<<endl;
		}
		if ( ( _SIZE_of_chosen_chn_real > 3 ) && _DH_flag ) {//if chain and with dihedral calculation;	
			//cout<<" index_chosen="<<_INDEX_chosen<<endl;
			if( _TYPE_atom_ind_real == -1 ) { //first...
				tempind=_INDEX_chosen+2; //left dihedral
				dihedral2H1(_INDEX_chosen);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				////getchar();
				//}
			} else if( _INDEX_atm_in_chn_ind_real==1 ) { //second... 0 <1> 2 3 ... 
				tempind=_INDEX_chosen+1; //left dihedral
				dihedral2H2(_INDEX_chosen-1);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				//}
				if (_SIZE_of_chosen_chn_real > 4) {
					tempind=_INDEX_chosen+2; //left dihedral
					dihedral2H1(_INDEX_chosen);
					_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
					//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				}
				////getchar();
			} else if ( _INDEX_atm_in_chn_ind_real==(_SIZE_of_chosen_chn_real-2) ) { //0 1 ... n-3 <n-2> n-1;
				tempind=_INDEX_chosen; //left dihedral
				dihedral2H3(_INDEX_chosen-2);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				//}
				if (_SIZE_of_chosen_chn_real > 4) {
					tempind=_INDEX_chosen-1; //left dihedral
					dihedral2H4(_INDEX_chosen-3);
					_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
					//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				}
				//getchar();
			} else if (_TYPE_atom_ind_real == 1) {
				tempind=_INDEX_chosen-1;
				dihedral2H4(_INDEX_chosen-3);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				//getchar();
				//}
			} else { // in the middle of the chain and chanlen > 4
			    //... , chosen , ...
				//when <chosen> chosen, important to note that _dis_x and - _dis_x;
				//but in this model this angle calculation is negligible; think!
				tempind=_INDEX_chosen;
			    dihedral2H3(_INDEX_chosen-2);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
				//}
				
				tempind=_INDEX_chosen+1;
				dihedral2H2(_INDEX_chosen-1);
				_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
				//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;

				if( _SIZE_of_chosen_chn_real>5 ) {
					if(_INDEX_atm_in_chn_ind_real==2) { //0 1 <2> 3 ... n-1;
						tempind=_INDEX_chosen+2;
						dihedral2H1(_INDEX_chosen);
						_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
						//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
						//}
					} else if ( _INDEX_atm_in_chn_ind_real==(_SIZE_of_chosen_chn_real-3) ) { //0 1 ... <n-3> n-2 n-1;
						tempind=_INDEX_chosen-1;
						dihedral2H4(_INDEX_chosen-3);
						_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
						//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
						//}
					} else {// not head not tail, not <2> not <n-3>, what else? when 6: 0[h] 1 <2> 3 <4> 5[t].
					  // the only result: size of chain >=5 and in the middle of a chain not tail not head not <1> not <n-2>
						tempind=_INDEX_chosen+2;
						dihedral2H1(_INDEX_chosen);
						_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
						//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;

						tempind=_INDEX_chosen-1;
						dihedral2H4(_INDEX_chosen-3);
						_ENER_DH_eachatom[tempind]=Energy_DH(_tempkd, _tempDH-_tempphi0);
						//cout<<" index "<<tempind<<": "<<_ENER_DH_eachatom[tempind]<<endl;
					}
				}
				//getchar();
			}
			//cout<<endl;
		}
	} else { // if not chain, L-J potential calculation between "neighboring atoms";
		//if( _SIZE_of_chosen_chn > 1 ) { // 2 atoms or more;

		////#pragma ompparallel for
		for(k=_Index_lneighbor_ind_real; k<_INDEX_chosen; k++) {
			//cout<<" here1, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
			_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
			/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
				for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
					for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
						_DIS_x=_PBL_X*x;
						_DIS_y=_PBL_Y*y;
						_DIS_z=_PBL_Z*z;//(a+da)^2+(b+db)^2+(c+dc)^2=a^2+b^2+c^2+.....
						_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
						                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
						                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
						                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
						_DIS6=TRIPLE(_DIS2);
						_DIS12=_DIS6*_DIS6; 
						if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
							_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=
								Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
									      _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
									      _DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
						}
					}
				}
			}*/
			
			_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
			//_DIS6=_DIS2*_DIS2*_DIS2;
			//_DIS12=_DIS6*_DIS6; 
			if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
				if(_ACCINDEX>0) { 
					index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
					_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
				} else {
					_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=
					Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							  _LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							  _EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							  _SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						      _SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						      _SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						      _SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						      _DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
				}
			}
		}
		////#pragma ompparallel for
		for(k=_INDEX_chosen+1; k<_Index_rneighbor_ind_real; k++) {
			//cout<<" here2, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
			_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
			/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
				for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
					for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
						_DIS_x=_PBL_X*x;
						_DIS_y=_PBL_Y*y;
						_DIS_z=_PBL_Z*z;
						_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
						                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
						                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
						                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
						_DIS6=TRIPLE(_DIS2);
						_DIS12=_DIS6*_DIS6;
						if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
							_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
							_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
						}
					}
				}
			}*/
						
			_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
			//_DIS6=_DIS2*_DIS2*_DIS2;
			//_DIS12=_DIS6*_DIS6;
			if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
				if(_ACCINDEX>0) { 
					index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
					_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
				} else {
					_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
					_SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
				}
			}
		}
		//}
	}

    //Lennard-Jones Potential Part;  inter chain;
	//cout<<" [k:";
	//#pragma ompparallel for
	for(k=1; k<_INDEX_chnhead_real; k++) {// must be '<', attention! 
		//cout<<" here3, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;

					if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
					}
				}
			}
		}*/

		_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		//_DIS6=_DIS2*_DIS2*_DIS2;
		//_DIS12=_DIS6*_DIS6;

		if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
			if(_ACCINDEX>0) { 
				index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
			} else {
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
			}
		}
	}
	// intra chain but not neighbor
	//#pragma ompparallel for
	for(k=_INDEX_chnhead_real; k<_Index_lneighbor_ind_real; k++) {// must be '<', attention! 
		//cout<<" here3, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;

					if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
					}
				}
			}
		}*/
		/*if(_INDEX_chn_or_not) {
			_DIS_x=_XX[k]-_XX[_INDEX_chosen];
			_DIS_y=_YY[k]-_YY[_INDEX_chosen];
			_DIS_z=_ZZ[k]-_ZZ[_INDEX_chosen];
			_DIS2=SQUARE(_DIS_x)+SQUARE(_DIS_y)+SQUARE(_DIS_z);
		} else {
			_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		}*/
		_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		//_DIS6=_DIS2*_DIS2*_DIS2;
		//_DIS12=_DIS6*_DIS6;

		if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
			if(_ACCINDEX>0) { 
				index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
			} else {
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
			}
		}
	}
	//cout<<" ("<<_INDEX_chosen<<")";
	// intra chain but not neighbor
	//#pragma ompparallel for
	for(k=_Index_rneighbor_ind_real; k<_INDEX_chntail_real; k++) {// must be '<', attention! 
		//cout<<" here4, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;
					if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
					}
				}
			}
		}*/
		
		/*if(_INDEX_chn_or_not) {
			_DIS_x=_XX[k]-_XX[_INDEX_chosen];
			_DIS_y=_YY[k]-_YY[_INDEX_chosen];
			_DIS_z=_ZZ[k]-_ZZ[_INDEX_chosen];
			_DIS2=SQUARE(_DIS_x)+SQUARE(_DIS_y)+SQUARE(_DIS_z);
		} else {
			_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		}*/
		_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		//_DIS6=_DIS2*_DIS2*_DIS2;
		//_DIS12=_DIS6*_DIS6;
		if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
			if(_ACCINDEX>0) { 
				index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
			} else {
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
			}
		}
	}

	//inter chain
	//#pragma ompparallel for
	for(k=_INDEX_chntail_real; k<_SIZE_memo; k++) {// must be '<', attention! 
		//cout<<" here4, _INDEX_chosen="<<_INDEX_chosen<<"; k="<<k<<endl;
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;
					if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
						_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
					}
				}
			}
		}*/

		_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k];
		//_DIS6=_DIS2*_DIS2*_DIS2;
		//_DIS12=_DIS6*_DIS6;
		if(_DIS2<_R_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]]) {
			if(_ACCINDEX>0) { 
				index_lj=int(_DIS2/_LJ_interval.pArray[_Res_chosen][_INDEX_RES_ATM[k]]);
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=_ELJ_ACC.pArray[_Res_chosen][_INDEX_RES_ATM[k]][index_lj];
			} else {
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_PPTYPE_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_LAMBDA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_EPSILON_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA12_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA6_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA2_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],
				_SIGMA_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]],_DIS2)-_E_cut_RR_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[k]];
			}
		}
	}
	//cout<<"]"<<endl;
	//cout<<" _ENER_eachatom["<<_INDEX_chosen<<"]="<<_ENER_eachatom[_INDEX_chosen]<<endl;
}
/////////////////////////
////////////////////////////
double cmc::calc_energy() {                       //the total energy;
	tell_procid(); cout<<" calculating the total energy start (calc_energy)"<<endl;
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	//int x=0;
	//int y=0;
	//int z=0;

	_ENER_total=0.0;

	int Size_ATM=0;
	int Size_ATM_real=0;
	int Size_ATM_A=0;

	int index_atm=0;
	int index_atm_a=0; 

	/*for(i=1;i<_SIZE_memo-1;i++) { // important, you have to do this before calculation of eachatom: <angle calculation>
		for(j=i+1;j<_SIZE_memo;j++) {
			if( fabs( _DIS_x_eachatom.pArray[i][j]-DIS_PBC_X(_XX[i],_XX[j]) ) >1e-6 
			||  fabs( _DIS_y_eachatom.pArray[i][j]-DIS_PBC_Y(_YY[i],_YY[j]) ) >1e-6 
			||	fabs( _DIS_z_eachatom.pArray[i][j]-DIS_PBC_Z(_ZZ[i],_ZZ[j]) ) >1e-6 ) {
				cout<<"i="<<i<<", j="<<j<<" x "<<_DIS_x_eachatom.pArray[i][j]<<" "<<_DIS_x_eachatom.pArray[j][i]<<" "<<DIS_PBC_X(_XX[i],_XX[j])<<endl;
				cout<<"i="<<i<<", j="<<j<<" y "<<_DIS_y_eachatom.pArray[i][j]<<" "<<_DIS_y_eachatom.pArray[j][i]<<" "<<DIS_PBC_Y(_YY[i],_YY[j])<<endl;
				cout<<"i="<<i<<", j="<<j<<" z "<<_DIS_z_eachatom.pArray[i][j]<<" "<<_DIS_z_eachatom.pArray[j][i]<<" "<<DIS_PBC_Z(_ZZ[i],_ZZ[j])<<endl;
			}
			_DIS_x_eachatom.pArray[i][j]=DIS_PBC_X(_XX[i],_XX[j]);
			_DIS_y_eachatom.pArray[i][j]=DIS_PBC_Y(_YY[i],_YY[j]);
			_DIS_z_eachatom.pArray[i][j]=DIS_PBC_Z(_ZZ[i],_ZZ[j]);
			_DIS2_eachatom.pArray[i][j]=_DIS_x_eachatom.pArray[i][j]*_DIS_x_eachatom.pArray[i][j]
					                   +_DIS_y_eachatom.pArray[i][j]*_DIS_y_eachatom.pArray[i][j]
					                   +_DIS_z_eachatom.pArray[i][j]*_DIS_z_eachatom.pArray[i][j];	
			_DIS_x_eachatom.pArray[j][i]=-_DIS_x_eachatom.pArray[i][j];
			_DIS_y_eachatom.pArray[j][i]=-_DIS_y_eachatom.pArray[i][j];
			_DIS_z_eachatom.pArray[j][i]=-_DIS_z_eachatom.pArray[i][j];
			_DIS2_eachatom.pArray[j][i]=_DIS2_eachatom.pArray[i][j];
		}
	}*/

	for(i=0; i<_NUM_chains; i++) {// intra chain neighbor (hook<->chain) (L.J.<->others);
		Size_ATM=_SIZE_of_chn[i];
		Size_ATM_real=_SIZE_of_chn_real[i];
		if( _INDEX_chn_or_not_chain[i] ) {//if chain
			if( ( Size_ATM_real > 1 ) && _BF_flag ) {
				for(j=1; j<=Size_ATM; j++) {// must be '<', attention! 
				//index+1, cos' _XX _YY _ZZ defined to do this! 
					tempind_x=index_atm+j;
					tempind_y=index_atm+j+1;
					if(tempind_x>=_INDEX_CHN_TAIL_real[i]) { break; }
					//cout<<tempind_x<<":"<<tempind_y<<endl;
					_ENER_total+=Energy_BF(_PARA_KB.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
						_DIS2_eachatom.pArray[tempind_x][tempind_y],
						_BOND_length.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
						_BOND_length2.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
						_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
						_SIGMAWELL2_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
					    _SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]);
					//enumer++;
					//cout<<index_atm+j<<"<->"<<index_atm+j+1<<endl;
					//cout<<tempind_x<<" "<<sqrt(_DIS2_eachatom.pArray[tempind_x][tempind_y])<<endl;
					//cout<<_ENER_total<<endl;
				}
			}
			if( ( Size_ATM_real > 2 ) && _AG_flag ) {
				for(j=1; j<=Size_ATM; j++) {
					tempind_x=index_atm+j;
					tempind_y=index_atm+j+1;//this one is the angle location;
					if(tempind_y>=_INDEX_CHN_TAIL_real[i]) { break; }
					//_dis is actually theta;
					_DIS=acos( ( (_DIS_x_eachatom.pArray[tempind_x][tempind_y]*_DIS_x_eachatom.pArray[tempind_y][tempind_y+1])
				           		+(_DIS_y_eachatom.pArray[tempind_x][tempind_y]*_DIS_y_eachatom.pArray[tempind_y][tempind_y+1]) 
				     	        +(_DIS_z_eachatom.pArray[tempind_x][tempind_y]*_DIS_z_eachatom.pArray[tempind_y][tempind_y+1]) )
				           	   /sqrt(_DIS2_eachatom.pArray[tempind_x][tempind_y]*_DIS2_eachatom.pArray[tempind_y][tempind_y+1]) );
					_DIS=_DIS-(_THETAZERO.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]+_THETAZERO.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])/2.0;
					//_DIS2=_DIS*_DIS;
					//if( fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]])>1e-12 &&
					//    fabs(_PARA_KA.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])>1e-12 ) {
					_ENER_total+=Energy_AG(
						(_PARA_KA.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]
						+_PARA_KA.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])/2.0,
						_DIS);
					/*cout<<tempind_x<<" "<<Energy_AG(
						(_PARA_KA.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]
						+_PARA_KA.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])/2.0,
						_DIS)<<endl;*/
						////getchar();
					//}
					/*cout<<" c="<<i<<" a="<<tempind_y<<" EA="<<Energy_AG( 
						(_PARA_KA.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]
						+_PARA_KA.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])/2.0,
						_DIS2)<<endl;*/
					//cout<<" tt> tempind_y:"<<tempind_y<<" ener="<<Energy_AG((_PARA_KA.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]] 
					//+_PARA_KA.pArray[_INDEX_RES_ATM[tempind_y]][_INDEX_RES_ATM[tempind_y+1]])/2.0,_DIS2)<<endl;
				}
			}
			//cout<<0<<" "<<2<<" "<<_ENER_DH_eachatom[2]<<endl;
			if( ( Size_ATM_real > 3 ) && _DH_flag ) {
				for(j=1; j<=Size_ATM; j++) {
					tempind_x=index_atm+j; // x, y, y+1, y+2;
					//tempind_y=index_atm+j+1;// y+1 is the angle location;
					if(tempind_x>=_INDEX_CHN_TAIL_real[i]-2) { break; }
					dihedral2H1(tempind_x);
					_ENER_total+=Energy_DH(_tempkd, _tempDH-_tempphi0);
					//cout<<tempind_x<<" "<<tempind_x+2<<" "<<_ENER_DH_eachatom[tempind_x+2]<<" "<<Energy_DH(_tempkd, _tempDH-_tempphi0)<<endl;
					//}
					
				}
			}
			//cout<<++tempind_x<<" "<<tempind_x+2<<" "<<_ENER_DH_eachatom[tempind_x+2]<<endl;
			//cout<<++tempind_x<<" "<<tempind_x+2<<" "<<_ENER_DH_eachatom[tempind_x+2]<<endl;
		} else {// if not chain 
			for(j=1; j<=Size_ATM-_NEIGHBOR_lj; j++) {//index+1, cos' _XX _YY _ZZ defined to do this! 
				for(k=1; k<=_NEIGHBOR_lj; k++) {
					tempind_x=index_atm+j;
					tempind_y=index_atm+j+k;
					/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
						for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
							for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
								_DIS_x=_PBL_X*x;
								_DIS_y=_PBL_Y*y;
								_DIS_z=_PBL_Z*z;
								_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
								                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
								_DIS6=TRIPLE(_DIS2);
								_DIS12=_DIS6*_DIS6;
								if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
									_ENER_total+=Energy_LJ(_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
										                   _SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
										                   _SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
										                   _DIS12,
										                   _DIS6)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
								}
							}
						}
					}*/
					
					_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y];
					//_DIS6=_DIS2*_DIS2*_DIS2;
					//_DIS12=_DIS6*_DIS6;
					if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
						if(_ACCINDEX>0) { 
							index_lj=int(_DIS2/_LJ_interval.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]);
							_ENER_total+=_ELJ_ACC.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]][index_lj];
						} else {
							_ENER_total+=Energy_LJ(_PPTYPE_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_LAMBDA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							                   _SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							                   _SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							                   _SIGMA2_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							                   _SIGMA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							                   _DIS2)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
						}
					}
					//enumer++;
				}
			}
		}
		index_atm+=Size_ATM;
	}
	//cout<<" check point: "<<_ENER_total<<endl;
	index_atm=0;
	for(i=0; i<_NUM_chains; i++) {// intra chain, but not neighbor, ( LJ - all );
		Size_ATM=_SIZE_of_chn[i];
		for(j=1; j<Size_ATM-_NEIGHBOR_lj; j++) {// must be '<', attention! 
            //index+1, cos' _XX _YY _ZZ defined to do this! 
			for(k=j+_NEIGHBOR_lj+1; k<=Size_ATM; k++) {// must be '<', attention! 
			//index+1, cos' _XX _YY _ZZ defined to do this! 
				/*tempind_x=index_atm+j;
				tempind_y=index_atm+k;
				for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
					for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
						for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
							//here add code for PME etc. which will take PBC into the consideration;
							_DIS_x=_PBL_X*x;
							_DIS_y=_PBL_Y*y;
							_DIS_z=_PBL_Z*z;
							_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
							                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
							                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
							                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
							_DIS6=TRIPLE(_DIS2); 
							_DIS12=_DIS6*_DIS6;
							if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
								_ENER_total+=Energy_LJ(_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
							}
						}
					}
				}*/
				tempind_x=index_atm+j;
				tempind_y=index_atm+k;
				
				/*if(_INDEX_chn_or_not) {
					_DIS_x=_XX[tempind_y]-_XX[tempind_x];
					_DIS_y=_YY[tempind_y]-_YY[tempind_x];
					_DIS_z=_ZZ[tempind_y]-_ZZ[tempind_x];
					_DIS2=SQUARE(_DIS_x)+SQUARE(_DIS_y)+SQUARE(_DIS_z);
				} else {
					_DIS2=_DIS2_eachatom.pArray[tempind_y][tempind_x];
				}*/
				//mcout<<tempind_x<<":"<<tempind_y<<" "<<sqrt(_DIS2)<<" v.s. "<<sqrt(_DIS2_eachatom.pArray[tempind_x][tempind_y])<<endl;
				_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y];
				//_DIS6=_DIS2*_DIS2*_DIS2; 
				//_DIS12=_DIS6*_DIS6;
				if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
						if(_ACCINDEX>0) { 
							index_lj=int(_DIS2/_LJ_interval.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]);
							_ENER_total+=_ELJ_ACC.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]][index_lj];
						} else {
							_ENER_total+=Energy_LJ(_PPTYPE_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_LAMBDA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_SIGMA2_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
							_SIGMA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],_DIS2)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
					}
				}
				//enumer++;
				//cout<<index_atm+j<<"<->"<<index_atm+k<<":"<<endl;
				//printf("%20.15f\n", _ENER_total);
			}
		}
		index_atm+=Size_ATM;
	}

	index_atm=0;
	for(i=0; i<_NUM_chains-1; i++) {// inter chain, ( LJ - all );
		Size_ATM=_SIZE_of_chn[i];
		index_atm_a=0;
		for(j=0; j<i+1; j++) {
			index_atm_a+=_SIZE_of_chn[j];
		}
		for(j=i+1; j<_NUM_chains; j++) {// must be '<', attention! 
			Size_ATM_A=_SIZE_of_chn[j];
			for(k=1; k<=Size_ATM; k++) {// must be '<', attention! 
			//index+1, cos' _XX _YY _ZZ defined to do this!  	
				tempind_x=index_atm+k;
				for(l=1; l<=Size_ATM_A; l++) {//index+1, cos' _XX _YY _ZZ defined to do this!
					tempind_y=index_atm_a+l;
					/*for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
						for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
							for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
								_DIS_x=_PBL_X*x;
								_DIS_y=_PBL_Y*y;
								_DIS_z=_PBL_Z*z;
								_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
								                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
								_DIS6=TRIPLE(_DIS2);
								_DIS12=_DIS6*_DIS6;
								if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
									_ENER_total+=Energy_LJ(_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
									_SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
									_SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
									_DIS12,_DIS6)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
								}
							}
						}
					}*/	
					if( fabs(tempind_x-tempind_y)==_NEIGHBOR_lj && ( _TYPE_atom_real[tempind_x]==0 || _TYPE_atom_real[tempind_y]==0 ) ) {
						continue;
					} else {
						_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y];
						//_DIS6=_DIS2*_DIS2*_DIS2;
						//_DIS12=_DIS6*_DIS6;
						if(_DIS2<_R_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]) {
							if(_ACCINDEX>0) { 
								index_lj=int(_DIS2/_LJ_interval.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]]);
								_ENER_total+=_ELJ_ACC.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]][index_lj];
							} else {
								_ENER_total+=Energy_LJ(_PPTYPE_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_LAMBDA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_EPSILON_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA12_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA6_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA2_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],
								_SIGMA_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]],_DIS2)-_E_cut_RR_eachres.pArray[_INDEX_RES_ATM[tempind_x]][_INDEX_RES_ATM[tempind_y]];
							}
						}
					}		
					//enumer++;
				}
			}
			index_atm_a+=Size_ATM_A;
		}
		index_atm+=Size_ATM;
	}
	tell_procid(); cout<<" calculating the total energy done! (calc_energy) E="<<_ENER_total<<endl;
	return _ENER_total;
}
////////////////////////////////
void cmc::init_energy() {
    // bfflag_bak for calculating neighboring _DIS2_eachatom;
	/*double bfflag_bak=false;
	if(_BF_flag==false) {
		cout<<" temporarily change _BF_flag to: [true] (neighbor distance calculation)"<<endl;
		bfflag_bak=true;
		_BF_flag=true;
		//_PARA_K=0.0;
	}*/

	tell_procid(); cout<<" *--> energy initializing start .. "<<endl;	
	///////////////   code start  //////////////////
	//enumer=0;
	calc_energy();//necessary for angle calculation since i-j x <-> j-i -x
	//chck_bond_len();
	//cout<<" calc_energy() succ!"<<endl;
	//cout<<" --------------------- enumer = "<<enumer<<endl;
	//enumer=0;
	calc_energy_each_atom();
	//chck_bond_len();
    double Total=get_energy();//necessary
	//cout<<" --------------------- enumer = "<<enumer<<endl;
	//cout<<" calc_energy_each_atom() succ!"<<endl;
	/* this is for the wall attraction */
	///////////////////// this is just for check at initial point!
	if( fabs(_ENER_total-Total) > 1e-6 ) {
		printf(" _ENER_total-total_ener=%20.15f\n", _ENER_total-Total);
		printf(" _ENER_total=%20.15f\n", _ENER_total);
		printf(" _ENER_get=%20.15f\n", Total);
		ErrorMSG("fabs(_ENER_total-total_ener) > 1e-6 ::cmc::init_energy");
		exit(LOGICERROR);
	}
    ///////////////   code start  //////////////////
    tell_procid(); cout<<" *--> energy initializing done! "<<endl;
	/*if(bfflag_bak==true) {
		cout<<" change _BF_flag back to: [false] (neighbor distance calculation)"<<endl;
		_BF_flag=false;
	}*/
	MPI_Barrier(MPI_COMM_WORLD);
}
/////////////////////////////////
void cmc::calc_energy_each_atom() {
	for(int i=1; i<_SIZE_memo; i++) {
		make_choice(i);
		calc_energy_chosen_atom();
		//cout<<" EA["<<i<<"]="<<_ENER_AG_eachatom[i]<<endl;
	}
	_Enery_initialized=true;
}
/////////////////////////////////
double cmc::get_energy() {
	double Ret_Ener=0.0;
	int i=0;
	int j=0;
	double Ret_Ener_LJ=0.0;
	for(i=1; i<_SIZE_memo-1; i++) {
		for(j=i+1; j<_SIZE_memo; j++) {
			Ret_Ener_LJ+=_ENER_LJ_eachatom.pArray[i][j];
		}
	}
	Ret_Ener+=Ret_Ener_LJ;
	tell_procid(); printf(" Ret_Ener_LJ=%.8f",Ret_Ener_LJ); cout<<" (get_energy) "<<endl;

	double Ret_Ener_BF=0.0;
	if(_BF_flag) {
		for(i=1; i<_SIZE_memo; i++) {
			Ret_Ener_BF+=_ENER_BF_eachatom[i];
		}
	}
	Ret_Ener+=Ret_Ener_BF;
	tell_procid(); printf(" Ret_Ener_BF=%.8f",Ret_Ener_BF); cout<<" (get_energy) "<<endl;

	double Ret_Ener_AG=0.0;
	if(_AG_flag) {
		for(i=1; i<_SIZE_memo; i++) {
			Ret_Ener_AG+=_ENER_AG_eachatom[i];
			//cout<<" each> index:"<<i<<" ener="<<_ENER_AG_eachatom[i]<<endl;
		}
	}
	Ret_Ener+=Ret_Ener_AG;
	tell_procid(); printf(" Ret_Ener_AG=%.8f",Ret_Ener_AG); cout<<" (get_energy) "<<endl;

	double Ret_Ener_DH=0.0;
	if(_DH_flag) {
		for(i=1; i<_SIZE_memo; i++) {
			Ret_Ener_DH+=_ENER_DH_eachatom[i];
			//cout<<" each> index:"<<i<<" ener="<<_ENER_DH_eachatom[i]<<endl;
		}
	}
	Ret_Ener+=Ret_Ener_DH;
	tell_procid(); printf(" Ret_Ener_DH=%.8f",Ret_Ener_DH); cout<<" (get_energy) "<<endl;
	
	////////////////////// this is just for check at initial point!
	return Ret_Ener;
}
void cmc::init_mpi(int pargc, char** pargv) {
	MPI_Init(&pargc, &pargv); //fixed mpi_initialization;
	MPI_Comm_rank(MPI_COMM_WORLD, &_PROC_ID); //_PROC_ID = identity of process;
	MPI_Comm_size(MPI_COMM_WORLD, &_PROC_size); //number of processes in mpi;
}
void cmc::term_mpi() {
	tell_procid(); cout<<" mpi terminated!"<<endl;
	MPI_Finalize();
}
//////
void cmc::initialization() {
	load_parameters( "_config.in" );//single process.
	load_conformation(); //done, only _system_ of p1 initialized.
	write_parameters( "_config.in" ); //done, as above
	broadcast_parameters(); //done.
	memo_allocation(); //done.
	memo_setzero(); //done.
	memo_evaluation();
    init_energy();
    chck_bond_len();//secure!!
    srand(_PROC_ID*(unsigned)time(NULL));
    run_with_stepnumber(_NUM_atoms*1000); //minimization
    init_statistic();
}
////////////////////////////////////
///////////////////////////////////
/*void cmc::cout_parameters() {
	cout<<setiosflags(ios::left)<<endl;
	//cout<<setw(30)<<" ::_system_.natoms: "<<_system_.natoms<<endl;
	cout<<setw(30)<<" ::ee: "<<ee<<endl;
	cout<<setw(30)<<" ::_NUM_atoms: "<<_NUM_atoms<<endl;
	cout<<setw(30)<<" ::_totalener: "<<_ENER_total<<endl;
	cout<<setw(30)<<" ::_IFVERBOSE: "<<_IFVERBOSE<<endl;
	//cout<<setw(30)<<" ::_START_fromzero: "<<_START_fromzero<<endl;
	cout<<resetiosflags(ios::left)<<endl;
	//}
}*/
////////////////////////////////////
void cmc::check_energy() {
	double E_temp;
	double E_temp_2;
	E_temp=get_energy();
	E_temp_2=_ENER_total;
	calc_energy();
	printf(" T[%8.4f]: %20.15f(calc)-%20.15f(get)(%20.15f(delta))=%e\n", _TEMPERATURE, _ENER_total, E_temp, E_temp_2, _ENER_total-E_temp);//before fout_conf_PBC
	if( fabs(_ENER_total-E_temp)>1e-6) {
		cout<<" there's something wrong with energy, mpi termed!"<<endl;
		exit(LOGICERROR);
	}

}
////////////////////////////////////
void cmc::run_with_stepnumber(const int stepnumber) {

	tell_procid(); cout<<endl<<" starting minimization run..... "<<endl<<endl;

	for(int i=0; i<stepnumber; i++) {
		//cout<<i<<endl;
		make_mcmove( int(rand_seed(iseed_index)*_realSize) );
		//check_energy();
	}
	tell_procid(); cout<<endl<<" minimization ended......"<<endl<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
}
/////////////////
void cmc::run_rem() {
	//tell_procid(); cout<<" passed1"<<endl;
	MPI_Gather(&_ENER_total, 1, MPI_DOUBLE, _E_rep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//tell_procid(); cout<<" passed2"<<endl;
    if( _PROC_ID == 0 ) {
		int i=0;
		//double TTemp=0.0;
		int index_i=0;
		int index_j=0;
		for(i=_NUM_replicas-1; i>0; i--) {		//pair(i, i-1)
			index_i=_Index_T_rep[i]; // t[i] in which rep
			index_j=_Index_T_rep[i-1];
			_B_delta=(_T_rep_eachrep[i]-_T_rep_eachrep[i-1])*(_E_rep[index_j]-_E_rep[index_i]);
			//+1.5_realSize*(_E_rep[index_j]-_E_rep[index_i])
			//cout<<"_B_delta="<<_B_delta<<endl;
			_NUM_rem[i-1]+=1;
			if( _B_delta>0.0 ) {	
				if(rand_seed(iseed_rand)>=exp(-_B_delta)) {	   //rejected
					continue;		//do nothing
				}
			}
			//exchange these two conformation pair(i, i-1)
		
			/*TTemp=_T_rep_eachrep[index_i]; 
			_T_rep_eachrep[index_i]=_T_rep_eachrep[index_j];	
			_T_rep_eachrep[index_j]=TTemp;*/

			_Index_T_rep[i]=index_j;	
			_Index_T_rep[i-1]=index_i;
			
			_Index_rep_T[index_i]=i-1; 
			_Index_rep_T[index_j]=i;	

			_NUM_rem_acc[i-1]+=1;
		}
		/*tell_id();cout<<" T="<<1.0/_T_rep_eachreplica[0]
			      <<" index_of_T="<<_Index_REP_eachreplica[0]
				  <<" index_of_R="<<_Index_T_eachreplica[_Index_REP_eachreplica[0]]<<endl;*/
	}
	MPI_Scatter(_Index_rep_T, 1, MPI_INT, &_INDEX_TEMPERATURE, 1, MPI_INT, 0, MPI_COMM_WORLD);
	_TEMPERATURE_REP=_T_rep_eachrep[_INDEX_TEMPERATURE];
	_TEMPERATURE=_T_eachrep[_INDEX_TEMPERATURE];
	//tell_procid(); cout<<" T["<<_INDEX_TEMPERATURE<<"]="<<_TEMPERATURE<<"; "<<endl;//T_rep="<<_TEMPERATURE_REP<<";"<<endl;
	//////////////////
	//prep_afterREM();//not important any more if temperature have been exchanged;
    //////////////////
	//ener_dispatch_r1();
	//_MC_NUM_SUC+=1;/*************  for check program ************/
	//_MC_NUM_TOT+=1;/*************  for check program ************/
}
/////////////////
inline bool cmc::ener_dispatch() {
	if( (_ENER_total < _E_lowest) || ( _ENER_total >= _E_highest ) ) {
		return false;
	}
	tempindex_judge=int((_ENER_total-_E_lowest)/_E_interval);
	_Probability.pArray[_INDEX_TEMPERATURE][tempindex_judge]+=1;
	return true;
	//_NUM_conf_eachtemp[_INDEX_TEMPERATURE]+=1;	
}
/////////////////
/////////////////
void cmc::run() {
	//int i=0;
	//int TEN_NUM_atoms=0;
	tell_procid(); cout<<" before the production run..."<<endl;
	//cout<<dihedral2(1.0, 0.0, 0.0, 0.0, 1.0, 0.0 ,1.0,0.0,0.0);
	//getchar();
	check_energy();
	rand_init(_PROC_ID);
	srand(_PROC_ID*(unsigned)time(NULL));

	cout<<" starting production run..... "<<endl<<endl;

	//_MC_NUM_TOT=0;
	//_MC_NUM_SUC=0;
	//_MC_NUM_FIL=0;

	MEMOSETZERO(_MC_NUM_SUC_stat, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_TOT_stat, sizeof(double)*_NUM_replicas);
	MEMOSETZERO(_MC_NUM_FIL_stat, sizeof(double)*_NUM_replicas);

	int i=0;
	if( _PROC_ID == 0 ) {
		bool file_exist=false;
		//cout<<" your running option is: "<<_RUN_options<<endl;
		ifstream check_file("accept_ratio.dat");
		if( check_file ) {
			file_exist=true;
		}
		check_file.close();
		if(file_exist) {
			//system("del accept_ratio.dat");
			if(system("rm accept_ratio.dat")){};
			cout<<" later you can check: [ accept_ratio.dat ] for rem accepted ratio ~ "<<endl;
		}
	}
	for(_I_eachstep=0, _I_totalnum=0; _I_totalnum < _RUNTIMES_totalnum; _I_eachstep++) {
		//cout<<_I_eachstep<<endl;
		if( _I_eachstep == _RUNTIMES_eachstep ) {
			//tell_procid(); cout<<_I_totalnum<<endl;
			if( ( (_I_totalnum+1)%_RUNTIMES_remgap ) == 0 ) {
				MPI_Barrier(MPI_COMM_WORLD); // for collecting data!!!
				run_rem();
				/*if(_I_totalnum*_RUNTIMES_eachstep>1e5) {
					rand_update(_PROC_ID);
				}*/
			}
			if( (_I_totalnum+1)%_RUNTIMES_output == 0 ) {
				tell_procid(); 
				cout<<" _I_totalnum/_RUNTIMES_totalnum="<<setw(6)<<_I_totalnum+1
					<<"/"<<setw(6)<<setiosflags(ios::left)<<_RUNTIMES_totalnum
					<<resetiosflags(ios::left)<<endl;
				//fout_conformation(false);
				if(_PROC_ID==0) {
					ofstream ofstream_acptrto("accept_ratio.dat", ios::app);
					if( ofstream_acptrto == NULL ) {
						ErrorMSG("can not open [ accept_ratio.dat ] :: cmc::run()");
					}
					ofstream_acptrto<<setw(4)<<setiosflags(ios::left)
									<<(_I_totalnum+1)/_RUNTIMES_output
									<<" output step:"<<resetiosflags(ios::left)<<endl;
					for(i=0; i<_NUM_replicas-1; i++) {
						ofstream_acptrto<<" accepted ratio: ["
							 <<setw(8)<<setiosflags(ios::left)
							 <<setprecision(5)<<setiosflags(ios::fixed)
							 <<_T_eachrep[i]<<", "
							 <<_T_eachrep[i+1]<<"]:  ";
						if (_NUM_rem[i]!=0) {
							 ofstream_acptrto<<double(_NUM_rem_acc[i])/double(_NUM_rem[i])
							 <<resetiosflags(ios::left)<<endl;
						} else {
							ofstream_acptrto<<0.0<<resetiosflags(ios::left)<<" -> when output, no rem yet."<<endl;
						}
					}
					ofstream_acptrto.close();
				}
			}
			_I_eachstep=0;
			_I_totalnum++;
			
			//trajectory_rec();			
		}
		//int tempdix=int(rand_seed(iseed_index)*_realSize);
		//if ( tempdix==63 )cout<<tempdix<<endl;
		if( make_mcmove(int( rand_seed(iseed_index)*_realSize) )>0  ) {
			if( ener_dispatch() ) {
				statistic();
				if(_RUNTIMES_recording!=0) {
					if ( (_I_eachstep+1)/_RUNTIMES_recording==0 ) {
						recording();
					}
				}
			}
		}
	}
	fout_conformation(1);
	MPI_Barrier(MPI_COMM_WORLD);
	check_energy();
	MPI_Barrier(MPI_COMM_WORLD);
	//tell_procid(); cout<<" succ rate : "<<double(_MC_NUM_SUC)/double(_MC_NUM_TOT)<<endl;
	//tell_procid(); cout<<" fail rate : "<<double(_MC_NUM_FIL)/double(_MC_NUM_TOT)<<endl;
	_I_totalnum=0;
	//MPI_Barrier(MPI_COMM_WORLD);
	getrealrate();
		//fout_conformation_eachrep(true);
	//}
}
void cmc::getrealrate() {
	//
	MPI_Reduce(_MC_NUM_SUC_stat, _MC_NUM_SUC_all, _NUM_replicas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_MC_NUM_TOT_stat, _MC_NUM_TOT_all, _NUM_replicas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_MC_NUM_FIL_stat, _MC_NUM_FIL_all, _NUM_replicas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(_PROC_ID!=0) {
		return;
	}
	cout<<endl;
	cout<<" this is the probability of the accepted trials if the move is not prohibited."<<endl;
	cout<<endl;
	cout<<endl;
	int i=0;
	for(i=0; i<_NUM_replicas; i++) {
		//printf(" @proc[%03d]", i); 
		if(i<10) {
			cout<<" @proc[00"<<i<<"]";
		} else if(i<100) {
			cout<<" @proc[0"<<i<<"]";
		} else {
			cout<<" @proc["<<i<<"]";
		}
		cout<<" succ rate : "<<_MC_NUM_SUC_all[i]/_MC_NUM_TOT_all[i]<<endl;
	}
	cout<<endl;
	cout<<" this is the probability of the prohibited trials."<<endl;
	cout<<endl;
	cout<<endl;
	for(i=0; i<_NUM_replicas; i++) {
		//printf(" @proc[%03d]", i); 
		if(i<10) {
			cout<<" @proc[00"<<i<<"]";
		} else if(i<100) {
			cout<<" @proc[0"<<i<<"]";
		} else {
			cout<<" @proc["<<i<<"]";
		}
		cout<<" fail rate : "<<_MC_NUM_FIL_all[i]/(_MC_NUM_TOT_all[i]+_MC_NUM_FIL_all[i])<<endl;
	}
}
/////////////////////////////////
void cmc::ensembler() {
	MPI_Reduce(_Probability.pArray[0], _Probability_tot.pArray[0], _NUM_replicas*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(_PROC_ID!=0) {
		return;
	}
	int i=0;
	int rs=0;
	char* FileName;
	FileName=new char[30];
	sprintf(FileName, "probability_all.dat");
	ofstream wham_stream_all(FileName);
	if(wham_stream_all==NULL) {
		cout<<" can not write into file: [ "<<FileName<<" ];"<<endl;
		exit(IOERROR);
	}
	//double tempprob=0.0;
	
	for(i=0; i<_E_totalnum; i++) {
		//tempprob=0.0;
		for(rs=0; rs<_NUM_replicas; rs++) {
			_Probability_all[i]+=_Probability_tot.pArray[rs][i];
		}
		wham_stream_all<<setiosflags(ios::left)
					   <<setw(10)<<_E_lowest+_E_interval*i<<" "
					   <<setw(20)<<_Probability_all[i]
					   <<resetiosflags(ios::left)<<endl;
		
	}
	cout<<" energy distribution info was written into: [ "<<FileName<<" ];"<<endl;
	wham_stream_all.close();
	
	ofstream wham_stream_each( "probability_each.dat" );
	if(wham_stream_each==NULL) {
		cout<<" can not write into file: [ probability_each.dat ];"<<endl;
		exit(IOERROR);
	}
	for(i=0; i<_E_totalnum; i++) {
		wham_stream_each<<setiosflags(ios::left)<<setw(10)<<_E_lowest+_E_interval*i;
		for(rs=0; rs<_NUM_replicas; rs++) {
			wham_stream_each<<" "<<setw(10)<<_Probability_tot.pArray[rs][i];	
		}
		wham_stream_each<<resetiosflags(ios::left)<<endl;
	}
	cout<<" energy distribution info was written into: [ probability_each.dat ];"<<endl;
	wham_stream_each.close();

	ofstream ensembler_stream( "ensembler.gpl" );
	if(ensembler_stream==NULL) {
		cout<<" can not write into file: [ ensembler.gpl ];"<<endl;
		exit(IOERROR);
	}
	ensembler_stream<<" nrep="<<_NUM_replicas<<endl;
	ensembler_stream<<" emin="<<_E_lowest<<endl;
	ensembler_stream<<" emax="<<_E_highest<<endl;
	ensembler_stream.close();
	if(system("gnuplot draw_prob.gpl")){};
	cout<<" u may check file: [ probability.eps ] now."<<endl;
	delete[] FileName;
	
}
/////////////////////////////////
void cmc::init_statistic() {
	if(!_Enery_initialized) {// no need any more... but put it here, I dont wanna change;
		cout<<" Energy not initialized, error: init_statistic() "<<endl;
		exit(LOGICERROR);
	}
	//tcharn=new char[30];
	//if(_FLAG_ener) {
		//_rh2_stream=new ofstream[_NUM_chains];
	//}
	//if(_FLAG_rg2) {
		//_rg2_stream=new ofstream[_NUM_chains];
	//}
	cout<<" !! initializating variables ..." <<endl; 
	double tempEtot=0.0;
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		stat_head=_INDEX_CHN_HEAD[stat_i];
		stat_tail=_INDEX_CHN_TAIL[stat_i]+1;
		if(_FLAG_ener) {
			for(stat_j=0; stat_j<=stat_i; stat_j++) {
				tempElj[_CINDEXMAP[stat_i][stat_j]]=0.0;
				if(stat_i==stat_j) {		
					stat_head_j=_INDEX_CHN_HEAD[stat_j];
					stat_tail_j=_INDEX_CHN_TAIL[stat_j]+1;
					for(stat_k=stat_head_j; stat_k<stat_tail_j; stat_k++) {
						for(stat_l=stat_k+1; stat_l<stat_tail_j; stat_l++) {
							tempElj[_CINDEXMAP[stat_i][stat_j]]+=_ENER_LJ_eachatom.pArray[stat_k][stat_l];
						}
					}
				} else {
					for(stat_k=_INDEX_CHN_HEAD[stat_i]; stat_k<=_INDEX_CHN_TAIL[stat_i]; stat_k++) {
						for(stat_l=_INDEX_CHN_HEAD[stat_j]; stat_l<=_INDEX_CHN_TAIL[stat_j]; stat_l++) {
							tempElj[_CINDEXMAP[stat_i][stat_j]]+=_ENER_LJ_eachatom.pArray[stat_k][stat_l];							
						}
					}
				}
				tempEtot+=tempElj[_CINDEXMAP[stat_i][stat_j]];
			}
			if(_BF_flag) {
				////add your code here;
				tempEbf[stat_i]=0.0;
				for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
					tempEbf[stat_i]+=_ENER_BF_eachatom[stat_j];
				}
				tempEtot+=tempEbf[stat_i];
			}
			if(_AG_flag) {
				////add your code here;
				tempEag[stat_i]=0.0;
				for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
					tempEag[stat_i]+=_ENER_AG_eachatom[stat_j];
				}
				/*if(stat_i==0 && fabs(tempEag[stat_i])>1e-12) {
					cout<<" chn0 Eag="<<tempEag[stat_i]<<endl;
					//getchar();
				}*/
				tempEtot+=tempEag[stat_i];
			}
			if(_DH_flag) {
				////add your code here;
				tempEdh[stat_i]=0.0;
				for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
					tempEdh[stat_i]+=_ENER_DH_eachatom[stat_j];
				}
				tempEtot+=tempEdh[stat_i];
			}
		}
		if(_FLAG_rg2) {
			//sprintf(tcharn,"rg2%03d%03d%05d.dat",(stat_i)+1,_NUM_chains,_NUM_atoms);
			//_rg2_stream[stat_i].open(string(tcharn).c_str());
			//cout<<" file: ["<<tcharn<<"] opened."<<endl;
			//initialization;
			_COM_x[stat_i]=0.0;
			_COM_y[stat_i]=0.0;
			_COM_z[stat_i]=0.0;
			for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
				_COM_x[stat_i]+=_XX[stat_j];
				_COM_y[stat_i]+=_YY[stat_j];
				_COM_z[stat_i]+=_ZZ[stat_j];
			}
			stat_com_x[stat_i]=_COM_x[stat_i]/stat_size;
			stat_com_y[stat_i]=_COM_y[stat_i]/stat_size;
			stat_com_z[stat_i]=_COM_z[stat_i]/stat_size;
			cout<<" chain "<<stat_i<<" com: "<<endl;
			cout<<"       "<<stat_com_x[stat_i]<<endl;
			cout<<"       "<<stat_com_y[stat_i]<<endl;
			cout<<"       "<<stat_com_z[stat_i]<<endl;
			if(_FLAG_dis) {
				for(stat_j=0; stat_j<stat_i; stat_j++) {
					_DISSTAT[_CINDEXMAP[stat_i][stat_j]]=DIS_PBC_X(stat_com_x[stat_i]-stat_com_x[stat_j])*DIS_PBC_X(stat_com_x[stat_i]-stat_com_x[stat_j])
						+DIS_PBC_Y(stat_com_y[stat_i]-stat_com_y[stat_j])*DIS_PBC_Y(stat_com_y[stat_i]-stat_com_y[stat_j])
						+DIS_PBC_Z(stat_com_z[stat_i]-stat_com_z[stat_j])*DIS_PBC_Z(stat_com_z[stat_i]-stat_com_z[stat_j]);
				}
			}

			//over
			if( _SIZE_of_chn[stat_i] > 1 ) {//necessary?
				_RG2_x[stat_i]=0.0;
				_RG2_y[stat_i]=0.0;
				_RG2_z[stat_i]=0.0;
				_RG2_ec[stat_i]=0.0;
				for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
					tempnum=_XX[stat_j]-stat_com_x[stat_i];
					_RG2_x[stat_i]+=tempnum*tempnum;
					//_RG2_x_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
					tempnum=_YY[stat_j]-stat_com_y[stat_i];
					_RG2_y[stat_i]+=tempnum*tempnum;
					//_RG2_y_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
					tempnum=_ZZ[stat_j]-stat_com_z[stat_i];
					_RG2_z[stat_i]+=tempnum*tempnum;
					//_RG2_z_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
				}
				_RG2_ec[stat_i]=(_RG2_x[stat_i]+_RG2_y[stat_i]+_RG2_z[stat_i])/stat_size;
				cout<<" chain "<<stat_i<<" RG2: "<<_RG2_ec[stat_i]<<endl;
				cout<<"       rg2x="<<_RG2_x[stat_i]/stat_size<<endl;
				cout<<"       rg2y="<<_RG2_y[stat_i]/stat_size<<endl;
				cout<<"       rg2z="<<_RG2_z[stat_i]/stat_size<<endl;
				if(_RG_totalnum>1) {
					_indexRGSTAT[stat_i]=int((_RG2_ec[stat_i]-_RG_lowest)/_RG_interval);
					if(_indexRGSTAT[stat_i]>=_RG_totalnum) {
						_indexRGSTAT[stat_i]=_RG_totalnum-1;
					} else if (_indexRGSTAT[stat_i]<0) {
						_indexRGSTAT[stat_i]=0;
					} else {
						_indexRGSTAT[stat_i]=_indexRGSTAT[stat_i];
					}
					cout<<"       "<<_indexRGSTAT[stat_i]<<endl;
				}
			}
		}
	}
	if( _FLAG_ener && _len_EINT>0 ) {
		_CNSTAT=0.0;
		tempindexstat=_EINTlist[0];
		stat_head=_INDEX_CHN_HEAD[tempindexstat];
		stat_tail=_INDEX_CHN_TAIL[tempindexstat]+1;
		for(stat_k=stat_head,stat_i=0; stat_k<stat_tail; stat_k++,stat_i++) {
			tempindexstat_ij=_EINTlist[1];
			stat_head_j=_INDEX_CHN_HEAD[tempindexstat_ij];
			stat_tail_j=_INDEX_CHN_TAIL[tempindexstat_ij]+1;
			for(stat_l=stat_head_j,stat_j=0; stat_l<stat_tail_j; stat_l++,stat_j++) {
				if(_DIS2_eachatom.pArray[stat_k][stat_l]<_SIGMADIS_eachres.pArray[_INDEX_RES_ATM[stat_k]][_INDEX_RES_ATM[stat_l]]) {
					_CN_backup.pArray[stat_i][stat_j]=1.0;
					_CNSTAT+=1.0;
				}
				if(_IFVERBOSE) cout<<" "<<_CN_backup.pArray[stat_i][stat_j];
			}
			if(_IFVERBOSE) cout<<endl;
		}
		//cout<<" _CNSTAT="<<_CNSTAT<<endl;
		//getchar();
	}
	cout<<" @Proc"<<_PROC_ID<<" Etot="<<tempEtot<<" v.s. CurrentE="<<_ENER_total<<endl;
	if(_FLAG_ener) {
		if( fabs(tempEtot-_ENER_total)>1e-6 || _ENER_total!=_ENER_total || tempEtot!=tempEtot ) {
			cout<<" energy wrong..."<<endl;
			exit(-1);
		}
	}
	cout<<" !! statistic initialized, all variables have to be evaluated, no doubt."<<endl<<endl;
	_Statistic_over=false;
}
////////////////
inline void cmc::statistic() { //every step!!
	//int i=0;
	//int j=0;
	//tell_procid(); cout<<" statistic() "<<endl;	
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		stat_head=_INDEX_CHN_HEAD[stat_i];
		stat_tail=_INDEX_CHN_TAIL[stat_i]+1;
		if(_FLAG_ener) {
			for(stat_j=0; stat_j<=stat_i; stat_j++) {
				_ENERLJ_stat.pArray[_CINDEXMAP[stat_i][stat_j]][tempindex_judge]+=tempElj[_CINDEXMAP[stat_i][stat_j]];
			}
			if(_BF_flag) {
				////add your code here;
				_ENERBF_stat.pArray[stat_i][tempindex_judge]+=tempEbf[stat_i];
				//tell_procid(); cout<<" statistic() BF over! "<<stat_i<<endl;
			}
			if(_AG_flag) {
				////add your code here;
				_ENERAG_stat.pArray[stat_i][tempindex_judge]+=tempEag[stat_i];
				//tell_procid(); cout<<" statistic() AG over! "<<stat_i<<endl;
			}
			if(_DH_flag) {
				////add your code here;
				_ENERDH_stat.pArray[stat_i][tempindex_judge]+=tempEdh[stat_i];
				//tell_procid(); cout<<" statistic() DH over! "<<stat_i<<endl;	
			}
		}

		if(_FLAG_rg2) {
			//cout<<_COM_x[stat_i]<<endl;
			if(tempnumer_ret>1 && stat_i==_INDEX_chn_ind) { //only current chn com and rg2 needs to be updated.
				stat_com_x[stat_i]=_COM_x[stat_i]/stat_size;
				stat_com_y[stat_i]=_COM_y[stat_i]/stat_size;
				stat_com_z[stat_i]=_COM_z[stat_i]/stat_size;
				//cout<<_COM_x[stat_i]/stat_size<<endl;
				if( _SIZE_of_chosen_chn>1 ) { //necessary?
					_RG2_x[stat_i]=0.0;
					_RG2_y[stat_i]=0.0;
					_RG2_z[stat_i]=0.0;
					for(stat_j=stat_head; stat_j<stat_tail; stat_j++) {
						tempnum=_XX[stat_j]-stat_com_x[stat_i];
						_RG2_x[stat_i]+=tempnum*tempnum;
						//_RG2_x_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
						tempnum=_YY[stat_j]-stat_com_y[stat_i];
						_RG2_y[stat_i]+=tempnum*tempnum;
						//_RG2_y_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
						tempnum=_ZZ[stat_j]-stat_com_z[stat_i];
						_RG2_z[stat_i]+=tempnum*tempnum;
						//_RG2_z_stat.pArray[stat_i][tempindex_judge]+=tempnum*tempnum;
					}

					if( _RG_totalnum>1 ) {
						_indexRGSTAT[stat_i]=int(((_RG2_x[stat_i]+_RG2_y[stat_i]+_RG2_z[stat_i])/stat_size-_RG_lowest)/_RG_interval);
						if(_indexRGSTAT[stat_i]>=_RG_totalnum) {
							_indexRGSTAT[stat_i]=_RG_totalnum-1;
						} /*else if (_indexRGSTAT[stat_i]<0) {//normally impossible since rg from 0.0.
							_indexRGSTAT[stat_i]=0;
						} *//*else {//unchanged index
							_indexRGSTAT[stat_i]=_indexRGSTAT[stat_i];
						}*/
						//cout<<"       "<<_indexRGSTAT[stat_i]<<endl;
					}
				}
				if(_FLAG_dis) {
					for(stat_j=0; stat_j<stat_i; stat_j++) { // att: i!=j; otherwise it's an error!!! this index not in array _cindexmap;
						tempnumerstat_x=DIS_PBC_X(stat_com_x[stat_i]-stat_com_x[stat_j]);
						tempnumerstat_y=DIS_PBC_Y(stat_com_y[stat_i]-stat_com_y[stat_j]);
						tempnumerstat_z=DIS_PBC_Z(stat_com_z[stat_i]-stat_com_z[stat_j]);
						_DISSTAT[_CINDEXMAP[stat_i][stat_j]]=tempnumerstat_x*tempnumerstat_x+tempnumerstat_y*tempnumerstat_y+tempnumerstat_z*tempnumerstat_z;
						//tell_procid(); cout<<" statistic() DIS over! "<<stat_i<<" "<<stat_j<<endl;	
					}	
					for(stat_j=stat_i+1; stat_j<_NUM_chains; stat_j++) {
						tempnumerstat_x=DIS_PBC_X(stat_com_x[stat_i]-stat_com_x[stat_j]);
						tempnumerstat_y=DIS_PBC_Y(stat_com_y[stat_i]-stat_com_y[stat_j]);
						tempnumerstat_z=DIS_PBC_Z(stat_com_z[stat_i]-stat_com_z[stat_j]);
						_DISSTAT[_CINDEXMAP[stat_i][stat_j]]=tempnumerstat_x*tempnumerstat_x+tempnumerstat_y*tempnumerstat_y+tempnumerstat_z*tempnumerstat_z;
						//tell_procid(); cout<<" statistic() DIS over! "<<stat_i<<" "<<stat_j<<endl;	
					}
				}
			}

			//if( stat_i==1 && fabs(stat_com_x[stat_i]-50)>1e-12 )
			/*if(tempindex_judge==4063 && stat_i==1) {
				cout<<" Eid["<<tempindex_judge<<"] Cid["<<stat_i<<"] com= "<<stat_com_x[stat_i]<<" rg2= "<<_RG2_x[stat_i]<<endl;
				cout<<" @PROC"<<_PROC_ID<<" "<<_Probability.pArray[_INDEX_TEMPERATURE][tempindex_judge]<<endl;
			}*/
			_COM_x_stat.pArray[stat_i][tempindex_judge]+=Coor_PBC_X(stat_com_x[stat_i]);
			_COM_y_stat.pArray[stat_i][tempindex_judge]+=Coor_PBC_Y(stat_com_y[stat_i]);
			_COM_z_stat.pArray[stat_i][tempindex_judge]+=Coor_PBC_Z(stat_com_z[stat_i]);

			//tell_procid(); cout<<" statistic() COM over! "<<stat_i<<endl;	

			if( _SIZE_of_chn[stat_i]>1 ) {//necessary?
				/*i=int((_RG2_x[stat_i]/stat_size-_RG_lowest)/_RG_interval);
				if(i>=_RG_totalnum) {
					_RG2_x_stat.pArray[stat_i][tempindex_judge][_RG_totalnum-1]+=1.0;
				} else if (i<0) {
					_RG2_x_stat.pArray[stat_i][tempindex_judge][0]+=1.0;
				} else {
					_RG2_x_stat.pArray[stat_i][tempindex_judge][i]+=1.0;
				}
				i=int((_RG2_y[stat_i]/stat_size-_RG_lowest)/_RG_interval);
				if(i>=_RG_totalnum) {
					_RG2_y_stat.pArray[stat_i][tempindex_judge][_RG_totalnum-1]+=1.0;
				} else if (i<0) {
					_RG2_y_stat.pArray[stat_i][tempindex_judge][0]+=1.0;
				} else {
					_RG2_y_stat.pArray[stat_i][tempindex_judge][i]+=1.0;
				}
				i=int((_RG2_z[stat_i]/stat_size-_RG_lowest)/_RG_interval);
				if(i>=_RG_totalnum) {
					_RG2_z_stat.pArray[stat_i][tempindex_judge][_RG_totalnum-1]+=1.0;
				} else if (i<0) {
					_RG2_z_stat.pArray[stat_i][tempindex_judge][0]+=1.0;
				} else {
					_RG2_z_stat.pArray[stat_i][tempindex_judge][i]+=1.0;
				}*/
				
				if(_RG_totalnum>1 ) {
					_RG2_stat.pArray[stat_i][tempindex_judge][_indexRGSTAT[stat_i]]+=1.0;
				}
				_RG2_actual_x.pArray[stat_i][tempindex_judge]+=_RG2_x[stat_i];
				_RG2_actual_y.pArray[stat_i][tempindex_judge]+=_RG2_y[stat_i];
				_RG2_actual_z.pArray[stat_i][tempindex_judge]+=_RG2_z[stat_i];
			}
			//if(_COM_x_stat.pArray[stat_i][tempindex_judge]>1e-6) cout<<stat_i<<" "<<_COM_x_stat.pArray[stat_i][tempindex_judge]<<" "<<_RG2_x_stat.pArray[stat_i][tempindex_judge]<<endl;
		}
	}
	if( _NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) { 
		for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
			for(stat_j=0; stat_j<stat_i; stat_j++) {
				tempindexstat_ij=_CINDEXMAP[stat_i][stat_j];
				tempnumerstat_x=_DISSTAT[tempindexstat_ij];
				_DIS_stat_actual.pArray[tempindexstat_ij][tempindex_judge]+=tempnumerstat_x;
				tempindexstat=int((tempnumerstat_x-_DISSTAT_lowest)/_DISSTAT_interval);
				if(tempindexstat>=_RG_totalnum) {
					_DIS_stat.pArray[tempindexstat_ij][tempindex_judge][_RG_totalnum-1]+=1.0;
				} /*else if (tempindexstat<0) {
					_DIS_stat.pArray[_CINDEXMAP[stat_i][stat_j]][tempindex_judge][0]+=1.0;
				} */else {
					_DIS_stat.pArray[tempindexstat_ij][tempindex_judge][tempindexstat]+=1.0;
				}
			}	
		}
	}

	//tell_procid(); cout<<" here1"<<endl;
	if( _FLAG_ener && _len_EINT>0 ) { 
		if(_EINT_totalnum>1) {
			tempindexstat=int((tempElj[_EINT_index]-_EINT_lowest)/_EINT_interval);
			if(tempindexstat>=_EINT_totalnum) {
				_EINT_stat.pArray[tempindex_judge][_EINT_totalnum-1]+=1.0;
			} else if(tempindexstat<0) {
				_EINT_stat.pArray[tempindex_judge][0]+=1.0;
			} else {
				_EINT_stat.pArray[tempindex_judge][tempindexstat]+=1.0;
			}
		}
		if( tempnumer_ret==2 ) { //only succ we update the CN;
			if(_INDEX_chn_ind==_EINTlist[0]) {
				tempindexstat=_EINTlist[1];
				stat_head=_INDEX_CHN_HEAD[tempindexstat];
				stat_tail=_INDEX_CHN_TAIL[tempindexstat]+1;
				for(stat_k=stat_head,stat_j=0; stat_k<stat_tail; stat_k++,stat_j++) {
					/*cout<<endl<<"_INDEX_chosen:"<<_INDEX_chosen<<endl;
					cout<<"stat_k:"<<stat_k<<endl;
					cout<<"_Res_chosen:"<<_Res_chosen<<endl;
					cout<<"_INDEX_RES_ATM[stat_k]:"<<_INDEX_RES_ATM[stat_k]<<endl;
					cout<<"dis2:"<<_DIS2_eachatom.pArray[_INDEX_chosen][stat_k]<<endl;
					cout<<"sigd:"<<_SIGMADIS_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[stat_k]]<<endl;
					cout<<"_CN_backup0:"<<_CN_backup.pArray[_INDEX_chosen][stat_k]<<endl;
					cout<<"_CNSTAT0:"<<_CNSTAT<<endl;*/
					if(_DIS2_eachatom.pArray[_INDEX_chosen][stat_k]<_SIGMADIS_eachres.pArray[_Res_chosen][_INDEX_RES_ATM[stat_k]]) {
						if(_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]<1.0) {
							_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]=1.0;
							_CNSTAT+=1.0;
						}
					} else {
						if(_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]>1e-6) {
							_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]=0.0;
							_CNSTAT-=1.0;
						}
					}
					//cout<<"_CN_backup1:"<<_CN_backup.pArray[_INDEX_chosen][stat_k]<<endl;
					//cout<<"_CNSTAT1:"<<_CNSTAT<<endl;
				}
			} else if(_INDEX_chn_ind==_EINTlist[1]) {
				//cout<<" in this case this is impossible!"<<endl;
				tempindexstat=_EINTlist[0];
				stat_head=_INDEX_CHN_HEAD[tempindexstat];
				stat_tail=_INDEX_CHN_TAIL[tempindexstat]+1;
				for(stat_k=stat_head,stat_j=0; stat_k<stat_tail; stat_k++,stat_j++) {
					if(_DIS2_eachatom.pArray[stat_k][_INDEX_chosen]<_SIGMADIS_eachres.pArray[_INDEX_RES_ATM[stat_k]][_Res_chosen]) {
						if(_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]<1.0) {
							_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]=1.0;
							_CNSTAT+=1.0;
						}
					} else {
						if(_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]>1e-6) {
							_CN_backup.pArray[_INDEX_atm_in_chn_ind][stat_j]=0.0;
							_CNSTAT-=1.0;
						}
					}
				}
			}
		}
		_CN_stat[tempindex_judge]+=_CNSTAT;
	}

	//tell_procid(); cout<<" here2"<<endl;
	/*double tempEtot=0.0;
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) { 
		//stat_size=_SIZE_of_chn[stat_i];
		//stat_head=_INDEX_CHN_HEAD[stat_i];
		//stat_tail=_INDEX_CHN_TAIL[stat_i];
		if(_FLAG_ener) {
			for(stat_j=0; stat_j<=stat_i; stat_j++) {
				cout<<" inter["<<stat_i<<"]["<<stat_j<<"]="<<tempElj[stat_i][stat_j]<<endl;
				tempEtot+=tempElj[stat_i][stat_j];
			}
		}
	}
	if( fabs(tempEtot-_ENER_total)>1e-6) {
		cout<<" Etot="<<tempEtot<<" v.s. CurrentE="<<_ENER_total<<endl;
		cout<<" energy wrong..."<<endl;
		exit(-1);
	}*/

}
////////////////
void cmc::output_statistic() { //every _runtimes_eachstep
	MPI_Reduce(_ENERLJ_stat.pArray[0], _ENERLJ_stat_tot.pArray[0], (_NUM_chains+1)*_NUM_chains/2*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_ENERAG_stat.pArray[0], _ENERAG_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_ENERBF_stat.pArray[0], _ENERBF_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_ENERDH_stat.pArray[0], _ENERDH_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_COM_x_stat.pArray[0], _COM_x_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_COM_y_stat.pArray[0], _COM_y_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_COM_z_stat.pArray[0], _COM_z_stat_tot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(_RG2_x_stat.pArray[0][0], _RG2_x_stat_tot.pArray[0][0], _NUM_chains*_E_totalnum*_RG_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(_RG2_y_stat.pArray[0][0], _RG2_y_stat_tot.pArray[0][0], _NUM_chains*_E_totalnum*_RG_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(_RG2_z_stat.pArray[0][0], _RG2_z_stat_tot.pArray[0][0], _NUM_chains*_E_totalnum*_RG_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_RG2_stat.pArray[0][0], _RG2_stat_tot.pArray[0][0], _NUM_chains*_E_totalnum*_RG_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if( _len_EINT>0 && _EINT_totalnum>1 ) {
		MPI_Reduce(_EINT_stat.pArray[0], _EINT_stat_tot.pArray[0], _E_totalnum*_EINT_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	if( _len_EINT>0 ) {
		MPI_Reduce(_CN_stat, _CN_stat_tot, _E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	MPI_Reduce(_RG2_actual_x.pArray[0], _RG2_actual_xtot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_RG2_actual_y.pArray[0], _RG2_actual_ytot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(_RG2_actual_z.pArray[0], _RG2_actual_ztot.pArray[0], _NUM_chains*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if( _NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) {
		MPI_Reduce(_DIS_stat.pArray[0][0], _DIS_stat_tot.pArray[0][0], (_NUM_chains-1)*_NUM_chains/2*_E_totalnum*_RG_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(_DIS_stat_actual.pArray[0], _DIS_stat_actual_tot.pArray[0], (_NUM_chains-1)*_NUM_chains/2*_E_totalnum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	if(_PROC_ID!=0) {
		return;
	}
	cout<<" start statistic output!! "<<endl;
	char tempfilename[50];
	//char totofilename[30];
	//ofstream  totoosELJ; // for check.
	ofstream* temposELJ;
	ofstream* temposEAG;
	ofstream* temposEBF;
	ofstream* temposEDH;
	ofstream* temposCOM;
	ofstream* temposRG2;
	//ofstream* temposRG2_plot;
	//ofstream* temposRG2_x_plot;
	//ofstream* temposRG2_y_plot;
	//ofstream* temposRG2_z_plot;
	ofstream* temposRG2_plot;
	
	temposELJ=new ofstream[(_NUM_chains+1)*_NUM_chains/2];
	temposEBF=new ofstream[_NUM_chains];
	temposEAG=new ofstream[_NUM_chains];
	temposEDH=new ofstream[_NUM_chains];
	temposCOM=new ofstream[_NUM_chains];
	temposRG2=new ofstream[_NUM_chains];
	//temposRG2_x_plot=new ofstream[_NUM_chains];
	//temposRG2_y_plot=new ofstream[_NUM_chains];
	//temposRG2_z_plot=new ofstream[_NUM_chains];
	temposRG2_plot=new ofstream[_NUM_chains];
	
	//double tempNum_x=0.0;
	//double tempNum_y=0.0;
	//double tempNum_z=0.0;
	//double tempSum_x=0.0;
	//double tempSum_y=0.0;
	//double tempSum_z=0.0;
	double tempNum_all=0.0;
	double tempSum_all=0.0;
	if(_NUM_chains>1 && _FLAG_rg2 && _FLAG_dis ) {
		ofstream* temposDIS;
		ofstream* temposDIS_plot;
		temposDIS=new ofstream[(_NUM_chains-1)*_NUM_chains/2];
		temposDIS_plot=new ofstream[(_NUM_chains-1)*_NUM_chains/2];
		for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
			for(stat_j=0; stat_j<stat_i; stat_j++) {
				sprintf(tempfilename, "DIS%03d%03d.dat", stat_i+1, stat_j+1);
				cout<<" now writing : "<<tempfilename<<" to stream ~ "<<_CINDEXMAP[stat_i][stat_j]<<" ;-) "<<endl;
				temposDIS[_CINDEXMAP[stat_i][stat_j]].open(tempfilename);
				cout<<" file "<<tempfilename<<" opened ."<<endl;
				sprintf(tempfilename, "dis%03d%03d_plot.dat", stat_i+1, stat_j+1);
				cout<<" now writing : "<<tempfilename<<" to stream ~ "<<_CINDEXMAP[stat_i][stat_j]<<" ;-) "<<endl;
				temposDIS_plot[_CINDEXMAP[stat_i][stat_j]].open(tempfilename);
				cout<<" file "<<tempfilename<<" opened ."<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) { //each E bin;
					tempNum_all=0.0;
					tempSum_all=0.0;
					for(stat_l=0; stat_l<_RG_totalnum; stat_l++) {
						tempNum_all+=_DIS_stat_tot.pArray[_CINDEXMAP[stat_i][stat_j]][stat_k][stat_l];
						tempSum_all+=_DIS_stat_tot.pArray[_CINDEXMAP[stat_i][stat_j]][stat_k][stat_l]*(_DISSTAT_lowest+_DISSTAT_interval*stat_l);
					}
					if( fabs(_Probability_all[stat_k])-tempNum_all>1e-12 ) {
						cout<<" something wrong: _Probability_all[stat_k]="<<_Probability_all[stat_k]<<" temp_all="<<tempNum_all<<endl;
						exit(LOGICERROR);
					}
					if(_Probability_all[stat_k]<1e-6) {
						temposDIS[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "<<setw(8)<<0.0<<" "<<setw(8)<<0.0<<endl;
					} else {
						temposDIS[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					             <<setw(10)<<_DIS_stat_actual_tot.pArray[_CINDEXMAP[stat_i][stat_j]][stat_k]/_Probability_all[stat_k]
					             <<setw(10)<<tempSum_all/_Probability_all[stat_k]<<endl;	
					}
					for(stat_l=0; stat_l<_RG_totalnum; stat_l++) {
						if( fabs(tempNum_all)>1e-6 ) {
						   	temposDIS_plot[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_DISSTAT_lowest+double(stat_l)*_DISSTAT_interval<<" "
						                  <<setw(8)<<double(_DIS_stat_tot.pArray[_CINDEXMAP[stat_i][stat_j]][stat_k][stat_l])/double(tempNum_all)<<endl;
						} else {
						    temposDIS_plot[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_DISSTAT_lowest+double(stat_l)*_DISSTAT_interval<<" "
						                  <<setw(8)<<0.0<<endl;
						}
					}
					temposDIS_plot[_CINDEXMAP[stat_i][stat_j]]<<" "<<endl;
				}
				temposDIS[_CINDEXMAP[stat_i][stat_j]].close();
				temposDIS_plot[_CINDEXMAP[stat_i][stat_j]].close();
			}
		}
	}//end of dis output;
	if(_len_EINT>0) {
		if(_FLAG_ener) {
			if(_EINT_totalnum>1) {
				ofstream temposEINT_plot;
				sprintf(tempfilename, "eint%03d%03d_plot.dat", _EINTlist[0]+1, _EINTlist[1]+1);
				cout<<" now writing : "<<tempfilename<<" to stream ~ "<<_EINT_index<<" ;-) "<<endl;
				temposEINT_plot.open(tempfilename);
				cout<<" file "<<tempfilename<<" opened ."<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) { //each E bin;
					tempNum_all=0.0;
					tempSum_all=0.0;
					for(stat_l=0; stat_l<_EINT_totalnum; stat_l++) {
						tempNum_all+=_EINT_stat_tot.pArray[stat_k][stat_l];
						tempSum_all+=_EINT_stat_tot.pArray[stat_k][stat_l]*(_EINT_lowest+_EINT_interval*stat_l);
					}
					if( fabs(_Probability_all[stat_k])-tempNum_all>1e-12 ) {
						cout<<" something wrong: _Probability_all[stat_k]="<<_Probability_all[stat_k]<<" temp_all="<<tempNum_all<<endl;
						exit(LOGICERROR);
					}
					for(stat_l=0; stat_l<_EINT_totalnum; stat_l++) {
						if( fabs(tempNum_all)>1e-6 ) {
						   	temposEINT_plot<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_EINT_lowest+double(stat_l)*_EINT_interval<<" "
						                  <<setw(8)<<double(_EINT_stat_tot.pArray[stat_k][stat_l])/double(tempNum_all)<<endl;
						} else {
						    temposEINT_plot<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_EINT_lowest+double(stat_l)*_EINT_interval<<" "
						                  <<setw(8)<<0.0<<endl;
						}
					}
					temposEINT_plot<<" "<<endl;
				}
				temposEINT_plot.close();
			}
			///add your code here/// for contact number;
			ofstream temposCN;
			sprintf(tempfilename, "CN%03d%03d.dat", _EINTlist[0]+1, _EINTlist[1]+1);
			cout<<" now writing : "<<tempfilename<<" to stream ~ "<<_EINT_index<<" ;-) "<<endl;
			temposCN.open(tempfilename);
			cout<<" file "<<tempfilename<<" opened ."<<endl;
			for(stat_k=0; stat_k<_E_totalnum; stat_k++) { //each E bin;
				if(_Probability_all[stat_k]<1e-6) {
					temposCN<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                  <<setw(8)<<0.0<<endl;
				} else {
					temposCN<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                  <<setw(8)<<_CN_stat_tot[stat_k]/_Probability_all[stat_k]<<endl;
				}
			}
			temposCN.close();
		}
	}//end of eint output;
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		cout<<" now writing chn: "<<stat_i+1<<" sz: "<<stat_size<<endl;
		//stat_head=_INDEX_CHN_HEAD[stat_i];
		//stat_tail=_INDEX_CHN_TAIL[stat_i];
		if(_FLAG_ener) {
			/*totoosELJ.open("EnerLJtot.dat");
			for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
				totoosELJ<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" ";
				for(stat_j=0; stat_j<=stat_i; stat_j++) {
					if(_Probability_all[stat_k]<1e-6) {
						totoosELJ<<setw(10)<<0.0<<" ";
					} else {
						totoosELJ<<setw(10)<<_ENERLJ_stat_tot.pArray[stat_i*_NUM_chains+stat_j][stat_k]/_Probability_all[stat_k]<<" ";
					}
				}
				totoosELJ<<endl;
			}*/
			for(stat_j=0; stat_j<=stat_i; stat_j++) {
				sprintf(tempfilename, "EnerLJ%03d%03d.dat", stat_i+1, stat_j+1);
				cout<<" now writing : "<<tempfilename<<" to stream ~ "<<_CINDEXMAP[stat_i][stat_j]<<" ;-) "<<endl;
				temposELJ[_CINDEXMAP[stat_i][stat_j]].open(tempfilename);
				cout<<" file "<<tempfilename<<" opened ."<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
					if(_Probability_all[stat_k]<1e-6) {
						temposELJ[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "<<setw(10)<<0.0<<endl;
					} else {
						temposELJ[_CINDEXMAP[stat_i][stat_j]]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						      <<setw(10)<<_ENERLJ_stat_tot.pArray[_CINDEXMAP[stat_i][stat_j]][stat_k]/_Probability_all[stat_k]<<endl;
					}
				}
				temposELJ[_CINDEXMAP[stat_i][stat_j]].close();
			}
			if(_BF_flag) {
				////add your code here;
				sprintf(tempfilename, "EnerBF%03d.dat", stat_i+1);
				temposEBF[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
					if(_Probability_all[stat_k]<1e-6) {
						temposEBF[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<0.0<<endl;
					} else { 
						temposEBF[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<_ENERBF_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<endl;
					}
				}
				temposEBF[stat_i].close();
			}
			if(_AG_flag) {
				////add your code here;
				sprintf(tempfilename, "EnerAG%03d.dat", stat_i+1);
				temposEAG[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
					if(_Probability_all[stat_k]<1e-6) {
						temposEAG[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<0.0<<endl;
					} else {
						//cout<< _ENERAG_stat_tot.pArray[stat_i][stat_k]<<endl;
						temposEAG[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<_ENERAG_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<endl;
					}
				}
				temposEAG[stat_i].close();
			}
			if(_DH_flag) {
				////add your code here;
				sprintf(tempfilename, "EnerDH%03d.dat", stat_i+1);
				temposEDH[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;
				for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
					if(_Probability_all[stat_k]<1e-6) {
						temposEDH[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<0.0<<endl;
					} else { 
						temposEDH[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
					                 <<setw(10)<<_ENERDH_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<endl;
					}
				}
			}
			temposEDH[stat_i].close();
			//totoosELJ.close();
		}
		if(_FLAG_rg2) {

			sprintf(tempfilename, "COM%03d.dat", stat_i+1);
			temposCOM[stat_i].open(tempfilename);
			cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;
			for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
				if(_Probability_all[stat_k]<1e-6) {
					temposCOM[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
				                 <<setw(10)<<0.0<<" "
				                 <<setw(10)<<0.0<<" "
				                 <<setw(10)<<0.0<<endl;
				} else { 
					temposCOM[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
				                 <<setw(10)<<_COM_x_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<" "
				                 <<setw(10)<<_COM_y_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<" "
				                 <<setw(10)<<_COM_z_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<endl;
				    /*if(stat_k==4063 && stat_i==1) {
				    	cout<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "<<_Probability_all[stat_k]<<" "
				                 <<setw(10)<<_COM_x_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<" "
				                 <<setw(10)<<_COM_y_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<" "
				                 <<setw(10)<<_COM_z_stat_tot.pArray[stat_i][stat_k]/_Probability_all[stat_k]<<endl;
				    }*/
				}
			}
			temposCOM[stat_i].close();

			//if( !_INDEX_chn_or_not_chain[stat_i] || (_SIZE_of_chn[stat_i] > 2) ) {//necessary? {
			if( _SIZE_of_chn[stat_i] > 1 ) {//necessary? {
				sprintf(tempfilename, "RG2%03d.dat", stat_i+1);
				temposRG2[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;

				/*sprintf(tempfilename, "rg2%03d_xplot.dat", stat_i+1);
				temposRG2_x_plot[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;

				sprintf(tempfilename, "rg2%03d_yplot.dat", stat_i+1);
				temposRG2_y_plot[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;

				sprintf(tempfilename, "rg2%03d_zplot.dat", stat_i+1);
				temposRG2_z_plot[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;*/

				sprintf(tempfilename, "rg2%03d_plot.dat", stat_i+1);
				temposRG2_plot[stat_i].open(tempfilename);
				cout<<" now writing : "<<tempfilename<<" to stream~"<<stat_i<<" ;-) "<<endl;

				for(stat_k=0; stat_k<_E_totalnum; stat_k++) {
					//tempNum_x=0.0;
					//tempNum_y=0.0;
					//tempNum_z=0.0;
					tempNum_all=0.0;
					//tempSum_x=0.0;
					//tempSum_y=0.0;
					//tempSum_z=0.0;
					tempSum_all=0.0;
					for(stat_j=0; stat_j<_RG_totalnum; stat_j++) {
						//tempNum_x+=_RG2_x_stat_tot.pArray[stat_i][stat_k][stat_j];
						//tempNum_y+=_RG2_y_stat_tot.pArray[stat_i][stat_k][stat_j];
						//tempNum_z+=_RG2_z_stat_tot.pArray[stat_i][stat_k][stat_j];
						tempNum_all+=_RG2_stat_tot.pArray[stat_i][stat_k][stat_j];
						//tempSum_x+=_RG2_x_stat_tot.pArray[stat_i][stat_k][stat_j]*(_RG_lowest+double(stat_j)*_RG_interval);
						//tempSum_y+=_RG2_y_stat_tot.pArray[stat_i][stat_k][stat_j]*(_RG_lowest+double(stat_j)*_RG_interval);
						//tempSum_z+=_RG2_z_stat_tot.pArray[stat_i][stat_k][stat_j]*(_RG_lowest+double(stat_j)*_RG_interval);
						tempSum_all+=_RG2_stat_tot.pArray[stat_i][stat_k][stat_j]*(_RG_lowest+double(stat_j)*_RG_interval);
					}
					
					/*if( fabs(_Probability_all[stat_k])-tempNum_x>1e-12 ||
						tempNum_x!=tempNum_y || tempNum_x!=tempNum_z || tempNum_y!=tempNum_z ||
						tempNum_x!=tempNum_all ) {*/
					if( _RG_totalnum>1 && fabs(_Probability_all[stat_k])-tempNum_all>1e-12 ) {
						/*cout<<" something wrong: trgx="<<tempNum_x<<" trgy="<<tempNum_y<<" trgz="<<tempNum_z
							<<" _Probability_all[stat_k]="<<" temp_all="<<tempNum_all<<endl;*/
						cout<<" something wrong: _Probability_all="<<_Probability_all[stat_k]<<" temp_all="<<tempNum_all<<endl;
						exit(LOGICERROR);
					}
					if(_Probability_all[stat_k]<1e-6) {
						temposRG2[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
									<<setw(8)<<0.0<<" "
									<<setw(8)<<0.0<<" "
									<<setw(8)<<0.0<<" "
									<<setw(8)<<0.0<<" "
									<<setw(8)<<0.0<<endl;
						              /*<<setw(8)<<0.0<<" "
						              <<setw(8)<<0.0<<" "
						              <<setw(8)<<0.0<<" "
						              <<setw(8)<<0.0<<endl;*/
					} else {
						temposRG2[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						<<setw(10)<<(_RG2_actual_xtot.pArray[stat_i][stat_k]
							        +_RG2_actual_ytot.pArray[stat_i][stat_k]
							        +_RG2_actual_ztot.pArray[stat_i][stat_k])/stat_size/_Probability_all[stat_k]<<" "
						<<setw(10)<<tempSum_all/_Probability_all[stat_k]<<" "
						<<setw(10)<<_RG2_actual_xtot.pArray[stat_i][stat_k]/stat_size/_Probability_all[stat_k]<<" "
						<<setw(10)<<_RG2_actual_ytot.pArray[stat_i][stat_k]/stat_size/_Probability_all[stat_k]<<" "
						<<setw(10)<<_RG2_actual_ztot.pArray[stat_i][stat_k]/stat_size/_Probability_all[stat_k]<<endl;
					             /*<<setw(10)<<tempSum_x/_Probability_all[stat_k]<<" "
					             <<setw(10)<<tempSum_y/_Probability_all[stat_k]<<" "
					             <<setw(10)<<tempSum_z/_Probability_all[stat_k]<<" "
					             <<setw(10)<<tempSum_all/_Probability_all[stat_k]<<endl;*/	
					}
					for(stat_j=0; stat_j<_RG_totalnum; stat_j++) {
						/*if( fabs(tempNum_x)>1e-6 ) {
							temposRG2_x_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<double(_RG2_x_stat_tot.pArray[stat_i][stat_k][stat_j])/double(tempNum_x)<<endl;
						} else {
							temposRG2_x_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<0.0<<endl;
					   	}
					   	if( fabs(tempNum_y)>1e-6 ) {
						   	temposRG2_y_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<double(_RG2_y_stat_tot.pArray[stat_i][stat_k][stat_j])/double(tempNum_y)<<endl;
						} else {
						    temposRG2_y_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<0.0<<endl;
						}
						if( fabs(tempNum_z)>1e-6 ) {
						   	temposRG2_z_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<double(_RG2_z_stat_tot.pArray[stat_i][stat_k][stat_j])/double(tempNum_z)<<endl;
						} else {
						    temposRG2_z_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<0.0<<endl;
						}*/
						if( fabs(tempNum_all)>1e-6 ) {
						   	temposRG2_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<double(_RG2_stat_tot.pArray[stat_i][stat_k][stat_j])/double(tempNum_all)<<endl;
						} else {
						    temposRG2_plot[stat_i]<<setw(8)<<(_E_lowest+_E_interval*stat_k)<<" "
						                  <<setw(8)<<_RG_lowest+double(stat_j)*_RG_interval<<" "
						                  <<setw(8)<<0.0<<endl;
						}
					}
					//temposRG2_x_plot[stat_i]<<" "<<endl;
					//temposRG2_y_plot[stat_i]<<" "<<endl;
					//temposRG2_z_plot[stat_i]<<" "<<endl;
					temposRG2_plot[stat_i]<<" "<<endl;
				}
				temposRG2[stat_i].close();
				//temposRG2_x_plot[stat_i].close();
				//temposRG2_y_plot[stat_i].close();
				//temposRG2_z_plot[stat_i].close();
				temposRG2_plot[stat_i].close(); 
			}
		}
	}
	delete[] temposCOM;
	delete[] temposELJ;
	delete[] temposEAG;
	delete[] temposEBF;
	delete[] temposEDH;
	delete[] temposRG2; 
	//delete[] temposRG2_x_plot; 
	//delete[] temposRG2_y_plot; 
	//delete[] temposRG2_z_plot; 
	delete[] temposRG2_plot; 
	//string tempstr;
	char aptempa[100];
	char aptempb[100];
	sprintf(aptempa, "echo \"xmin=%f\" > range.gpl", _E_lowest);
	sprintf(aptempb, "echo \"xmax=%f\" >> range.gpl", _E_highest);
	//tempstr="echo \"set xrange ["+string(aptempa)+string(":")+string(aptempb)+"]\" > range.gpl";
	if(system(aptempa)) {};
	if(system(aptempb)) {};
	
	if(system("cp range.gpl rg2maprange.gpl")){};
	sprintf(aptempa, "echo \"ymin=%f\" >> rg2maprange.gpl", _RG_lowest);
	sprintf(aptempb, "echo \"ymax=%f\" >> rg2maprange.gpl", _RG_highest);
	//tempstr="echo \"set yrange ["+string(aptempa)+string(":")+string(aptempb)+"]\" >> rg2maprange.gpl";
	if(system(aptempa)){};
	if(system(aptempb)){};

	sprintf(aptempa, "echo \"dmin=%f\" >> dismaprange.gpl", _DISSTAT_lowest);
	sprintf(aptempb, "echo \"dmax=%f\" >> dismaprange.gpl", _DISSTAT_highest);
	//tempstr="echo \"set yrange ["+string(aptempa)+string(":")+string(aptempb)+"]\" >> rg2maprange.gpl";
	if(system(aptempa)){};
	if(system(aptempb)){};

	sprintf(aptempa, "echo \"cmin=%f\" >> rg2maprange.gpl", _EINT_lowest);
	sprintf(aptempb, "echo \"cmax=%f\" >> rg2maprange.gpl", _EINT_highest);
	//tempstr="echo \"set yrange ["+string(aptempa)+string(":")+string(aptempb)+"]\" >> rg2maprange.gpl";
	if(system(aptempa)){};
	if(system(aptempb)){};
}
////////////////
void cmc::close_statistic() { 
	if(_PROC_ID!=0) {
		return;
	}
	/*char xtemp[30];
	sprintf(xtemp,"%d",_RUNTIMES_totalnum);
	string tempstr;
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		if(_FLAG_ener) {
			//sprintf(tcharn,"rh2%03d%03d%05d",(stat_i)+1,_NUM_chains,_NUM_atoms);
			ifstream tempstream( (string(tcharn)+string(".dat")).c_str() );
			if(getline(tempstream, tempstr)) {
				//system( (string("echo \"fn=\'") + tcharn + string("\';xlow=1;xhigh=")+string(xtemp)+string("\" > rh2.par")).c_str() );
				//system("gnuplot < rh2.gpl");
			}
			cout<<" file: ["<<tcharn<<".dat] closed, please check ["<<tcharn<<".eps]."<<endl;
			//_rh2_stream[stat_i].close();
		}
		if(_FLAG_rg2) {
			//sprintf(tcharn,"rg2%03d%03d%05d",(stat_i)+1,_NUM_chains,_NUM_atoms);
			ifstream tempstream( (string(tcharn)+string(".dat")).c_str() );
			if(getline(tempstream, tempstr)) {
				system( (string("echo \"fn=\'") + tcharn + string("\';xlow=1;xhigh=")+string(xtemp)+string("\" > rg2.par")).c_str() );
				system("gnuplot < rg2.gpl");
			}
			//cout<<" file: ["<<tcharn<<".dat] closed, please check ["<<tcharn<<".eps]."<<endl;
			//_rg2_stream[stat_i].close();
		}
	}   
	delete[] tcharn;
	if(_FLAG_ener) {
		//delete[] _rh2_stream;
	}
	if(_FLAG_rg2) {
		//delete[] _rg2_stream;
	}*/
	_Statistic_over=true;
}
////////////////
void cmc::translate2box() {
	int i=0;
	/*double com_x=0.0;
	double com_y=0.0;
	double com_z=0.0;*/
	for(i=1; i<_SIZE_memo; i++){
		_XX_rec[i]=Coor_PBC_X(_XX_rec[i]);
		_YY_rec[i]=Coor_PBC_Y(_YY_rec[i]);
		_ZZ_rec[i]=Coor_PBC_Z(_ZZ_rec[i]);
		/*com_x+=_XX[i];
		com_y+=_YY[i];
		com_z+=_ZZ[i];*/
	}
	/*com_x/=_NUM_atoms;
	com_y/=_NUM_atoms;
	com_z/=_NUM_atoms;
	for(i=1; i<_SIZE_memo; i++){
		_XX[i]-=com_x;
		_YY[i]-=com_y;
		_ZZ[i]-=com_z;
	}*/
}









////////////////////////////////
