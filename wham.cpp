

#include "rwpara.h"
#include "aglcalc.h"
#include "energy.h"


#define _CONVERGENCE 1e-12


void smoothArray(const double *in, double *out, long N, int n) {
    long k, n1, n2;
    int i;
    for(k=0; k<N; k++) {
        n1=k-n;
        if(n1<0) n1=0;
        n2=k+n;
        if(n2>N-1) n2=N-1; 
        out[k]=0;
        for(i=n1; i<=n2; i++) {
            out[k]+=in[i]/(n2-n1+1);
        }
    }
}

int main(int argc, char** argv) {
	if(argc<6) {
		cout<<" cli: wham smoothingsteps epsilon bondlen sigma N [usedrepN]"<<endl;
		exit(-1);
	}
	int smoothsteps=atoi(argv[1]);
	double epsilon=atof(argv[2]);
	double bondlen=atof(argv[3]);
	double sigma=atof(argv[4]);
	double NumAtoms=atof(argv[5]);
	int usedrepN=0;
	if(argc==7) {
		usedrepN=atoi(argv[6]);
	}
	int ie=0;
	int rs=0;
	int nrep=0;
	int nene=0;
	
	string tempstr;
	vector<string> tempvec;

	string FILENAME_para;
	FILENAME_para=string("_config.in");
	ifstream para_ifstream(FILENAME_para.c_str());
	if(para_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
		exit(-1);
	}
	cout<<" Loading parameters from file [ "<<FILENAME_para<<" ] ..."<<endl;
	
	while( getline(para_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		if(tempvec.size()!=0) {
			readparameter(tempstr, string("_NUM_replicas"), nrep);
			readparameter(tempstr, string("_E_totalnum"), nene);
		}
	}
	cout<<" there are originally "<<nrep<<" replicas, and "<<nene<<" energy bins."<<endl;
	cout<<" and we are gonna use "<<usedrepN<<" replicas of it."<<endl;
	para_ifstream.close();

	if(argc>=7) {
		if( usedrepN>nrep || usedrepN<1 ) {
			cout<<" Error: usedrepN setup wrong, should be [1,"<<nrep<<"] "<<endl;
			exit(-1);
		}
		//nrep=usedrepN;
	} else {
		usedrepN=nrep;
	}

	///////////
	CMyArray<double> ProbEach;
	ProbEach.Build(usedrepN,nene);
	ProbEach.SetZero();
	///////////
	///////////
	double* EBins;
	EBins=new double[nene];
	MEMOSETZERO(EBins, sizeof(double)*nene);	
	///////////
	double* ProAll;
	ProAll=new double[nene];
	MEMOSETZERO(ProAll, sizeof(double)*nene);	

	FILENAME_para=string("probability_each.dat");
	ifstream probeach(FILENAME_para.c_str());
	if(probeach==NULL) {
		cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
		exit(IOERROR);
	}
	cout<<" Loading data_each from file [ "<<FILENAME_para<<" ] ..."<<endl;
	ie=0;
	int tempsz=0;
	while( getline(probeach, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		tempsz=tempvec.size();
		if( tempsz!=(nrep+1) ) {
			cout<<" error, tempvec.size()="<<tempvec.size()<<" nrep="<<nrep<<endl;
			exit(-1);
		} else {
			//cout<< tempstr<<endl;
			EBins[ie]=atof(tempvec[0].c_str());
			ProAll[ie]=0.0;
			for(rs=0; rs<usedrepN; rs++) {
				ProbEach.pArray[rs][ie]=atof(tempvec[rs+1].c_str());
				ProAll[ie]+=ProbEach.pArray[rs][ie];
				//if(ProAll[ie]>1e-6) cout<<ProAll[ie]<<endl;
			}
			ie++;
		}
	}
	probeach.close();
	cout<<" all data loaded."<<endl;
	
	double TempMax=0.0;
	//double TempMin=0.0;
	int* _INDEX_maxenerdis_eachrep;
	int* _INDEX_minenerdis_eachrep;
	_INDEX_maxenerdis_eachrep=new int[usedrepN];
	_INDEX_minenerdis_eachrep=new int[usedrepN];
	for(rs=0; rs<usedrepN; rs++) {
		TempMax=0.0;
		//TempMin=0.0;
		_INDEX_maxenerdis_eachrep[rs]=0;
		for(ie=0; ie<nene; ie++) {
			if( ProbEach.pArray[rs][ie] > TempMax ) {
				//cout<<"ProbEach["<<rs<<"]["<<ie<<"]="<<ProbEach.pArray[rs][ie]<<endl;
				TempMax=ProbEach.pArray[rs][ie];
				_INDEX_maxenerdis_eachrep[rs]=ie;
			}
		}
		//cout<<" maxdis["<<rs<<"] at "<<_INDEX_maxenerdis_eachrep[rs]<<" w/ pe="<<ProbEach.pArray[rs][_INDEX_maxenerdis_eachrep[rs]]<<endl;
		_INDEX_minenerdis_eachrep[rs]=0;
		for(ie=_INDEX_maxenerdis_eachrep[rs]; ie>=0; ie--) {
			if( ProbEach.pArray[rs][ie]<1e-6 ) {
				_INDEX_minenerdis_eachrep[rs]=ie+1;
				break;
			}
		}
		if(_INDEX_minenerdis_eachrep[rs]>_INDEX_maxenerdis_eachrep[rs]) {
			_INDEX_minenerdis_eachrep[rs]=_INDEX_maxenerdis_eachrep[rs];
		}
		cout<<" mindis["<<rs<<"] at "<<_INDEX_minenerdis_eachrep[rs]<<" "<<ProbEach.pArray[rs][_INDEX_minenerdis_eachrep[rs]]<<endl;
		for(ie=_INDEX_maxenerdis_eachrep[rs]; ie<nene; ie++) {
			if( ProbEach.pArray[rs][ie]<1e-6 ) {
				_INDEX_maxenerdis_eachrep[rs]=ie-1;
				break;
			} else {
				_INDEX_maxenerdis_eachrep[rs]=ie;
			}
		}
		cout<<" maxdis["<<rs<<"] at "<<_INDEX_maxenerdis_eachrep[rs]<<" "<<ProbEach.pArray[rs][_INDEX_maxenerdis_eachrep[rs]]<<endl;
	}
	
	//////////////// index_max_dis begin /////////////////
	int _INDEX_maxenerdis=_INDEX_maxenerdis_eachrep[0];
	for(rs=1; rs<usedrepN; rs++) {
		//cout<<rs<<endl;
		if(_INDEX_maxenerdis_eachrep[rs]>_INDEX_maxenerdis) {
			_INDEX_maxenerdis=_INDEX_maxenerdis_eachrep[rs];
		}
	}
	//////////////// index_max_dis end /////////////////
	//cout<<"again..."<<endl;
	//////////////// index_min_dis begin ///////////////////
	int _INDEX_minenerdis=_INDEX_minenerdis_eachrep[0];
	for(rs=1; rs<usedrepN; rs++) {
		if( _INDEX_minenerdis_eachrep[rs]!=0 ) {
			_INDEX_minenerdis=_INDEX_minenerdis_eachrep[rs];
			break;
		}
	}
	//cout<<" calc: min_ener_dist: "<<_INDEX_minenerdis<<"::"<<ProAll[_INDEX_minenerdis]<<endl;
	for(rs=1; rs<usedrepN; rs++) {
		if( _INDEX_minenerdis_eachrep[rs]<_INDEX_minenerdis && _INDEX_minenerdis_eachrep[rs]!=0 ) {
			_INDEX_minenerdis=_INDEX_minenerdis_eachrep[rs];
		}
	}
	////////////////////////////////////////////////////	
	cout<<" calc: min_ener_dist: "<<_INDEX_minenerdis<<"::"<<ProAll[_INDEX_minenerdis]<<" E="<<EBins[_INDEX_minenerdis]<<":"<<EBins[_INDEX_minenerdis]/epsilon<<endl;
	cout<<" calc: max_ener_dist: "<<_INDEX_maxenerdis<<"::"<<ProAll[_INDEX_maxenerdis]<<" E="<<EBins[_INDEX_maxenerdis]<<":"<<EBins[_INDEX_minenerdis]/epsilon<<endl;
	delete[] _INDEX_maxenerdis_eachrep;
	delete[] _INDEX_minenerdis_eachrep;

	
	int MinIdx=_INDEX_minenerdis;
	int MaxIdx=_INDEX_maxenerdis;
	
	//MaxIdx=nene-1;//for beautiful graphics. shit. I hate this.
	/*cout<<" max @-> P(E["<<MaxIdx<<"]="<<EBins[MaxIdx]<<")="<<ProAll[MaxIdx]<<endl;
	
	cout<<" proall["<<ie-1<<"]="<<ProAll[ie-1]<<endl;
	cout<<" proall["<<ie<<"]="<<ProAll[ie]<<endl;
	cout<<" proall["<<ie+1<<"]="<<ProAll[ie+1]<<endl;
	MinIdx=ie+1;
	cout<<" min @-> P(E["<<MinIdx<<"]="<<EBins[MinIdx]<<")="<<ProAll[MinIdx]<<endl;*/
	//////////////// index_min_dis end ///////////////////*/

	///////start to calc S(E)///////
	double* Fm;
	Fm=new double[usedrepN];
	MEMOSETZERO(Fm, sizeof(double)*usedrepN);	

	double* Temperatures;
	Temperatures=new double[usedrepN];
	MEMOSETZERO(Temperatures, sizeof(double)*usedrepN);		

	double* _SE;
	_SE=new double[nene];
	MEMOSETZERO(_SE, sizeof(double)*nene);	

	double* _NUM_conf;
	_NUM_conf=new double[usedrepN];
	MEMOSETZERO(_NUM_conf, sizeof(double)*usedrepN);
	for(rs=0; rs<usedrepN; rs++) {
		for(ie=0; ie<nene; ie++) {
			_NUM_conf[rs]+=ProbEach.pArray[rs][ie];
		}
	}

	double SE_MAX=0.0;
	double OldRatio=0.0;
	double NewRatio=0.0;
	//int Numberator=0;
	int MiddleIndex_1=MinIdx+(MaxIdx-MinIdx)/3;
	int MiddleIndex_2=MinIdx+2*(MaxIdx-MinIdx)/3;

	for(ie=0; ie<nene; ie++) {
		_SE[ie]=0.0;   // actually log(_SE[]);
	}

	FILENAME_para=string("_temperaturelist.pls");
	ifstream temperatures_ifstream(FILENAME_para.c_str());
	if(temperatures_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
		exit(-1);
	}
	cout<<" Loading temperatures from file [ "<<FILENAME_para<<" ] ..."<<endl;
	
	tempsz=0;
	while( getline(temperatures_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		if(tempvec.size()!=0) {
			Temperatures[tempsz]=atof(tempvec[0].c_str());
		}
		cout<<" T["<<tempsz+1<<"]="<<Temperatures[tempsz]<<endl;
		tempsz++;
		if(tempsz==usedrepN) {
			break;
		}
	}
	cout<<" load "<<usedrepN<<" temperatures for replicas."<<endl;
	temperatures_ifstream.close();
	//def: for speeding iteration, use the following double** ;
	double** T_rep_E_stp;
	T_rep_E_stp=new double*[usedrepN];

	for(rs=0; rs<usedrepN; rs++) {
		T_rep_E_stp[rs]=new double[nene];
		for(ie=0; ie<nene; ie++) {
			T_rep_E_stp[rs][ie]=1/Temperatures[rs]*EBins[ie];
		}
	}
	//def: end;
	cout<<" now the iteration start, waiting ..."<<endl;

	//iteration: start, to calculate the n(E) using the self-consistent iter;
	while(true) {
		for(rs=0; rs<usedrepN; rs++) {
			Fm[rs]=_MIN_DOUBLE;
			for(ie=0; ie<nene; ie++) {
				Fm[rs]=log_cc( Fm[rs], _SE[ie]-T_rep_E_stp[rs][ie] );
			}
			//cout<<"Fm["<<rs<<"]="<<Fm[rs]<<endl;
		}

		for(ie=0; ie<nene; ie++) {
			_SE[ie]=_MIN_DOUBLE;   // actually log(_SE[]);
			for(rs=0; rs<usedrepN; rs++) {
				if( _NUM_conf[rs] > 0 ) {
					_SE[ie]=log_cc( _SE[ie], log(_NUM_conf[rs])-Fm[rs]-T_rep_E_stp[rs][ie] );
				}
			}
			//cout<<T_rep_E_stp[0][ie]<<" : "<<_E_distribution_allrep[ie]<<"/exp("<<_SE[ie]<<")->";
			if( fabs(ProAll[ie])<1e-6 ) {
				_SE[ie]=_MIN_DOUBLE-_SE[ie];
			} else {
				_SE[ie]=log(ProAll[ie])-_SE[ie];
			}
			/*cout<<_SE[ie]<<endl;
			if(fabs(_SE[ie])>1e-6)
			{
			}*/
		}

		SE_MAX=_SE[0];
		//int tempmaxindex=0;
		for(ie=1; ie<nene; ie++) {
			if(SE_MAX<_SE[ie]) {
				SE_MAX=_SE[ie];
				//tempmaxindex=ie;
			}
		}
		//cout<<ie<<" :::: "<<NE_MAX<<endl;

		for(ie=0; ie<nene; ie++) {
			_SE[ie]=_SE[ie]-SE_MAX; //normalize all _SE[];
		}
		
		//printf(" _SE[%4d]=%030.25f\n", MiddleIndex_2, _SE[MiddleIndex_2] );
		//printf(" _SE[%4d]=%030.25f\n", MiddleIndex_1, _SE[MiddleIndex_1] );
		
		NewRatio=_SE[MiddleIndex_2]-_SE[MiddleIndex_1];
		//printf(" NewRatio=%020.15f OldRatio=%020.15f\n", NewRatio, OldRatio );
		if( fabs(NewRatio-OldRatio) < _CONVERGENCE ) {				
			printf(" fabs(NewRatio-OldRatio)=%030.25f\n", fabs(NewRatio-OldRatio) );
			cout<<" converged!"<<endl;
			break;
		} else {
			OldRatio=NewRatio;
		}
	}
	//iteration: end;

	cout<<" fout entropy to [ entropy.dat ] ... ";
	ofstream entropy_stream( "entropy.dat" );
	if(entropy_stream==NULL) {
		cout<<" can not write file: [ entropy.dat ]. "<<endl;
		exit(-1);
	}
	for(ie=0; ie<nene; ie++) {
		entropy_stream<<setiosflags(ios::left)<<setw(10)
			<<EBins[ie]<<" "<<setw(20)<<_SE[ie]<<resetiosflags(ios::left)<<endl;
	}
	entropy_stream.close();		
	cout<<"done!"<<endl;

	ofstream entgpl_stream( "paraentropy.gpl" );
	if(entgpl_stream==NULL) {
		cout<<" can not write file: [ paraentropy.gpl ]. "<<endl;
		exit(-1);
	}
	double tempmin=0.0;
	double tempmax=0.0;
	for(ie=MinIdx; ie<MaxIdx; ie++) {
		if(ie==MinIdx) {
			tempmin=_SE[ie];
			tempmax=_SE[ie];
		} else {
			if(tempmin>_SE[ie]) {
				tempmin=_SE[ie];
			}
			if(tempmax<_SE[ie]) {
				tempmax=_SE[ie];
			}
		}
	}
	entgpl_stream<<" emin="<<EBins[MinIdx]<<endl;
	entgpl_stream<<" emax="<<EBins[MaxIdx]<<endl;
	entgpl_stream<<" smin="<<tempmin<<endl;
	entgpl_stream<<" smax="<<tempmax<<endl;
	entgpl_stream<<" epsilon="<<epsilon<<endl;
	entgpl_stream<<" bondlen="<<bondlen<<endl;
	entgpl_stream<<" sigma="<<sigma<<endl;
	entgpl_stream<<" numatoms="<<NumAtoms<<endl;
	entgpl_stream.close();
	if(system("gnuplot draw_entropy.gpl")) {
		cout<<" now you may check [ entropy.eps ]. "<<endl;
	}

	///////////// beta calc start //////////////
	double* _PARA_beta;
	_PARA_beta=new double[nene];
	MEMOSETZERO(_PARA_beta, sizeof(double)*nene);

	double _E_interval=EBins[1]-EBins[0];	

	for(ie=MaxIdx; ie<nene; ie++) {
		_PARA_beta[ie]=1.0/Temperatures[0];   // the highest temperature is imposed on node0;
		//why this value, because this can ensure _beta to be positive, 
		//if _beta is negtive, bigener*_beta+_alpha will be very small, and accepted easily,
		//then problems coming, like to many bigener sample, and "factor==0.5" things happen!
	}
	for(ie=MaxIdx-1; ie>MinIdx; ie--) {
		_PARA_beta[ie]=(_SE[ie+1]-_SE[ie])/_E_interval;
	}
	for(ie=MinIdx; ie>=0; ie--) {
		_PARA_beta[ie]=_PARA_beta[MinIdx];
	}
	double* _PARA_beta_T;
	_PARA_beta_T=new double[nene];
	MEMOSETZERO(_PARA_beta_T, sizeof(double)*nene);


	cout<<" smoothing beta, smooth steps: "<<smoothsteps;
	smoothArray(_PARA_beta, _PARA_beta_T, nene, smoothsteps);
	cout<<" done."<<endl;
	/*for(ie=0; ie<nene; ie++) {
		_PARA_beta[ie]=_PARA_beta_T[ie];
	}*/

	cout<<" fout beta to [ beta.dat ] ... ";
	ofstream beta_stream( "beta.dat" );
	if(beta_stream==NULL) {
		cout<<" can not write file: [ beta.dat ]. "<<endl;
		exit(-1);
	}
	for(ie=0; ie<nene; ie++) {
		beta_stream<<setiosflags(ios::left)<<setw(10)
			<<EBins[ie]<<" "<<setw(20)<<_PARA_beta[ie]<<resetiosflags(ios::left)<<endl;
	}
	beta_stream.close();		
	cout<<"done!"<<endl;

	cout<<" fout beta to [ beta_smoothed.dat ] ... ";
	ofstream betasm_stream( "beta_smoothed.dat" );
	if(betasm_stream==NULL) {
		cout<<" can not write file: [ beta_smoothed.dat ]. "<<endl;
		exit(-1);
	}
	for(ie=0; ie<nene; ie++) {
		betasm_stream<<setiosflags(ios::left)<<setw(10)
			<<EBins[ie]<<" "<<setw(20)<<_PARA_beta_T[ie]<<resetiosflags(ios::left)<<endl;
	}
	betasm_stream.close();		
	cout<<"done!"<<endl;

	ofstream betgpl_stream( "parabeta.gpl" );
	if(betgpl_stream==NULL) {
		cout<<" can not write file: [ parabeta.gpl ]. "<<endl;
		exit(-1);
	}
	tempmin=0.0;
	tempmax=0.0;
	for(ie=MinIdx; ie<MaxIdx; ie++) {
		if(ie==MinIdx) {
			tempmin=_PARA_beta_T[ie];
			tempmax=_PARA_beta_T[ie];
		} else {
			if(tempmin>_PARA_beta_T[ie]) {
				tempmin=_PARA_beta_T[ie];
			}
			if(tempmax<_PARA_beta_T[ie]) {
				tempmax=_PARA_beta_T[ie];
			}
		}
	}
	beta_stream.close();
	betgpl_stream<<" emin="<<EBins[MinIdx]<<endl;
	betgpl_stream<<" emax="<<EBins[MaxIdx]<<endl;
	betgpl_stream<<" bmax="<<tempmax<<endl;
	//betgpl_stream<<" bmax="<<(_PARA_beta[MiddleIndex_1]-_PARA_beta[nene-1])*1.5+_PARA_beta[nene-1]<<endl;
	betgpl_stream<<" bmin="<<tempmin<<endl;
	betgpl_stream<<" epsilon="<<epsilon<<endl;
	betgpl_stream<<" bondlen="<<bondlen<<endl;
	betgpl_stream<<" sigma="<<sigma<<endl;
	betgpl_stream<<" numatoms="<<NumAtoms<<endl;
	betgpl_stream.close();
	if(system("gnuplot draw_beta.gpl")) {
		cout<<" now you may check [ beta.eps ]. "<<endl;
	}

	//free: temp memory:
	delete[] Fm;
	Fm=NULL;

	for(rs=0; rs<usedrepN; rs++) {
		delete[] T_rep_E_stp[rs];
		T_rep_E_stp[rs]=NULL;
	}
	delete[] T_rep_E_stp;
	T_rep_E_stp=NULL;

	delete[] _SE;
	_SE=NULL;

	delete[] _PARA_beta;
	_PARA_beta=NULL;

	delete[] _PARA_beta_T;
	_PARA_beta_T=NULL;

	delete[] _NUM_conf;
	_NUM_conf=NULL;

	delete[] Temperatures;
	Temperatures=NULL;
	//free: end


	ProbEach.Release();
	delete[] EBins;
	EBins=NULL;
	delete[] ProAll;
	ProAll=NULL;
	//////////////
	return 0;
}
