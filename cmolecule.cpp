//structure.io.cpp

#include "cmolecule.h"
////////////////////////////////////

/////////////////////////////////////////////////////
void ErrorMSG(string Error_STR) {
	cout<<endl<<string(" Error: ")<<Error_STR<<endl;
}
/////////////////////////////////////////////////////
//////// atom
CAtom::CAtom() {
	residueiname=0;
	x_coordinate=0.0;
	y_coordinate=0.0;
	z_coordinate=0.0;
}
/////////////////////////
CAtom::~CAtom() {	
	this->memo_free();
}
void CAtom::memo_free() {	
}
/////////////////////////
void CAtom::add_arbitrary_info(string someatom_type,
							   int    someatom_index,
							   string someatom_name_chemical,
							   string someatom_name_special,
							   string someresidue_name,
							   string somechain_name,
							   int someresidue_iname, //useless and be careful
							   double xx_coordinate, double yy_coordinate, double zz_coordinate,
							   double ZOOMFACTOR) {
	type=someatom_type;
	atomindex=someatom_index;
	atomname_chemical=someatom_name_chemical;
	atomname_special=someatom_name_special;
	atomname=atomname_chemical+atomname_special;
	residuename=someresidue_name;
	chainname=somechain_name;
	residueiname=someresidue_iname;
	x_coordinate=xx_coordinate*ZOOMFACTOR;
	y_coordinate=yy_coordinate*ZOOMFACTOR;
	z_coordinate=zz_coordinate*ZOOMFACTOR;
}
/////////////////////////
void CAtom::readpdbinfo(const string atominfo, const string atominfo_xyz) {
	type=BFilter( atominfo.substr(0, 6) );
	atomindex=atoi( BFilter(atominfo.substr(6, 5)).c_str() );
	atomname_chemical=BFilter(atominfo.substr(12, 2));
	atomname_special=BFilter(atominfo.substr(14, 2));
	atomname=atomname_chemical+atomname_special;
	residuename=BFilter( atominfo.substr(17, 3) );
	chainname=BFilter( atominfo.substr(21, 1) );
	residueiname=atoi( BFilter( atominfo.substr(22, 4) ).c_str() );
	x_coordinate=atof( BFilter( atominfo.substr(30, 8) ).c_str() );
	y_coordinate=atof( BFilter( atominfo.substr(38, 8) ).c_str() );
	z_coordinate=atof( BFilter( atominfo.substr(46, 8) ).c_str() );

	x_coordinate=atof( BFilter( atominfo_xyz.substr(0,  30) ).c_str() );
	y_coordinate=atof( BFilter( atominfo_xyz.substr(31, 30) ).c_str() );
	z_coordinate=atof( BFilter( atominfo_xyz.substr(62, 30) ).c_str() );

}
/////////////////////////
void CAtom::writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz) {
	//string -> char* is necessary, cos' some MPI can not support <string> well! (like 102...)
	targetstream<<setw(6)<<setiosflags(ios::left)<<type.c_str()<<resetiosflags(ios::left)
		<<setw(5)<<atomindex<<" "
		<<setw(2)<<atomname_chemical.c_str()
		<<setw(2)<<setiosflags(ios::left)<<atomname_special.c_str()<<resetiosflags(ios::left)<<" "
		<<setw(3)<<residuename.c_str()<<" "
		<<chainname.c_str()
		<<setw(4)<<residueiname<<setw(4)<<" "
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<x_coordinate-double(int(x_coordinate)/1000)*1000.0
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<y_coordinate-double(int(y_coordinate)/1000)*1000.0
		<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<z_coordinate-double(int(z_coordinate)/1000)*1000.0
		<<setw(26)<<" "<<endl;
	targetstream_xyz<<setw(30)<<setiosflags(ios::fixed)<<setprecision(22)<<x_coordinate<<" "
		            <<setw(30)<<setiosflags(ios::fixed)<<setprecision(22)<<y_coordinate<<" "
					<<setw(30)<<setiosflags(ios::fixed)<<setprecision(22)<<z_coordinate<<endl;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
////////end atom

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////residue
CResidue::CResidue() {		
	residueiname=0;
	natoms=0;
	natoms_subchain=0;
	/////////////////////
}
/////////////////////////
CResidue::~CResidue() {
	this->memo_free();
	//natoms=0;
}
/////////////////////////
void CResidue::add_arbitrary_info(string someatom_type,
								  int    someatom_index,
								  string someatom_name_chemical,
								  string someatom_name_special,
								  string someresidue_name,
								  string somechain_name,
								  int someresidue_iname,
								  double xx_coordinate, double yy_coordinate, double zz_coordinate,
								  double ZOOMFACTOR,
								  bool sortingflag) {
	CAtom tempatom;
	tempatom.add_arbitrary_info(someatom_type,
		someatom_index,
		someatom_name_chemical,
		someatom_name_special,
		someresidue_name,
		somechain_name,
		someresidue_iname,
		xx_coordinate, yy_coordinate, zz_coordinate,
		ZOOMFACTOR);
	atoms.push_back(tempatom);
	natoms+=1;
	if(natoms==1) {
		residueiname=someresidue_iname;
		chainname=somechain_name;
		residuename=someresidue_name;
	}
	if(sortingflag) {
		std::sort(atoms.begin(),atoms.end());
	}
}
/////////////////////////
void CResidue::add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag) {
	add_arbitrary_info(someatom.type,
		               someatom.atomindex,
		               someatom.atomname_chemical,
					   someatom.atomname_special,
					   someatom.residuename,
					   someatom.chainname,
					   someatom.residueiname,
					   someatom.x_coordinate, 
					   someatom.y_coordinate,
					   someatom.z_coordinate,
					   ZOOMFACTOR,
					   sortingflag);
}
/////////////////////////
void CResidue::readpdbinfo(const string atominfo, const string atominfo_xyz, double ZOOMFACTOR, bool sortingflag) {
	CAtom tempatom;
	tempatom.readpdbinfo(atominfo, atominfo_xyz);
	add_arbitrary_info(tempatom, ZOOMFACTOR, sortingflag);
}
/////////////////////////
void CResidue::writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz) {
	int temp_int=atoms.size();
	for(int i=0; i!=temp_int; i++) {
		atoms[i].writepdbinfo(targetstream, targetstream_xyz);
	}
}
/////////////////////////
void CResidue::memo_free() {
	if(natoms!=0){
		for(int i=0; i<natoms; i++) {
			atoms[i].memo_free();
		}
		atoms.clear(); 
	}
	natoms=0;
}
///////////////////////////////////////////////////////
/////////////////end residue

//////////////start chain/////////////////////////
CChain::CChain() {
	//natoms=0;
	nresidues=0;
}
/////////////////////////
CChain::~CChain() {
	this->memo_free();
	//nresidues=0;
	//natoms=0;
}
/////////////////////////
void CChain::add_arbitrary_info(string someatom_type,
								int    someatom_index,
								string someatom_name_chemical,
								string someatom_name_special,
								string someresidue_name,
								string somechain_name,
								int someresidue_iname,
								double xx_coordinate, double yy_coordinate, double zz_coordinate,
								double ZOOMFACTOR,
								bool sortingflag) {
	int index_residue=0;
	bool flag=false;
	for(index_residue=0; index_residue!=nresidues; index_residue++) {
		if( ( residues[index_residue].residuename==someresidue_name ) 
			&& ( residues[index_residue].residueiname==someresidue_iname )
			&& ( residues[index_residue].chainname==somechain_name ) ) {
			flag=true;
			break;
		}
	}
	if(flag==true) {
		residues[index_residue].add_arbitrary_info(someatom_type, 
			                                       someatom_index,
			                                       someatom_name_chemical, 
												   someatom_name_special,
												   someresidue_name,
												   somechain_name,
												   someresidue_iname,
												   xx_coordinate, yy_coordinate, zz_coordinate,
												   ZOOMFACTOR,
												   sortingflag);
		//natoms+=1;
	} else {
		CResidue TempResidue;
		TempResidue.add_arbitrary_info(someatom_type,
			                           someatom_index,
			                           someatom_name_chemical,
									   someatom_name_special, 
									   someresidue_name,
									   somechain_name,
									   someresidue_iname, 
									   xx_coordinate, yy_coordinate, zz_coordinate,
									   ZOOMFACTOR,
									   sortingflag);
		residues.push_back(TempResidue);
		nresidues+=1;
		//natoms+=1;
	}
	if( nresidues==1 && residues[0].natoms==1 ) {
		chainname=somechain_name;
		index_chn_real=chainname.c_str()[0]-'A';
	}
	if(sortingflag) {
		std::sort(residues.begin(),residues.end());
	}
}
/////////////////////////
void CChain::add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag) {
	add_arbitrary_info(someatom.type,
		               someatom.atomindex,
		               someatom.atomname_chemical,
					   someatom.atomname_special,
					   someatom.residuename,
					   someatom.chainname,
					   someatom.residueiname,
					   someatom.x_coordinate*ZOOMFACTOR, 
					   someatom.y_coordinate*ZOOMFACTOR,
					   someatom.z_coordinate*ZOOMFACTOR,
					   ZOOMFACTOR,
					   sortingflag);
}
/////////////////////////
void CChain::readpdbinfo(const string atominfo, const string atominfo_xyz, double ZOOMFACTOR, bool sortingflag) {
	CAtom tempatom;
	tempatom.readpdbinfo(atominfo, atominfo_xyz);
	add_arbitrary_info(tempatom, ZOOMFACTOR, sortingflag);
}
/////////////////////////
void CChain::writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz) {
	int i=0; 
	//int j=0;
	//int totoalnum=0;
	int temp_int=residues.size();
	if(temp_int!=nresidues) {
		cout<<" there is some error: residues.size()="<<residues.size()<<" != nresidues="<<nresidues<<endl;
		cout<<" @ CChain::writepdbinfo."<<endl;
		exit(SIZEERROR);
	}
	for(i=0; i!=nresidues; i++) {
		residues[i].writepdbinfo(targetstream, targetstream_xyz);
		//totoalnum+=residues[i].natoms;
	}
	/*if(totoalnum!=natoms) {
		cout<<" there is some error: totoalnum="<<totoalnum<<" != natoms="<<natoms<<endl;
		cout<<" @ CChain::writepdbinfo."<<endl;
		exit(SIZEERROR);
	}*/
}
void CChain::memo_free() {
	if(nresidues!=0) {
		for(int i=0; i<nresidues; i++) {
			residues[i].memo_free();
		}
		residues.clear(); 
	}
	nresidues=0;
}
//////////////////////////////////////////////////////
///////////////end chain//////////////////////////

////////////////////////////////////
//////////// molecule //////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CMolecule::CMolecule() {
	//natoms=0;
	nchains=0;
	nresidues=0;
	//indexatmres.clear();
	/*MemoAllocationFlag=false;
	MemoEvaluationFlag=false;
	_XX=NULL;
	_YY=NULL;
	_ZZ=NULL;
	_Index_CHN=NULL;
	_Index_CHN_real=NULL;
	_Index_ATM_in_RES=NULL;
	_Index_ATM_in_CHN=NULL;
	_Index_RES_in_CHN=NULL;
	_Index_RES_in_MOL=NULL;
	_Type_ATOM=NULL;*/
	//atomseries=NULL;
}
/////////////////////////
CMolecule::~CMolecule() {
	this->memo_free();
}
/////////////////////////
void CMolecule::add_arbitrary_info(string someatom_type,
								   int    someatom_index,
								   string someatom_name_chemical,
								   string someatom_name_special,
								   string someresidue_name,
								   string somechain_name,
								   int someresidue_iname,
								   double xx_coordinate, double yy_coordinate, double zz_coordinate,
								   double ZOOMFACTOR,
								   bool sortingflag) {
	int index_chain=0;
	bool flag=false;
	CAtom someatom;
	someatom.add_arbitrary_info(someatom_type, 
			                    someatom_index,
			                    someatom_name_chemical, 
							    someatom_name_special,
							    someresidue_name,
							    somechain_name,
							    someresidue_iname,
							    xx_coordinate, yy_coordinate, zz_coordinate,
							    ZOOMFACTOR);
	for(index_chain=0; index_chain<nchains; index_chain++) {
		if( chains[index_chain].chainname==somechain_name ) {
			flag=true;
			break;
		}
	}
	if(flag==true) {
		chains[index_chain].add_arbitrary_info(someatom, ZOOMFACTOR, sortingflag);
		//natoms+=1;
	} else {
		CChain TempChain;
		TempChain.add_arbitrary_info(someatom, ZOOMFACTOR, sortingflag);
		chains.push_back(TempChain);
		nchains+=1;
		//natoms+=1;
	}

	if(sortingflag) {
		std::sort(chains.begin(),chains.end());
	}
	allatoms.push_back(someatom);

	//cout<<" natoms="<<allatoms.size()<<endl;
	//cout<<" nchains="<<chains.size()<<endl;
	flag=false;
	int index_res=0;
	for(index_res=0; index_res<nresidues; index_res++) {
		if( residues[index_res].residuename==someresidue_name && residues[index_res].residueiname==someresidue_iname) {
			flag=true;
			break;
		}
	}
	if(flag==true) {
		residues[index_res].add_arbitrary_info(someatom, ZOOMFACTOR, sortingflag);
		//cout<<" new atom added to the old res"<<endl;
	} else {
		CResidue TempResidue;
		TempResidue.add_arbitrary_info(someatom, ZOOMFACTOR, sortingflag);
		residues.push_back(TempResidue);
		nresidues+=1;
		//cout<<" new atom added to the new res"<<endl;
	}
	//cout<<" nresidues="<<residues.size()<<endl;
	//cout<<" idxnum="<<index_res<<endl;
	//cout<<" pushing "<<index_res<<" into the reslist..."<<endl;
	indexatmres.push_back(index_res);
	//cout<<" current resindex="<<indexatmres[indexatmres.size()-1]<<endl;
}
/////////////////////////
void CMolecule::add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag) {
	add_arbitrary_info(someatom.type,
		               someatom.atomindex,
		               someatom.atomname_chemical,
					   someatom.atomname_special,
					   someatom.residuename,
					   someatom.chainname,
					   someatom.residueiname,
					   someatom.x_coordinate, 
					   someatom.y_coordinate,
					   someatom.z_coordinate,
					   ZOOMFACTOR,
					   sortingflag);
}
/////////////////////////
void CMolecule::readpdbinfo_str(const string atominfo, const string atominfo_xyz, double ZOOMFACTOR, bool sortingflag) {
	CAtom tempatom;
	tempatom.readpdbinfo(atominfo, atominfo_xyz);
	add_arbitrary_info(tempatom, ZOOMFACTOR, sortingflag);
	//allatoms.push_back(tempatom);
}
/////////////////////////
void CMolecule::writepdbinfo(const char* ftarget, bool ifverbose) {
	if(ifverbose) {
		cout<<endl<<" Warning: consider the distance between each atom pair, "<<endl
		          <<"          there will be memory-problem if too short, check it!"<<endl;	
	}
	ofstream outputfile(ftarget);
	if(outputfile==NULL) {
		cout<<" can not open file: "<<ftarget<<endl;
		exit(IOERROR);
	}
	string TTStr=string(ftarget);
	int TTSize=TTStr.size();
	ofstream outputfile_xyz( (TTStr.substr(0, TTSize-3)+string("xyz")).c_str() );
	if(outputfile_xyz==NULL) {
		cout<<" can not open file: [ "<<TTStr<<" ];"<<endl;
		exit(IOERROR);
	}
	int i=0;
	int j=0;
	int k=0;
	int Num_RES=0;
	int Num_ATM=0;
	int tempnumer=0;
	if(ifverbose==true) {
		cout<<" ( show XYZ info )"<<endl;
		cout<<"        "<<setw(6)<<"index"
			<<setw(9)<<" X_coor"<<setw(9)<<" Y_coor"<<setw(9)<<" Z_coor"<<endl;
		cout<<" nchain: "<<nchains<<" sizeofatoms:"<<allatoms.size()<<endl;
		for(i=0; i!=nchains; i++) {
			cout<<"i="<<i;
			for(j=0; j!=chains[i].nresidues; j++) {
				cout<<"j="<<j;
				for(k=0; k!=chains[i].residues[j].natoms; k++) {
					cout<<"k="<<k;
					/*cout<<"  index:"<<setw(6)<<i+1
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<chains[i].residues[j].atoms[k].x_coordinate
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<chains[i].residues[j].atoms[k].y_coordinate
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<chains[i].residues[j].atoms[k].z_coordinate
						<<endl;*/
					cout<<"  index:"<<setw(6)<<i+1
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<allatoms[tempnumer].x_coordinate
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<allatoms[tempnumer].y_coordinate
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
						<<allatoms[tempnumer].z_coordinate
						<<endl;
					tempnumer++;
				}
			}
		}
	}
	/*for(i=0; i!=nchains; i++) {
		std::ios::sync_with_stdio(false);
		chains[i].writepdbinfo(outputfile, outputfile_xyz);
	}*/
	int tempnumtot=allatoms.size();
	for(i=0; i!=tempnumtot; i++) {
		std::ios::sync_with_stdio(false);
		allatoms[i].writepdbinfo(outputfile, outputfile_xyz);
	}
	tempnumer=0;
	for(i=0; i!=nchains; i++) {
		Num_RES=chains[i].nresidues;
		for(j=0; j<Num_RES; j++) {
			Num_ATM=chains[i].residues[j].natoms;
			for(k=0; k<Num_ATM; k++) {
				if(tempnumer==tempnumtot-1) {
					break;
				}
				//cout<<" writing the "<<tempnumer+1<<"th bonds info..."<<endl;
				if( allatoms[tempnumer].chainname==allatoms[tempnumer+1].chainname ) {
					/*outputfile<<"CONECT"<<setw(5)<<chains[i].residues[j].atoms[k].atomindex
				    	                <<setw(5)<<chains[i].residues[j+1].atoms[0].atomindex<<endl;*/
					outputfile<<"CONECT"<<setw(5)<<allatoms[tempnumer].atomindex;
					tempnumer++;
					outputfile<<setw(5)<<allatoms[tempnumer].atomindex<<endl;
				} else {
					tempnumer++;
				}
			}
		}
	}
	outputfile.close();
	outputfile_xyz.close();
	cout<<" Molecule write into [ "<<ftarget<<" ]! Check the infomation u need!"<<endl;
	cout<<" Molecule write into [ "<<TTStr.substr(0, TTSize-3)+string("xyz")
		<<" ]! Check the infomation u need!"<<endl;
	if(system( (string("python pdb2psf.py ")+string(ftarget)).c_str() )) {};
	writelmpinfo((string(ftarget)+string(".dat")).c_str());
}
void CMolecule::writelmpinfo(const char* ftarget) {
	FILE* fq;
	fq = fopen(ftarget,"w");
	if(fq==NULL) {
		cout<<" can not open file: "<<ftarget<<endl;
		exit(IOERROR);
	}
	printf("entered\n");
	fprintf(fq,"# Model for PE\n\n");
	int allsize=allatoms.size();
	fprintf(fq,"%10d     atoms\n",allsize);
	fprintf(fq,"%10d     bonds\n",allsize-nchains);
	
	fprintf(fq, "\n");
	fprintf(fq,"%10d     atom types\n",nchains);
	fprintf(fq,"%10d     bond types\n",1);
	
	fprintf(fq, "\n");
	float xlo=0.0;
	float xhi=0.0;
	float ylo=0.0;
	float yhi=0.0;
	float zlo=0.0;
	float zhi=0.0;
	fprintf(fq,"%10.4f%10.4f xlo xhi\n",xlo,xhi);
	fprintf(fq,"%10.4f%10.4f ylo yhi\n",ylo,yhi);
	fprintf(fq,"%10.4f%10.4f zlo zhi\n\n",zlo,zhi);
	fprintf(fq,"Masses\n\n");
	float mass=56.0;
	int i=0;
	for(i=0; i<nchains; i++) {
		fprintf(fq,"%10d %14.2f\n",i+1,mass); 
	}
	fprintf(fq,"\nAtoms\n\n");
	int j=0;
	int k=0;
	int tempnumer=0;
	cout<<" nchain: "<<nchains<<" sizeofatoms:"<<allatoms.size()<<endl;
	for(i=0; i!=nchains; i++) {
		//cout<<"i="<<i;
		for(j=0; j!=chains[i].nresidues; j++) {
			//cout<<"j="<<j;
			for(k=0; k!=chains[i].residues[j].natoms; k++) {
				//cout<<"k="<<k;
				fprintf(fq,"%10d%10d%10d%10.4f%10.4f%10.4f\n",tempnumer+1,j+1,i+1,
					 allatoms[tempnumer].x_coordinate,allatoms[tempnumer].y_coordinate,allatoms[tempnumer].z_coordinate);
				tempnumer++;
			}
		}
	}
	cout<<"writing bonds"<<endl;
	fprintf(fq,"\nBonds \n\n");
	int count=1;
	int atom_1=1;
	int atom_2=2;
	//int k=0; //chain index;
	for(i=0;i<allsize-1;i++) { 
  		fprintf(fq,"%10d%10d%10d%10d\n",count,1,atom_1,atom_2);
  		count++;
  		atom_1++;
 		atom_2++;
 		//cout<<count-1<<" "<<count<<endl;
  		if(allatoms[count-2].chainname!=allatoms[count-1].chainname) {
  			atom_1++;
 			atom_2++;
 			i++;
 		}
  		
	}
	cout<<" Molecule write into [ "<<ftarget<<" ] for lammps! Check the infomation u need!"<<endl;
	fclose(fq);
}
/////////////////////////
void CMolecule::readpdbinfo(const char* fmolname, double ZOOMFACTOR, bool sortingflag, bool ifverbose) {
	memo_free();
	cout<<" Reading mol info from [ "<<fmolname<<" ]"<<endl;
	ifstream fmolecule(fmolname);
	if(fmolecule==NULL) {
		cout<<" Error: can not open file: "<<fmolname<<endl;
		exit(IOERROR);
	}
	string TTStr=string(fmolname);
	int TTSize=TTStr.size();
	bool xyzinfoflag=false;
	ifstream fmolecule_xyz( (TTStr.substr(0, TTSize-3)+string("xyz")).c_str() );
	if(fmolecule_xyz==NULL) {
		fmolecule_xyz.close();
		cout<<" with no XYZ info"<<endl;
	} else {
		xyzinfoflag=true;
		cout<<" with XYZ info: [ "<<TTStr.substr(0, TTSize-3)+string("xyz")<<" ]"<<endl;
	}
	int i=0;
	string tempstr;
	string tempstr_xyz;
	//CAtom tempatom;
	vector<string> tempvec;
	vector<string> templist;
	while(getline(fmolecule, tempstr)) {
		
		if( tempstr.substr(0, 4)==string("ATOM") || tempstr.substr(0, 6)==string("HETATM") ) {
			//cout<<tempstr<<endl;
			if(getline(fmolecule_xyz, tempstr_xyz)) {
			} else {
				cout<<" < * > error, no enought xyz infor! exit!"<<endl;
				exit(-1);
			}
			readpdbinfo_str(tempstr, tempstr_xyz, ZOOMFACTOR, sortingflag);
		}
		/*if( tempstr.substr(0, 6)==string("CONECT") ) {
			templist.push_back(tempstr);
		}*/
	}
	/*int templen=templist.size();
	int j=0;
	int tempressz;
	if(templen!=0) {
		bonds=new int*[templi];
	} else {
		templen=0;
		for(i=0;i<nchains;i++){
			tempressz=chains[i].nresidues;
			for(j=0;j<tempressz;j++) {
				templen+=chains[i].residues[j].natoms;
			}
			templen-=2;
		}
	}*/

	int temp_int=chains.size();
	if( nchains!=temp_int ) {
		cout<<" Error: there is some error: nchains="<<nchains<<" != chains.size="<<chains.size()<<endl;
		cout<<" @ CMolecule::readpdbinfo."<<endl;
		exit(SIZEERROR);
	}
	int atomnumber=0;	//for check!
	int j=0;
	for(i=0; i!=nchains; i++) {
		for(j=0; j!=chains[i].nresidues; j++) {
			atomnumber+=chains[i].residues[j].natoms;
		}
	}
	/*if( natoms!=atomnumber ) {		
		cout<<" there is some error: natoms="<<natoms<<" != atomnumber="<<atomnumber<<endl;
		cout<<" @ CMolecule::readpdbinfo."<<endl;
		exit(-2);
	}*/
	int nresidues=0;
	for(i=0; i!=nchains; i++) {
		nresidues+=chains[i].nresidues;
	}
	/*memo_allocation();
	memo_evaluation();*/
	cout<<" This molecule has "<<nchains<<" chains, "
		<<nresidues<<" residues, and "<<atomnumber<<" atoms:"<<endl;
	if(sortingflag) {
		std::sort(allatoms.begin(), allatoms.end());
	}
	//atomseries=new CAtom*[atomnumber];
	int Num_ATM=0;
	int Num_RES=0;
	//int k=0;
	//int numert=0;
	for(i=0; i<nchains; i++) {
		Num_RES=chains[i].nresidues;
		for(j=0; j<Num_RES; j++) {
			Num_ATM=chains[i].residues[j].natoms;
			printf(" [:]%1s[:]%04d[:]%3s[:]%05d\n", 
				   chains[i].chainname.c_str(),
				   j+1,
				   chains[i].residues[j].residuename.c_str(), 
				   Num_ATM);
			//for(k=0; k<Num_ATM; k++) {
				//atomseries[numert++]=&(chains[i].residues[j].atoms[k]);
			//}
		}
	}
	/*cout<<" ----------- [ "<<fmolname<<" ] loaded -----------"<<endl;
	cout<<" Building XYZ array data ..."<<endl;
	if(ifverbose==true) {
		cout<<" ( show XYZ info )"<<endl;
		cout<<"        "<<setw(6)<<"index"<<setw(9)<<" X_coor"<<setw(9)<<" Y_coor"<<setw(9)<<" Z_coor"<<endl;
		for(i=0, atomnumber=1; i!=nchains; i++) {
			Num_RES=chains[i].nresidues;
			for(j=0; j!=Num_RES; j++) {
				Num_ATM=chains[i].residues[j].natoms;
				for(k=0; k!=Num_ATM; k++) {
					cout<<"  index:"<<setw(6)<<atomnumber
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)<<_XX[atomnumber]
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)<<_YY[atomnumber]
						<<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)<<_ZZ[atomnumber]
						<<endl;
					atomnumber++;
				}
			}
		}
	}
	if(xyzinfoflag) {
		fmolecule_xyz.close();
	}*/
	//cout<<" XYZ array build!"<<endl;
	fmolecule.close();
}
///////////////////////////
void CMolecule::memo_free() {
	if(nchains!=0) {
		for(int i=0; i<nchains; i++) {
			chains[i].memo_free();
		}
		chains.clear(); 
		allatoms.clear();
	}
	nchains=0;
}
////// end molecule /////
/////////////////////////

