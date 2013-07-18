//cmolecule.h

#ifndef _C_MOLECULE_H
#define _C_MOLECULE_H

//#pragma warning ( disable: 4786 )
#include <fstream>
//#include <string>
//#include <vector>
#include <iomanip>
#include <cstdio>
#include "rwpara.h"
#include "aglcalc.h"
using namespace std;

////////////////////////////////////
void ErrorMSG(string Error_STR);
////////////////////////////////////////
////////////////////////////////////
///////////// atom /////////////////
//////////////////////////////////////////
class CAtom {
public:
	CAtom();
	~CAtom();
	void add_arbitrary_info(string someatom_type,
		                    int someatom_index,
							string someatom_name_chemical,
							string someatom_name_special,
							string someresidue_name,
							string somechain_name,
							int someresidue_iname,
							double xx_coordinate, double yy_coordinate, double zz_coordinate,
							double ZOOMFACTOR);
	void readpdbinfo(const string atominfo);
	void writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz);
	bool operator<(const CAtom &tatom) const {//for sorting atoms.
		return atomindex<tatom.atomindex;
	}
	void memo_free();
	string type;
	int    atomindex;
	string residuename;
	string chainname;
	int    residueiname; //acctually index+1.
	string atomname;
	string atomname_chemical;
	string atomname_special;
	double x_coordinate;
	double y_coordinate;
	double z_coordinate;
};
///////////////////////////////////////////////////////
// residue //
/////////////////////////////////////////////////////
class CResidue {
public:
	CResidue();
	~CResidue();
	void memo_free();
	void add_arbitrary_info(string someatom_type,
		                    int    someatom_index,
							string someatom_name_chemical,
							string someatom_name_special,
							string someresidue_name,
							string somechain_name,
							int someresidue_iname,
							double xx_coordinate, double yy_coordinate, double zz_coordinate,
							double ZOOMFACTOR,
							bool sortingflag);
	void add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag);
	void readpdbinfo(const string atominfo, double ZOOMFACTOR, bool sortingflag);
	bool operator<(const CResidue &tresidue) const {//for sorting residues.
		return residueiname<tresidue.residueiname;
	}
	void writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz);
	//bool operator<(\CAtom &tatomb);  //for sorting atoms.
 	string residuename;
	string chainname;
	int residueiname; //acctually residue index+1;
	int natoms;
	int natoms_subchain; //normally useless ... 
	vector<CAtom> atoms;
	//////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////
class CChain {
public:
	CChain();
	~CChain();
	void memo_free();
	void add_arbitrary_info(string someatom_type,
		                    int    someatom_index,
							string someatom_name_chemical,
							string someatom_name_special,
							string someresidue_name,
							string somechain_name,
							int someresidue_iname,
							double xx_coordinate, double yy_coordinate, double zz_coordinate,
							double ZOOMFACTOR,
							bool sortingflag);
	void add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag);
	void readpdbinfo(const string atominfo, double ZOOMFACTOR, bool sortingflag);
	void writepdbinfo(ofstream &targetstream, ofstream &targetstream_xyz);
	bool operator<(const CChain &tchain) const {//for sorting chains.
		return chainname<tchain.chainname;
	}
	string chainname;
	int index_chn_real;
	int nresidues;
	//int natoms;
	vector<CResidue> residues;	
	//////////////////////////////////////////////////////
};
//////////////////////////////////////////
///// molecule ////////
//////////////////////////////////////////
class CMolecule {
public:
	CMolecule();
	~CMolecule();
	void add_arbitrary_info(string someatom_type,
		                    int    someatom_index,
		                    string someatom_name_chemical,
							string someatom_name_special,
							string someresidue_name,
							string somechain_name,
							int someresidue_iname,
							double xx_coordinate, double yy_coordinate, double zz_coordinate,
							double ZOOMFACTOR,
							bool sortingflag);
	void add_arbitrary_info(CAtom someatom, double ZOOMFACTOR, bool sortingflag);
	void writepdbinfo(const char* ftarget, bool ifverbose);
	void writelmpinfo(const char* ftarget);
	//void sort();
	void readpdbinfo_str(const string atominfo, double ZOOMFACTOR, bool sortingflag);
	void readpdbinfo(const char* fmolname, double ZOOMFACTOR, bool sortingflag, bool ifverbose);
	//void mapping_atom_2PDB(const int index_in_array);
	//void mapping_atom_2CHN();
	//void memo_allocation();
	//void memo_evaluation();
	void memo_free();
	//bool MemoAllocationFlag;
	//bool MemoEvaluationFlag;
///////////////////////////////////////////////////////////////////////////////////////
	int nchains;
	int nresidues;
	//int natoms; //necessary here.
	vector<CChain> chains;
	vector<CResidue> residues;
	vector<CAtom> allatoms;
	std::vector<int> indexatmres;
	//vector<CAtom> allatoms;

  	//double* _XX;
	//double* _YY;
	//double* _ZZ;
	////////////////// Array start /////////////////////////
	//int*    _Index_CHN;               // the chain index of chosen atom;
	//int*    _Index_CHN_real;          // the is_chain_or_not index of chosen atom: 'A' for 0 ... 'Z' for 25;
	//int*    _Index_ATM_in_RES;        // the index of chosen atom in residue;
	//int*    _Index_ATM_in_CHN;        // the index of chosen atom in chain;
	//int*    _Index_RES_in_CHN;        // the residue index of chosen atom in chain;
	//int*    _Index_RES_in_MOL;        // the residue index of chosen atom in molecule;
	//int*    _Type_ATOM;               // at the front / in the middle / at the end of the chain/residue?
	////////////////// Array end //////////////////////////
private:	
	///////////////////////////////////////////////////////////////////////////////////////	
	///////////////////////////////////////////////////////////////////////////////////////	
};

#endif

