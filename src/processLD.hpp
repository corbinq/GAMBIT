#ifndef PROCESSLD_HPP
#define PROCESSLD_HPP

#include "tabixpp/tabix.hpp"
#include "tsl/hopscotch_map.h"
#include "tsl/hopscotch_set.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Eigenvalues"

using namespace std;

//extern Tabix tfile;
//extern string chr_global;
//extern string vcf_path;

//extern unordered_map<string, int> iCHR;
//extern vector<string> sCHR;

class TabixList
{
	private:
		static const int N_CHR = 25;
	public:
		bool opened[N_CHR];
		Tabix ts[N_CHR];
		TabixList() {};
};

class snpinfo
{
	public:
		vector<string> chr;
		vector<int> pos;
		int min_pos;
		int max_pos;
		vector<string> rsid;
		vector<string> ref;
		vector<string> alt;
		vector<int> ld_index;
		int size();
		void pop(vector<int>&);
		void push(string&, int&, string&, string&, string&, int&);
		void push(string&, int&, string&, string&, string&);
		void print();
};

class hdata
{
	public:
		vector<vector<int>> iids;
		vector<vector<int>> hcts;
		vector<int> nhap;
		vector<vector<int>> map;
		void push_block (vector<int>, vector<int>, int) ;
		void push_map (vector<int>) ;
};

class ld_datum
{
	public:
		string id;
		string chr;
		int pos;
		string ref;
		string alt;
		vector<int> carriers;
		vector<bool> genotypes;
		int dir;
		int mac;
		int nhaps;
		int block;
		int ld_count;
		bool missing;
		bool fetched;
		ld_datum();
		ld_datum(string&, int&, string&, string&);
		void fetchGeno();
};

struct pair_int_hash
{
  size_t operator()(const pair<int,int>&x)const{
    return hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
  }
};

class ld_data
{
	public:
		int nhaps;
		int nsnps;
		int ichr;
		string schr;
		vector<ld_datum> pos_data;
		
		tsl::hopscotch_map<pair<int,int>,double,pair_int_hash> ld_vals;
		
		tsl::hopscotch_map<string, int> snp_index;
		ld_datum & operator [](int a) {return pos_data[a];};
		ld_data(int);
		void fetch(int&, string&, string&, int&);
		// void push(string, vector<int>, vector<bool>, int, int, int, bool);
		// void getGeno(int&, string&, string&);
		double gcorr(int, int);
		void print_LD( vector<int>& );
		Eigen::MatrixXd LD( vector<int>& );
		int nmiss( vector<int>& );
		vector<int> missing( vector<int>& );
		int nmiss();
};

/*
class ldref
{
	public:
		vector<unordered_map<string, gpair>>
		void push (vector<int>&, vector<bool>&, int&, int&, int&) ;
		void push_blank();
		double gcorr (int, int);
		void print_LD();
		Eigen::MatrixXd LD();
};*/

class ldref
{
	public:
		vector<ld_data> chr_data;
		// vector<vector<int>> carriers;
		// vector<vector<bool>> genotypes;
		// vector<int> dir;
		// vector<int> mac;
		// vector<int> block;
		// void push (vector<int>&, vector<bool>&, int&, int&, int&);
		// void push_blank();
		ld_data & operator [](int a) {return chr_data[a];};
		void print_LD(string&, vector<int>&);
		void print_LD(snpinfo&);
		Eigen::MatrixXd LD(string&, vector<int>&);
		Eigen::MatrixXd LD(snpinfo&);
		int nmiss(string&, vector<int>&);
		int nmiss(snpinfo&);
		vector<int> missing(string&, vector<int>&);
		vector<int> missing(snpinfo&);
		double gcorr (string&, int&, int&);
		// void print_LD();
		// Eigen::MatrixXd LD(string&, vector<int>&);
		void lazy(string&, int&, string&, string&, int&);
		void fetch(string&, int&, string&, string&, int&);
		void fetch(string&, int&);
};

//extern ldref LDREF;

void setJump (int);

void setMemoizeLD(bool);
void setPreload(bool);

string gsubstr(string , string , string );

vector<int> getRegion (string str);

int allele_check(string&, string&, string&, string&);

string asRegion (string chr, int pos, int end);

void initCHR();

void setPath(string&);

string subchr(string, string);

void getgeno(string vcf_path, int ichr, int ipos, string iref, string ialt, string irsid);

Eigen::MatrixXd extractLD(string vcf_path, snpinfo& kinfo, ostream& outstr);

void checkup();

double sizeIS(vector<int> &iid1, int &mac1, vector<bool> &genotypes2, int &n );

#endif
