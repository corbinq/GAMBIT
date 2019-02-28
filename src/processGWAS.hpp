#ifndef PROCESSGWAS_HPP
#define PROCESSGWAS_HPP

#include "eigenmvn.hpp"
#include "processLD.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <random>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

using namespace std;

class snpdata
{
	public:
		snpinfo info;
		
		vector<double> z;
		vector<double> n;
		vector<double> w;
		
		string anno_class;
		string anno_subclass;

		string gene_list;
		
		string flag;
		
		int nvar_total;
		int nvar;
		
		double zstat;
		double w_adj;
		double theta;
		double qstat;
		double logdet;
		double bf_partial;
		double bf;
		double pval;
		double hsq;
		double ns;
		
		int size();
		
		void pop(vector<int>&);
		void pop_missing();
		void print();
		void skat();
		void burden();
		void cauchy();
		void cauchy_w(bool);
		
		string print_range(bool);
};

class dparam
{
	public:
		double p;
		double tau;
};
typedef tsl::hopscotch_map<string, dparam> dparam_l;

class unitdata
{
	public:
		string gene_name;
		string chrom;
		int start_pos;
		int end_pos;
		
		bool omit_gene;
	
		tsl::hopscotch_map<string, snpdata> sstats;
		tsl::hopscotch_map<string, double> pvals;
		
		//tsl::hopscotch_map<string, snpdata> regel_sstats;
		
		double global_pval;
		double Lform_MMVN_pval;
		double min_pval;
		double naive_pval;
		
		double pval_bonf, pval_mmvn, pval_acat, pval_skat;
		
		string Lform_pvals;
		
		string top_class; 
		string top_subclass; 
		
		set<int> uniq_pos;
		
		dparam_l par;
		
		vector<string> groups;
		vector<string> tissues;
		vector<string> regels;
		
		vector<string> Lform_groups;
		
		int n_tiss_uniq;
		int n_uniq_tests;
		unitdata() {};
		snpdata & operator[](string g) {return sstats[g];};
		
		//void skat(string);
		//void burden(string);
		//void cauchy(string);
		//void cauchy_w(string,bool);
		
		void runTest(string, char);
		
		void multiBurden();
		void globalPval();
		
		string print_summary();
		string print_all_groups();
		
		snpinfo sinfo(string);
		void push(string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, double w, bool is_tissue, int& index);
		//void print_range(string);
		string print_range(string&);
		void update(dparam_l, vector<string>);
		void Lform_covar( double&, string&, string& );
		Eigen::MatrixXd get_Lform_covar();
		Eigen::MatrixXd get_Lform_covar(vector<double>&);
		bool exists(string);
};


class listLD
{
	public:
		string gene1;
		string gene2;
		vector<string> subclass1;
		vector<string> subclass2;
		Eigen::MatrixXd LDmat;
};

class chrom_index 
{
	public:
		vector<int> start;
		vector<int> end;
	
		chrom_index();
		
		void init();
		void update(string&, int);
		int startIndex(string&);
		int endIndex(string&);
};

class tssdat
{
	public:
		int k_s, k_e, m;
		
		chrom_index chr_idx;
		
		vector<string> genes;
		vector<string> chr;
		vector<int> start;
		vector<int> end;

		tssdat();
		tssdat(string&, string);
		void readTSS(string&, string);
};

class regel
{
	public:
		int k, m;
		
		chrom_index chr_idx;
		
		vector<string> names;
		vector<string> chr;
		vector<int> start;
		vector<int> end;
		vector<string> elid;
		vector<map<string, double>> genes;
		map<string, string> gene_maps;
		map<string, vector<string>> gene_names;
		vector<map<string, double>> groups;
		map<string, string> group_maps;
		regel();
		regel(string&, string);
		void readElements(string&, string);
};

class eweight
{
	public:
		tsl::hopscotch_set<string> tissues;
		int k, m;
		vector<string> chr;
		
		chrom_index chr_idx;
		
		vector<int> pos;
		vector<string> ref;
		vector<string> alt;
		vector<vector<string>> tissue;
		vector<vector<string>> gene;
		vector<vector<double>> beta;
		eweight() {};
		eweight(string&, string, int);
		eweight(string&, string&, string, int);
		void readBetas(string&, string, int);
		void readBetas(string&, string, string, int);
		// void next();
};

class gwasdata
{
	public:
		tsl::hopscotch_map<string, unitdata> data;
		
		regel regel_data;
		eweight eweight_data;
		tssdat tss_data;
		
		dparam_l omega;
		double p_c;
		vector<string> elements;
		vector<string> groups;
		vector<string> tissues;
		vector<string> genes;
		unitdata & operator[](string& g) {return data[g];};
		
		bool cross_regel(string&, int&, string&, string&, string&, double&, double&, int&, set<string>&);
		bool cross_eweight(string&, int&, string&, string&, string&, double&, double&, int&, set<string>&);
		bool cross_tss(string&, int&, string&, string&, string&, double&, double&, int&, set<string>&);
		
		void push(string& gene, string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, int& index);
		void push(string& gene, string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, double w, int& index, bool is_tissue);
		void pushElement(string& name, string& elid, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, int& index);
		
		void runTest(string& , string&, bool, char );
		double EM_round(int);
};

class annodef
{
	public:
		bool nodef;
		tsl::hopscotch_map<string, vector<string>> defs;
		tsl::hopscotch_map<string, string> meta;
		vector<string> & operator[](string a) {return defs[a];};
		bool defn(string);
		void push(string, string);
};

string pretty(int n);

void excludeMissing();

void debug();

void setAlphaTSS(vector<double>);
void setWindowSizeTSS(int);
void setVerbosityTSS(double);
void setMultiForm(string);


bool read_snp_info(string &file_path, snpinfo &sinfo);

bool read_gwas(string& zpath, string& region, gwasdata &gwas, annodef &adef, int& fetch_mode);

bool read_defs(string &file_path, annodef &adef);

bool read_sstats(string &file_path, unitdata &data);

double skat_pval(vector<double>& Z, Eigen::MatrixXd& V, double& logdet);

#endif
