#include "processGWAS.hpp"

#include "geneBasedTests.hpp"
#include "common.hpp"

using namespace std;

Tabix zfile;

double tau_0 = 0.001;
double tau_lower = 0.000000000001;
double tau_upper = 1;
double theta_lower = 0.00000000001;
double p_0 = 0.005;
double p_lower = 1e-4;
double p_upper = 1;
double epsilon = 0.000000000001;

int N_Z_ZERO = 0;

double MIN_NZ_Z = 1;
double MIN_NZ_QT = 1;
random_device R_DEV;
mt19937 R_GEN(R_DEV());

uniform_real_distribution<double> NZ_RUNIF;

vector<double> TSS_ALPHA;
int TSS_WINDOW;
double TSS_VERBOSITY_PVAL;

bool DEBUG = false;
bool NOMISS = false;

bool ENH_CAUCHY = true;

bool PRINT_ALL_REGEL_GENES = true;

bool GET_GENE_SETS = false;

map<string, set<string>> gene_sets;

int JUMP = 100000;

annodef ANNO_DEFS;

void setDistNZ(){
	MIN_NZ_QT = 1 - pchisq(MIN_NZ_Z*MIN_NZ_Z, 1);
	cerr << "minimum non-zero abs(z-score) = " << MIN_NZ_Z << " (qt = " << MIN_NZ_QT << ")\n";
	cerr << N_Z_ZERO << " variants with abs(z-score) == 0 (assumed censored due to precision) \n\n";
	NZ_RUNIF = uniform_real_distribution<double>(0.00, MIN_NZ_QT );
}

double simChiNZ(){
	//double u = NZ_RUNIF(R_GEN);
	double out = qchisq(
		NZ_RUNIF(R_GEN), 1.00
	);
	//cout << u << "\t" << out << "\n";
	return out;
}

void debug() {
	DEBUG = true;
}

void setAlphaTSS(vector<double> a){
	TSS_ALPHA = a;
}

void setWindowSizeTSS(int window){
	TSS_WINDOW = window;
}

void setVerbosityTSS(double pv){
	TSS_VERBOSITY_PVAL = pv;
}

void excludeMissing() {
	NOMISS = true;
}

void updateGeneSets(vector<string>& genes){
	if ( GET_GENE_SETS ){
		for(string& a : genes){
			for(string& b : genes){
				gene_sets[a].insert(b);
			}
		}
	}
}

void updateGeneSets(set<string>& genes){
	if ( GET_GENE_SETS ){
		for(const auto& a : genes){
			for(const auto& b : genes){
				gene_sets[a].insert(b);
			}
		}
	}
}

void printGeneSets(){
	if ( GET_GENE_SETS ){
		for ( const auto& gset : gene_sets ){
			cout << gset.first << "\t";
			for( const auto& gene : gset.second ){
				cout << gene << ",";
			}
			cout << "\n";
		}
	}
}

string pretty(int n) {
	string ns = to_string(n);
	int i = ns.length() - 3;
	while (i > 0) {
		ns.insert(i, ",");
		i -= 3;
	}
	return ns;
}

string scien(double& x, int n = 4){
	ostringstream out;
	out << scientific;
	out.precision(n);
	out << x;
	return (string) out.str();
}


string round_to_string(double& x, int n = 4){
	ostringstream out;
	out.precision(n);
	out << x;
	return (string) out.str();
}

void processAnnoFlags(string& line, map<string,string>& subclass_map, bool& map_mode){
	if( line[0] == '#' ) {
		size_t p_dlim = line.find("=");
		string tflag = line.substr(0, p_dlim);
		line = line.substr(p_dlim+1);
		tflag.erase(remove(tflag.begin(), tflag.end(), ' '), tflag.end());
		tflag.erase(remove(tflag.begin(), tflag.end(), '#'), tflag.end());
		if( tflag=="SUBCLASS_IDS" || tflag=="SUBCLASS_KEYS" || tflag=="TISSUE_IDS" || tflag=="TISSUE_KEYS" ){
			map_mode = true;
			string parse_tis;
			istringstream iss_gene(line);
			while( getline(iss_gene, parse_tis, ',') ) {
				p_dlim = parse_tis.find(":");
				string key = parse_tis.substr(0, p_dlim);
				parse_tis = parse_tis.substr(p_dlim+1);
				subclass_map[key] = parse_tis;
			}

		}
	}
}

void pruneCovar(Eigen::MatrixXd& cmat, int prune_id)
{
	unsigned int m = cmat.rows();
	
	if( m != cmat.cols() ){
		cerr << "\n\nFATAL ERROR:: pruneCovar called on non-square matrix ... \n\n";
		abort();
	}
	if( m < prune_id ){
		cerr << "\n\nFATAL ERROR:: invalid prune_id supplied to pruneCovar ... \n\n";
		abort();
	}

    cmat.block(prune_id,0,m-prune_id,m) = cmat.bottomRows(m-prune_id);
	cmat.block(0,prune_id,m,m-prune_id) = cmat.rightCols(m-prune_id);

    cmat.conservativeResize(m-1,m-1);
	
	return;
}

void pruneLD(Eigen::MatrixXd& cmat, double thres = 1){
	unsigned int m = cmat.rows();
	if( m > 1 ){
		for( int i = 0; i < m-1; i++ ){
			for( int j = i+1; j < m; j++ ){
				if( abs(cmat(i,j)) >= thres ){
					pruneCovar(cmat, j);
					pruneLD(cmat);
					return;
				}
			}
		}
	}
}

void pruneLD(Eigen::MatrixXd& cmat, vector<double>& x, double thres = 1){
	unsigned int m = cmat.rows();
	if( x.size() != m ){
		cerr << "FATAL ERROR: pruneLD called on vector and matrix of differing sizes...\n";
		abort();
	}
	if( m > 1 ){
		for( int i = 0; i < m-1; i++ ){
			for( int j = i+1; j < m; j++ ){
				if( abs(cmat(i,j)) >= thres ){
					x.erase( x.begin() + j );
					pruneCovar(cmat, j);
					pruneLD(cmat, x);
					return;
				}
			}
		}
	}
}

bool read_snp_info(string &file_path, snpinfo &sinfo) {
	string line;
	ifstream file_stream(file_path);
	if( file_stream.is_open() ) {
		while ( getline(file_stream, line) ) {
			string chr;
			int pos;
			string ref, alt, rsid;
			istringstream iss(line);
			iss >> chr >> pos >> rsid >> ref >> alt;
			sinfo.push(chr, pos, ref, alt, rsid);
		}
		file_stream.close();
		return true;
	}
	return false;
}


void snpdata::pop(vector<int>& rm){
	info.pop(rm);
	// Note: rm elements are in reverse order
	for(int& i : rm){
		n.erase(n.begin() + i);
		z.erase(z.begin() + i);
		w.erase(w.begin() + i);
	}
}


void snpdata::pop_missing(){
	if( z.size() > 1 ){
		vector<int> to_pop = LDREF.missing(info);
		pop(to_pop);
	}
}

tssdat::tssdat(){
	m = 0;
	k_s = 0;
	k_e = 0;
}


tssdat::tssdat(string& file_path, string region) {
	readTSS(file_path, region);
}


bool gwasdata::cross_tss(string& chr_i, int& pos_i, string& ref_i, string& alt_i, string& rsid_i, double& n, double& z, int& index, set<string>& all_genes) {
	if( tss_data.chr.size() == 0 ) return false;
	
	vector<string>& chr = tss_data.chr;
	vector<int>& start = tss_data.start;
	vector<int>& end = tss_data.end;
	vector<string>& tss_genes = tss_data.genes;
	
	int& k_s = tss_data.k_s;
	int& k_e = tss_data.k_e;
	int& m = tss_data.m;
	
	string out_group = "TSS";
	
	bool out = false;
	bool contin = true;
	//if( iCHR[chr_i] != iCHR[chr[k_s]] ){
	//		cerr << chr_i << "\t" << pos_i << "\n";
	//		cerr << chr[k_s] << "\t" << start[k_s] << "\t" << end[k_s] << "\n";
	//}
	while( iCHR[chr_i] != iCHR[chr[k_s]] && contin ) {
		if( iCHR[chr_i] > iCHR[chr[k_s]] ) {
			if( k_s < m - 1 ) {
				k_s++;
			}else {
				contin = false;
			}
		}
		if( iCHR[chr_i] < iCHR[chr[k_s]] ) {
			//if( k > 0 ) {
			//	k--;
			//}
			//else {
				contin = false;
			//}
		}
		k_e = k_s;
	}
	
	int k = k_s;
	
	while( iCHR[chr_i] == iCHR[chr[k]] && contin ) {
		if( (pos_i >= start[k] && pos_i <= end[k]) || min(abs(start[k]-pos_i), abs(end[k]-pos_i)) <= TSS_WINDOW  ) {
			out = true;
			all_genes.insert(tss_genes[k]);
			
			double tss_weight = min(abs(start[k]-pos_i), abs(end[k]-pos_i));
			if(pos_i >= start[k] && pos_i <= end[k]){
				tss_weight = 1;
			}
			
			push(tss_genes[k], out_group, chr_i, pos_i, ref_i, alt_i, rsid_i, n, z, tss_weight, index, false);

			k++;
		}else if( pos_i > end[k] + TSS_WINDOW ) {
			if( k < m - 1 ) {
				k++;
				k_s = k;
			}
			else {
				contin = false;
			}
		}else if( pos_i < start[k] - TSS_WINDOW ) {
			k_e = k - 1;
			contin = false;
		}else{
			cerr << chr_i << "\t" << pos_i << "\n";
			cerr << chr[k] << "\t" << start[k] << "\t" << end[k] << "\n";
		}
	}
	return out;
}


void tssdat::readTSS(string& file_path, string region) {
	
	Tabix tmp(file_path);
	tmp.setRegion(region);
	
	int istart, iend;
	string line, ichr, igene;
	while( tmp.getNextLine(line) ) {
		if( line[0] != '#' && line.length() > 0 ) {
			istringstream iss(line);
			iss >> ichr >> istart >> iend >> igene;

			chr.push_back(ichr);
			start.push_back(istart);
			end.push_back(iend);
			genes.push_back(igene);

			m++;
			if( m % 500 == 0 ) {
				cerr << "\rprocessed " << pretty(m) << " genes from TSS bed file...";
			}
		}
	}
	cerr << "\rprocessed " << pretty(m) << " genes from TSS bed file...\n\n";
}

regel::regel(){
	m = 0;
}

bool gwasdata::cross_regel(string& chr_i, int& pos_i, string& ref_i, string& alt_i, string& rsid_i, double& n, double& z, int& index, set<string>& all_genes) {
	if( regel_data.chr.size() == 0 ) return false;
	
	vector<string>& chr = regel_data.chr;
	vector<int>& start = regel_data.start;
	vector<int>& end = regel_data.end;
	
	vector<string>& names = regel_data.names;
	vector<string>& elid = regel_data.elid;
	map<string, vector<string>>& gene_names = regel_data.gene_names;
	
	vector<map<string,double>>& gene_w = regel_data.genes;
	
	int& k = regel_data.k;
	int& m = regel_data.m;
	
	bool out = false;
	bool contin = true;
	while( iCHR[chr_i] != iCHR[chr[k]] && contin ) {
		if( iCHR[chr_i] > iCHR[chr[k]] ) {
			if( k < m - 1 ) {
				k++;
			}
			else {
				contin = false;
			}
		}
		if( iCHR[chr_i] < iCHR[chr[k]] ) {
			//if( k > 0 ) {
			//	k--;
			//}
			//else {
				contin = false;
			//}
		}
	}
	while( iCHR[chr_i] == iCHR[chr[k]] && contin ) {
		if( pos_i >= start[k] && pos_i <= end[k] ) {
			out = true;
			
			for(string& gene : gene_names[names[k]] ){
				all_genes.insert(gene);
			}
			
			//if( ENH_CAUCHY ){
			
			//}else{
				pushElement(names[k], elid[k], chr_i, pos_i, ref_i, alt_i, rsid_i, n, z, index);
			//}
			
			if( k < m - 1 ) {
				int j = k;
				while( j < m - 1 && iCHR[chr_i] == iCHR[chr[j]] && start[j] <= pos_i ) {
					j++;
					if( pos_i >= start[j] && pos_i <= end[j] ) {
						//if( ENH_CAUCHY ){
			
						//}else{
							pushElement(names[j], elid[j], chr_i, pos_i, ref_i, alt_i, rsid_i, n, z, index);
						//}
					}
				}
			}
			if( k > 0 ) {
				int j = k;
				while( j > 0 && iCHR[chr_i] == iCHR[chr[j]] && end[j] >= pos_i ) {
					j--;
					if( pos_i >= start[j] && pos_i <= end[j] ) {
						//if( ENH_CAUCHY ){
			
						//}else{
							pushElement(names[j], elid[j], chr_i, pos_i, ref_i, alt_i, rsid_i, n, z, index);
						//}
					}
				}
			}
			contin = false;
		}
		else if( pos_i > end[k] ) {
			if( k < m - 1 ) {
				k++;
			}
			else {
				contin = false;
			}
		}
		else if( pos_i < start[k] ) {
			contin = false;
		}
	}
	return out;
}


void regel::readElements(string& file_path, string region) {
	k = 0;
	m = 0;
	Tabix tmp(file_path);
	tmp.setRegion(region);
	
	bool map_mode = false;
	map<string, string> tissue_map;
	
	int istart, iend;
	string line, ichr, ename, ielid, gene_map, group_map;
	while( tmp.getNextLine(line) ) {
		if( line[0] == '#' ){
			bool last_map_mode = map_mode;
			processAnnoFlags(line, tissue_map, map_mode);
			if( map_mode && !last_map_mode ){
				cerr << "processed annotation keys for regulatory elements... \n\n";
			}
		}else if( line[0] != '#' && line.length() > 0 ) {
			istringstream iss(line);
			// cerr << line << "\n\n"
			iss >> ichr >> istart >> iend >> ename >> ielid >> gene_map >> group_map;
			// cerr <<  ichr << "\t" <<  istart << "\t" <<  iend << "\t" <<  ielid << "\t" <<  gene_map << "\t" <<  group_map << "\n\n";
			map<string, double> igenes;
			vector<string> igene_names;
			map<string, double> igroups;
			vector<string> igroup_names;
			string parse_gene;
			gene_maps[ielid] = gene_map;
			istringstream iss_gene(gene_map);
			while( getline(iss_gene, parse_gene, '|') ) {
				size_t pos_t = parse_gene.find(":");
				double iweight = stod(parse_gene.substr(pos_t+1));
				string igene = parse_gene.substr(0, pos_t);
				igenes[igene] = iweight;
				igene_names.push_back(igene);
				//cout << igene << "\n";
			}
			
			string parse_group;
			istringstream iss_group(group_map);
			while( getline(iss_group, parse_group, '|') ) {
				size_t pos_t = parse_group.find(":");
				double iweight = stod(parse_group.substr(pos_t+1));
				string igroup = parse_group.substr(0, pos_t);
				if(map_mode){
					igroup = tissue_map[igroup];
					igroup_names.push_back(igroup);
				}
				igroups[igroup] = iweight;
			}
			if( map_mode ){
				group_map = "";
				for( string& group : igroup_names ){
					if( group_map != "" ){
						group_map += ',';
					}
					group_map += (group + ":" + round_to_string(igroups[group],1)); 
				}
			}
			group_maps[ielid] = group_map;
			
			chr.push_back(ichr);
			start.push_back(istart);
			end.push_back(iend);
			names.push_back(ename);
			elid.push_back(ielid);
			genes.push_back(igenes);
			gene_names[ielid] = igene_names;
			groups.push_back(igroups);
			m++;
			if( m % 250 == 0 ) {
				cerr << "\rprocessed " << pretty(m) << " regulatory elements from bed file...";
			}
		}
	}
	cerr << "\rprocessed " << pretty(m) << " regulatory elements from bed file...\n\n";
}


regel::regel(string& file_path, string region) {
	readElements(file_path, region);
}


bool gwasdata::cross_eweight(string& chr_i, int& pos_i, string& ref_i, string& alt_i, string& rsid_i, double& n, double& z, int& index, set<string>& all_genes) {
	if( eweight_data.chr.size() == 0 ) return false;

	vector<string>& chr = eweight_data.chr;
	vector<int>& pos = eweight_data.pos;
	vector<vector<string>>& gene = eweight_data.gene;
	int& k = eweight_data.k;
	int& m = eweight_data.m;
	
	bool out = false;
	bool contin = true;
	while( iCHR[chr_i] > iCHR[chr[k]] && k < chr.size()-1 ) {
		k++;
	}
	while( iCHR[chr_i] >= iCHR[chr[k]] && pos_i >= pos[k] && contin ) {
		if( iCHR[chr_i] >= iCHR[chr[k]] && pos_i == pos[k] ) {
			int dir = allele_check(ref_i, alt_i, eweight_data.ref[k], eweight_data.alt[k]);
			if( abs(dir) > 0.0 ) {
				for( int j=0; j < gene[k].size(); j++) {
					double dbeta = eweight_data.beta[k][j] * dir;
					push(gene[k][j], eweight_data.tissue[k][j], chr_i, pos_i, ref_i, alt_i, rsid_i, n, z, dbeta, index, true);
					all_genes.insert(gene[k][j]);
				}
				out = true;
			}
		}
		if( k < chr.size()-1 ) {
			k++;
		}
		else {
			contin = false;
		}
	}
	return out;
}


void eweight::readBetas(string& file_path, string region, int tmerge) {
	k = 0;
	m = 0;
	
	map<string, string> tissue_map;
	bool map_mode = false;
	
	Tabix tmp(file_path);
	tmp.setRegion(region);
	string line;
	string ichr;
	int ipos;
	string irsid, iref, ialt, mapping;
	if(tmerge) tissues.insert("eSNPs");
	while( tmp.getNextLine(line) ) {
		if( line[0] == '#' ){
			bool last_map_mode = map_mode;
			processAnnoFlags(line, tissue_map, map_mode);
			if( map_mode && !last_map_mode ){
				cerr << "processed eSNP tissue keys ... \n\n";
			}
		}else if( line[0] != '#' && line.length() > 0 ) {
			istringstream iss(line);
			iss >> ichr >> ipos >> irsid  >> iref >> ialt >> mapping;
			// cerr << ichr << "\t" << ipos << "\t" << iref << "\t" << ialt << "\n";
			vector<string> itissues;
			vector<string> igenes;
			vector<double> ibetas;
			string parse_gene;
			istringstream iss_gene(mapping);
			while( getline(iss_gene, parse_gene, '|') ) {
				size_t pos_g = parse_gene.find("=");
				string igene = parse_gene.substr(0, pos_g);
				// cerr << "\t" << igene << "\n";
				parse_gene = parse_gene.substr(pos_g+1);

				double tbeta = 0;
				int nt = 0;

				string parse_tissue;
				istringstream iss_tissue(parse_gene);
				while( getline(iss_tissue, parse_tissue, ';') ) {
					size_t pos_t = parse_tissue.find("@");
					double ibeta = stod(parse_tissue.substr(0, pos_t));
					tbeta += ibeta;
					nt++;
					string itissue = parse_tissue.substr(pos_t+1);
					if( map_mode ){
						itissue = tissue_map[itissue];
					}
					// cerr << "\t\t" << itissue << "\t" << to_string(ibeta) << "\n";
					if( !tmerge ) {
						itissues.push_back(itissue);
						igenes.push_back(igene);
						ibetas.push_back(ibeta);
						tissues.insert(itissue);
					}
				}
				if( tmerge ) {
					itissues.push_back("eSNPs");
					igenes.push_back(igene);
					ibetas.push_back(tbeta/((double) nt));
				}
			}
			chr.push_back(ichr);
			pos.push_back(ipos);
			ref.push_back(iref);
			alt.push_back(ialt);
			tissue.push_back(itissues);
			gene.push_back(igenes);
			beta.push_back(ibetas);
			m++;
			if( m % 1000 == 0 ) {
				cerr << "\rprocessed " << pretty(m) << " eSNPs across "<< pretty(tissues.size()) <<" tissues...";
			}
		}
	}
	cerr << "\rprocessed " << pretty(m) << " eSNPs across "<< pretty(tissues.size()) <<" tissues...\n\n";
}


void eweight::readBetas(string& file_path, string region, string tlist_path, int tmerge) {
	string line, tis;
	int nt = 0;
	
	map<string, string> tissue_map;
	bool map_mode = false;
	
	if( tlist_path.find(".") != string::npos){
		ifstream tlist_stream(tlist_path);
		while ( getline(tlist_stream, tis) ) {
			tissues.emplace(tis.substr(0, tis.find("\t")));
			nt++;
		}
		tlist_stream.close();
	}else{
		istringstream tlist_stream(tlist_path);
		while( getline(tlist_stream, tis, ',') ){
			tissues.emplace(tis);
			nt++;
		}		
	}
	//cerr << "Retaining eSNPs for " << nt << " specified tissues  ... \n\n";
	k = 0;
	m = 0;
	Tabix tmp(file_path);
	tmp.setRegion(region);
	string ichr;
	int ipos;
	string iref, ialt, irsid, mapping;
	while( tmp.getNextLine(line) ) {
		if( line[0] == '#' ){
			bool last_map_mode = map_mode;
			processAnnoFlags(line, tissue_map, map_mode);
			if( map_mode && !last_map_mode ){
				cerr << "processed eSNP tissue keys ... \n\n";
			}
		}else if( line[0] != '#' && line.length() > 0 ) {
			istringstream iss(line);
			iss >> ichr >> ipos >> irsid >> iref >> ialt >> mapping;
			bool retain = false;
			vector<string> itissues;
			vector<string> igenes;
			vector<double> ibetas;
			string parse_gene;
			istringstream iss_gene(mapping);
			while( getline(iss_gene, parse_gene, '|') ) {
				size_t pos_g = parse_gene.find("=");
				string igene = parse_gene.substr(0, pos_g);
				parse_gene = parse_gene.substr(pos_g+1);

				double tbeta = 0;
				int nt = 0;

				bool retain_gene = false;

				string parse_tissue;
				istringstream iss_tissue(parse_gene);
				while( getline(iss_tissue, parse_tissue, ';') ) {
					size_t pos_t = parse_tissue.find("@");
					double ibeta = stod(parse_tissue.substr(0, pos_t));
					string itissue = parse_tissue.substr(pos_t+1);
					if( map_mode ){
						itissue = tissue_map[itissue];
					}
					//cerr << "\n" << itissue << "\n";
					if( tissues.find(itissue) != tissues.end() ) {
						retain = true;
						retain_gene = true;
						if( !tmerge ) {
							itissues.push_back(itissue);
							igenes.push_back(igene);
							ibetas.push_back(ibeta);
						}
						tbeta += ibeta;
						nt++;
					}
				}
				if( tmerge && retain_gene ) {
					itissues.push_back("eQTLs");
					igenes.push_back(igene);
					ibetas.push_back(tbeta/((double) nt));
				}
			}
			if( retain ) {
				chr.push_back(ichr);
				pos.push_back(ipos);
				ref.push_back(iref);
				alt.push_back(ialt);
				tissue.push_back(itissues);
				gene.push_back(igenes);
				beta.push_back(ibetas);
				m++;
				if( m % 1000 == 0 ) {
					cerr << "\rprocessed " << pretty(m) << " eSNPs across "<< pretty(tissues.size()) <<" tissues...";
				}
			}
		}
	}
	cerr << "\rprocessed " << pretty(m) << " eSNPs across "<< pretty(tissues.size()) <<" tissues...\n\n";
}


eweight::eweight(string& file_path, string region, int tmerge) {
	readBetas(file_path, region, tmerge);
}


eweight::eweight(string& file_path, string& tlist_path, string region, int tmerge) {
	readBetas(file_path, tlist_path, region, tmerge);
}


// void eweights::open(string tlist_path) {
// string line, tissue, tissue_path;
// ifstream tlist_stream(tlist_path);
// int nt = 0.0;
// if( tlist_stream.is_open() ) {
// while ( getline(tlist_stream, line) ) {
// istringstream iss(line);
// iss >> tissue >> tissue_path;
// eweight tmp(tissue, tissue_path);
// tstreams[tissue] = tmp;
// tissues.push_back(tissue);
// // cout << tissue << "\n";
// nt++;
// }
// tlist_stream.close();
// }
// // cerr << "\n--------------------------------------------------------------------\n\n";
// cerr << "  -  done processing " << pretty(nt) << " tissues \n\n";
// cerr << "  -  beginning to read GWAS data\n\n";
// }

bool processAnno(string &in, vector<string> &annos, vector<string> &genes) {
	if( in.find(':') == string::npos ) {
		return false;
	}
	istringstream iss(in);
	string line;
	while(getline(iss, line, ',')){
		// cout << line << "\n";
		if( line.find(':') != string::npos ){
			string gene;
			string anno = line.substr(0, line.find(":"));
			line = line.substr(line.find(":") + 1);
			istringstream lss(line);
			// cout << anno << "\n";
			while ( getline(lss, gene, '|') ) {
				// cout << gene << "\n";
				annos.push_back(anno);
				genes.push_back(gene);
			}
		}
	}
	return true;
}


void annodef::push(string child, string group) {
	if( defs.find(child) == defs.end() ) {
		defs[child] = vector<string>();
	}
	defs[child].push_back(group);
	// cout << child << " :: ";
	// for( int i = 0; i < defs[child].size(); i++){
	// cout << defs[child][i] << ", ";
	// }
	// cout << "\n";
}


bool annodef::defn(string child) {
	if( nodef ){
		return false;
	}else if( defs.find(child) == defs.end() ) {
		return false;
	}
	return true;
}


bool read_defs(string &file_path, annodef &adef) {
	string line;
	if( file_path == "" ){
		adef.nodef = true;
		return true;
	}
	adef.nodef = false;
	ifstream file_stream(file_path);
	if( file_stream.is_open() ) {
		while ( getline(file_stream, line) ) {
			string meta, group, child;
			istringstream iss(line);
			getline(iss, meta, '\t');
			getline(iss, group, '\t');
			adef.meta[group] = meta;
			while ( getline(iss, child, ',') ) {
				adef.push(child, group);
			}
		}
		ANNO_DEFS = adef;
		return true;
	}
	return false;
}



bool fetchLine(Tabix& tfile, string& line, vector<string>& regions, int& k, int& fetch_mode){
	if( fetch_mode ){
		if( k > 0 ){
			if( tfile.getNextLine(line) ){
				return true;
			}else if( k < regions.size() ){
				tfile.setRegion( regions[k] );
				k++;
				return fetchLine(tfile, line, regions, k, fetch_mode);
			}else{
				return false;
			}
		}else{
			tfile.setRegion( regions[k] );
			k++;
			return fetchLine(tfile, line, regions, k, fetch_mode);
		}
	}else{
		return tfile.getNextLine(line);
	}
}


void printProgress(int& nsnps, int& neqtl, int& nrege, int ngene, bool& use_esnps, bool& use_rsnps){
	if( use_esnps && use_rsnps ){
		cerr << "\rprocessed " << pretty(nsnps) << " SNPs (" << pretty(neqtl) << " eSNPs, " << pretty(nrege) << " in regulatory elements) and " << pretty(ngene) << " genes";
	}else if( use_esnps ){
		cerr << "\rprocessed " << pretty(nsnps) << " SNPs (" << pretty(neqtl) << " eSNPs) and " << pretty(ngene) << " genes";
	}else if( use_rsnps ){
		cerr << "\rprocessed " << pretty(nsnps) << " SNPs (" << pretty(nrege) << " in regulatory elements) and " << pretty(ngene) << " genes";
	}else{
		cerr << "\rprocessed " << pretty(nsnps) << " SNPs and " << pretty(ngene) << " genes";
	}
}

bool read_gwas(string& zpath, string& region, gwasdata &gwas, annodef &adef, int &fetch_mode) {
	string line;
	zfile.open(zpath);
	zfile.setRegion(region);
	// ifstream file_stream(file_path);
	// if( file_stream.is_open() ){
	int nsnps = 0;
	int nrege = 0;
	int nkept = 0;
	int neqtl = 0;
	int ngenes = 0;
	int npassed = 0;
	int wk = 0;
	vector<string> regions;
	
	bool use_esnps = (gwas.eweight_data.chr.size() > 0);
	bool use_rsnps = (gwas.regel_data.chr.size() > 0);
	
	if( fetch_mode ){
		string chr_s = gwas.eweight_data.chr[0];
		int pos_s = gwas.eweight_data.pos[0];
		int pos_e = gwas.eweight_data.pos[0];
		for( int i = 1; i < gwas.eweight_data.chr.size(); i++ ){
			if( gwas.eweight_data.pos[i] - pos_e > JUMP || gwas.eweight_data.chr[i] != chr_s ){
				string region = chr_s + ":" + 
					to_string( 
						(pos_s > 1) ? pos_s-1 : 0
					) + "-" + 
					to_string(pos_e+1);
				regions.push_back(region);
				chr_s = gwas.eweight_data.chr[i];
				pos_s = gwas.eweight_data.pos[i];
				pos_e = gwas.eweight_data.pos[i];
			}else{
				pos_e = gwas.eweight_data.pos[i];
			}
		}
		regions.push_back( chr_s + ":" + to_string(pos_s-1) + "-" + to_string(pos_e+1) );
	}
	while ( fetchLine(zfile, line, regions, npassed, fetch_mode) ) {
		vector<string> genes;
		set<string> gene_set_genes;
		vector<string> groups;
		// while ( getline(file_stream, line) ){
		if( line[0] != '#' && line.length() > 0 ) {
			string chr;
			int pos;
			string ref, alt, rsid, anno;
			double n, z, w;
			
			if( abs(z) > 0 ){
				MIN_NZ_Z = min(MIN_NZ_Z, abs(z));
			}else{
				N_Z_ZERO++;
			}
			
			istringstream iss(line);
			iss >> chr >> pos >> ref >> alt >> rsid >> n >> z >> anno;
			// cout << "MODE" << fetch_mode << "\t" << chr << ":" << pos << "\t" << regions[npassed-1] << "\n";
			// cout << rsid << "\t" << anno << "\n";
			bool kept = false;
			int index = 0;
			if( gwas.cross_eweight(chr, pos, ref, alt, rsid, n, z, index, gene_set_genes) ) {
				kept = true;
				neqtl++;
			}
			if( gwas.cross_regel(chr, pos, ref, alt, rsid, n, z, index, gene_set_genes) ) {
				kept = true;
				nrege++;
			}
			if( gwas.cross_tss(chr, pos, ref, alt, rsid, n, z, index, gene_set_genes) ) {
				kept = true;
			}
			if( !adef.nodef ){
				genes.clear();
				groups.clear();
				if( processAnno(anno, groups, genes) ) {
					// cerr << groups.size() << "\t" << genes.size() << "\n";
					for( int i = 0; i < groups.size(); i++ ){
						if( adef.defn(groups[i]) ) {
							kept = true;
							for( int j = 0; j < adef[groups[i]].size(); j++) {
								gwas.push(genes[i], adef[groups[i]][j], chr, pos, ref, alt, rsid, n, z, index);
							}
						}
						gene_set_genes.insert(genes[i]);
					}
				}
			}
			updateGeneSets(gene_set_genes);
			if( kept ) {
				nkept++;
			}
			nsnps++;

			if( nsnps % 1000 == 0 ) {
				printProgress(nsnps, neqtl, nrege, gwas.genes.size(), use_esnps, use_rsnps);
			}
		}
	}
	zfile.close();
	
	//printGeneSets();
	
	vector<string> omitted_genes;
	vector<string> all_genes = gwas.genes;
	for( int j = all_genes.size()-1; j >= 0; j-- ){
		if( gwas[all_genes[j]].omit_gene ){
			omitted_genes.push_back( all_genes[j] );
			gwas.genes.erase( gwas.genes.begin() + j );
			gwas.data.erase(all_genes[j]);
		}
	}
	
	printProgress(nsnps, neqtl, nrege, gwas.genes.size(), use_esnps, use_rsnps);
	cerr << "\n\n";
	
	setDistNZ();
	
	if( omitted_genes.size() > 0 ){
		
		cerr << "omitted " << pretty(omitted_genes.size()) << " genes due to inconsistent chromosome: ";
		
		for( string& gene : omitted_genes ){
			cerr << gene << ",";
		}

		cerr << "\n\n";
		
	}
	
	
	double n_groups = gwas.groups.size();
	double n_tissues = gwas.tissues.size();
	double n_annos = n_groups - n_tissues;

	cerr << "  -  done processing GWAS data and annotations\n\n";
	cerr << "  -  retained " << pretty(nkept) << " out of " << pretty(nsnps) << " SNPs  \n\n";
	cerr << "  -  starting analysis of " << pretty(gwas.genes.size()) << " genes\n\n";
	cerr << "  -  estimating LD from " << subchr(vcf_path, "*") << "\n\n";
	cerr << "  -  " << pretty(n_groups) << " total functional groups (" << pretty(n_tissues) << " tissues)\n\n";

	for( int i = 0; i < n_groups; i++ ) {
		gwas.omega[ gwas.groups[i] ].p = 1.0/(1.0 + n_annos);
		gwas.omega[ gwas.groups[i] ].tau = tau_0;
	}
	for( int i = 0; i < n_tissues; i++ ) {
		gwas.omega[ gwas.tissues[i] ].p *= 1.0/n_tissues;
	}

	gwas.p_c = p_0;

	for( int i = 0; i < gwas.genes.size(); i++ ) {
		gwas[ gwas.genes[i] ].par = gwas.omega;
	}

	return true;
}


snpinfo unitdata::sinfo(string group) {
	if( sstats.find(group) == sstats.end() ) {
		snpdata temp;
		sstats[group] = temp;
	}
	return sstats[group].info;
}


string snpdata::print_range(bool bmode){
	string out = info.chr[0] + "\t" + to_string(info.pos[0]);
	if( bmode ){
		out += "\t" + to_string(info.pos[info.pos.size()-1]);
	}else if( info.pos.size() >= 2){
		if( info.pos.size() >= 2){
			if( info.pos.size() == 2 ){
				out += "," + to_string(info.pos[info.pos.size()-1]);
			}else{
				out += "-" + to_string(info.pos[info.pos.size()-1]);
			}
		}
	}
	return out;
}

// NEEDS_MOD
// void unitdata::print_range(string group) {
	// if( sstats.find(group) != sstats.end() ) {
		// sstats[group].print_range(false);
	// }
// }


string unitdata::print_range(string& group){
	if( sstats.find(group) != sstats.end() ) {
		return sstats[group].print_range(true);
	}else{
		return "";
	}
}


double Lform_cov(snpdata& sstats1, snpdata& sstats2) {
	double r = 0;
	for( int i=0; i < sstats1.w.size(); i++ ){
		for( int j=0; j < sstats2.w.size(); j++ ){
			double r_ij = 0.0;
			if( sstats1.info.chr[i] == sstats2.info.chr[j] ){
				r_ij = LDREF.gcorr( sstats1.info.chr[i], sstats1.info.ld_index[i] , sstats2.info.ld_index[j] );
			}else{
				// cerr << "\n WARNING: MISMATCHED CHR: " << sstats1.info.chr[i] << "(" << ts1 << ")" << " versus " << sstats2.info.chr[j] << "(" << ts2 << ")" << "\n";
			}
			r += (r_ij * sstats1.w[i] * sstats2.w[j] );
		}
	}
	
	r = r / sqrt( sstats1.w_adj * sstats2.w_adj );
	
	if( abs(r) > 1.0001 ){
		cout << "FATAL ERROR:: IMPOSSIBLE CORR VALUE = " << r << " \n\n";
		r = 0;
		for( int i=0; i < sstats1.w.size(); i++ ){
			for( int j=0; j < sstats2.w.size(); j++ ){
				double r_ij = 0.0;
				if( sstats1.info.chr[i] == sstats2.info.chr[j] ){
					r_ij = LDREF.gcorr( sstats1.info.chr[i], sstats1.info.ld_index[i] , sstats2.info.ld_index[j] );
				}
				cout << ( r_ij * sstats1.w[i] * sstats2.w[j] ) << "(" << r_ij << "," << sstats1.w[i] << "," <<  sstats2.w[j] << ")\t";
				r += (r_ij * sstats1.w[i] * sstats2.w[j] );
			}
			cout << "\n\n";
		}
		LDREF.print_LD(sstats1.info);
		cout << "\n";
		LDREF.print_LD(sstats2.info);
		
		cout << r << "\t" << sqrt(sstats1.w_adj)*sqrt(sstats2.w_adj) << "(" << sqrt(sstats1.w_adj) << "," << sqrt(sstats2.w_adj) << ")\n";
		
		cout << "\n\n";
		
		for( int i=0; i < sstats1.w.size(); i++ ){
				cout << sstats1.info.pos[i] << "\t";
		}
		
		cout << "\n";
		
		for( int i=0; i < sstats2.w.size(); i++ ){
				cout << sstats2.info.pos[i] << "\t";
		}
		
		cout << "\n\n";
		
		abort();
	}else{
		if( abs(r) > 1.00 ){
			if(r < 0){
				r = (-1.00);
			}else{
				r = 1.00;
			}
		}
	}
	
	return r;
}

Eigen::MatrixXd get_Lform_cov_mat(vector<string>& id1, vector<string>& id2, unitdata& dat1, unitdata& dat2) {
	Eigen::MatrixXd out( id1.size() , id2.size() );
	for (int i=0; i < id1.size(); ++i){
		for(int j=0; j < id2.size(); ++j){
			out(i, j) = Lform_cov( dat1.sstats[id1[i]], dat2.sstats[id2[j]] );
		}
	}
}

listLD get_Lform_cov_list(vector<string>& id1, vector<string>& id2, unitdata& dat1, unitdata& dat2) {
	listLD out;
	out.gene1 = dat1.gene_name;
	out.gene2 = dat2.gene_name;
	out.subclass1 = id1;
	out.subclass2 = id2;
	out.LDmat( id1.size() , id2.size() );
	for (int i=0; i < id1.size(); ++i){
		for(int j=0; j < id2.size(); ++j){
			out.LDmat(i, j) = Lform_cov( dat1.sstats[id1[i]], dat2.sstats[id2[j]] );
		}
	}
}

double unitdata::Lform_covar(string& ts1, string& ts2) {
	Lform_cov(sstats[ts1], sstats[ts2]);
}

Eigen::MatrixXd unitdata::get_Lform_covar() {
	Eigen::MatrixXd out( Lform_groups.size() , Lform_groups.size() );
	for (int i=0; i < Lform_groups.size(); ++i){
		out(i,i) = 1.00;
		for(int j=i+1; j < Lform_groups.size(); ++j){
			out(i, j) = Lform_covar( Lform_groups[i], Lform_groups[j] );
			out(j, i) = out(i, j);
			
		}
	}
	
	// cout << "\nLform corr matrix:\n" << out << "\n\n";
	// cout << "Lform z-scores:\n\n";
	
	// for (int i=0; i < Lform_groups.size(); ++i){
	// 	cout << sstats[Lform_groups[i]].qstat << "\t";
	// }
	// cout << "\n\n";
	
	pruneLD(out, 1);
	// cout << "\nPRUNED Lform corr matrix:\n" << out << "\n\n";
	
	return out;
}

Eigen::MatrixXd unitdata::get_Lform_covar(vector<double>& x) {
	Eigen::MatrixXd out( Lform_groups.size() , Lform_groups.size() );
	
	if( x.size() != Lform_groups.size() ){
		cerr << "FATAL : L-form Cov dimension does not match z-score vector size ... \n";
		abort();
	}
	
	for (int i=0; i < Lform_groups.size(); ++i){
		out(i,i) = 1.00;
		for(int j=i+1; j < Lform_groups.size(); ++j){
			out(i, j) = Lform_covar( Lform_groups[i], Lform_groups[j] );
			out(j, i) = out(i, j);
			
		}
	}
	
	// cout << "\nLform corr matrix:\n" << out << "\n\n";
	// cout << "Lform z-scores:\n\n";
	
	// for (int i=0; i < tissues.size(); ++i){
	// 	cout << sstats[tissues[i]].qstat << "\t";
	// }
	// cout << "\n\n";
	
	pruneLD(out, x, 1);
	// cout << "\nPRUNED Lform corr matrix:\n" << out << "\n\n";
	
	return out;
}

void unitdata::push(string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, double w, bool is_tissue, int& index) {
	if( tissues.size() == 0 & groups.size() == 0 ){
		chrom = chr; 
		start_pos = pos; 
		end_pos = pos; 
		omit_gene = false;
	}else{
		if( chr != chrom ){
			omit_gene = true;
		}
		start_pos = min(pos, start_pos); 
		end_pos = max(pos, end_pos); 
	}
	if( sstats.find(group) == sstats.end() ) {
		snpdata temp;
		sstats[group] = temp;
		if( is_tissue ) {
			sstats[group].anno_class = "eQTL";
			sstats[group].anno_subclass = group;
			tissues.push_back(group);
		}else {
			sstats[group].anno_class = ANNO_DEFS.meta[group];
			sstats[group].anno_subclass = group;
			groups.push_back(group);
		}
	}
	// cout << "unitdata::" << chr << ":" << pos << ":" << ref << ":" << alt << "\n";
	sstats[group].info.push(chr, pos, ref, alt, rsid, index);
	sstats[group].n.push_back(n);
	sstats[group].z.push_back(z);
	sstats[group].w.push_back(w);
}


bool unitdata::exists(string group) {
	if ( sstats.find(group) != sstats.end() ) {
		if(  sstats[group].z.size() > 0.0 ) {
			return true;
		}
	}
	return false;
}


void unitdata::update(dparam_l omega, vector<string> agroups) {
	double mmax = 0.0;
	for(string& group : agroups) {
		double incr = log(omega[group].p);
		if( exists( group ) ) {
			if( sstats[group].hsq < 0.0 ) {
				// sstats[group].bf = sstats[group].bf_partial - 0.5 * sstats[group].nvar * log( omega[group].tau );
				sstats[group].bf = sstats[group].bf_partial - 0.5 * log( omega[group].tau );
			}
			else {
				// sstats[group].bf_partial -= 0.5 * log(1 + sstats[group].ns*sstats[group].hsq * omega[group].tau);
				sstats[group].bf_partial -= 0.5 * log(1 + sstats[group].ns * omega[group].tau);
			}
			incr += sstats[group].bf;
			if( incr > mmax ) {
				mmax = incr;
			}
		}
		else {
			sstats[group].bf = 0.0;
		}
		par[group].p = incr;
	}
	double denom = exp(-mmax) * (1 - p_0) / p_0;
	for(string& group : agroups) {
		par[group].p = exp( par[group].p - mmax );
								 // || isnormal(par[*group].p) < 1 ){
		if( par[group].p < epsilon ) {
			par[group].p = epsilon;
		}
		denom += par[group].p;
	}
	for(string& group : agroups) {
		par[group].p /= denom;
		// cout << "\t" << *group << " BF = " << sstats[*group].bf;
		// cout << ", pp = " << par[*group].p << "\n";
	}
	// double total_p = 0.0;
	// for(vector<string>::iterator group = agroups.begin(); group != agroups.end(); ++group){
	// total_p += par[*group].p;
	// }
	// cerr << "\n\nOVERALL P::    " << total_p << "\n\n";
}

/*
void unitdata::skat(string group) {
	//if( sstats.find(group) == sstats.end() ) {
	//	return;
	//}else{
		sstats[group].skat();
		pvals[group] = sstats[group].pval;
	//}
}*/

void unitdata::runTest(string group, char test_type) {
	if( group == "TSS" ){
		sstats[group].anno_class = "TSS";
		sstats[group].cauchy_w(true);
	}else{
		switch( test_type ){
			case 'L' : sstats[group].burden(); break;
			case 'Q' : sstats[group].skat(); break;
			case 'C' : sstats[group].cauchy(); break;
			case 'T' : sstats[group].cauchy_w(true); sstats[group].anno_class = "TSS";
		}
	}
	pvals[group] = sstats[group].pval;
	for( int& pos : sstats[group].info.pos ){
		uniq_pos.insert(pos);
	}
} 


void gwasdata::runTest(string& name, string& group, bool is_element = false, char test_type = 'C'){
	
	data[name].runTest(group, test_type);
	
	if( is_element ){
		//string& enh_chrom = data[name].sstats[group].info.chr[0];
		
		//if( ENH_CAUCHY ){
			
		//}else{
		
			string gene = regel_data.gene_names[group][0];

			if( data.find(gene) == data.end() ) {
				unitdata temp;
				
				temp.gene_name = gene;
				temp.chrom = data[name].sstats[group].info.chr[0];
				temp.start_pos = data[name].sstats[group].info.pos[0];
				temp.end_pos = data[name].sstats[group].info.pos[0];
				
				data[gene] = temp;

				genes.push_back(gene);
			}

			//for( string &gene : regel_data.gene_names[group] ){
				//string& gene_chrom = data[gene].sstats[group].info.chr[0];
				//if( enh_chrom ==data[gene].chrom ){
					data[gene].regels.push_back(group);

					for( const int& ps : data[name].sstats[group].info.pos ){
						data[gene].start_pos = min(data[gene].start_pos, ps);
						data[gene].end_pos = max(data[gene].end_pos, ps);
						data[gene].uniq_pos.insert(ps);
					}

					data[gene].sstats[group] = data[name].sstats[group];
					data[gene].sstats[group].anno_class = name;
					data[gene].sstats[group].anno_subclass = regel_data.group_maps[group];
					data[gene].pvals[group] = data[name].pvals[group];
					data[gene].sstats[group].gene_list = regel_data.gene_maps[group];
				//}else{
				//	cerr << "WARNING: Mismatched chromosome: " << name << " " << group << "(chr " << enh_chrom << ") assigned to " << gene << "(chr " << data[gene].chrom << ")\n";
				//	cerr << "         " << name << "-gene assignment ignored. \n"; 
				//}
			//}
		//}
	}
}

string unitdata::print_summary(){
	string n_classes = "E=" + to_string(tissues.size()) + ",R=" + to_string(regels.size()) + ",O=" + to_string(groups.size());
	string out = chrom + "\t" + to_string(start_pos) + "-" + to_string(end_pos) + "\t" + gene_name + "\t" + to_string( uniq_pos.size() ) + "\t" + n_classes + "\t" + top_class + "\t" + top_subclass + "\t" + scien(min_pval,3)  + "\t" + scien(naive_pval,3) + "\t" + scien(global_pval,3);
	out += "\n";
	return out;
}

string unitdata::print_all_groups(){
	string out = "";
	for(string &x : groups){
		out += sstats[x].info.chr[0] + "\t" + to_string(sstats[x].info.min_pos) + "-" + to_string(sstats[x].info.max_pos) + "\t" + gene_name + "\t" + sstats[x].anno_class + "\t" + sstats[x].anno_subclass + "\t" +  to_string(sstats[x].nvar) + "\t" + to_string(sstats[x].qstat) + "\t" + scien(pvals[x])  + "\t" + sstats[x].flag + "\n";
	}
	for(string &x : tissues){
		out += sstats[x].info.chr[0] + "\t" + to_string(sstats[x].info.min_pos) + "-" + to_string(sstats[x].info.max_pos) + "\t" + gene_name + "\t" + sstats[x].anno_class + "\t" + sstats[x].anno_subclass + "\t" +  to_string(sstats[x].nvar) + "\t" + to_string(sstats[x].qstat) + "\t" + scien(pvals[x])  + "\t" + sstats[x].flag + "\n";
	}
	for(string &x : regels){
		if( PRINT_ALL_REGEL_GENES ){
			out += sstats[x].info.chr[0] + "\t" + to_string(sstats[x].info.min_pos) + "-" + to_string(sstats[x].info.max_pos) + "\t" + sstats[x].gene_list + "\t" + sstats[x].anno_class + "\t" + sstats[x].anno_subclass + "\t" +  to_string(sstats[x].nvar) + "\t" + to_string(sstats[x].qstat) + "\t" + scien(pvals[x])  + "\t" + sstats[x].flag + "\n";
		}else{
			out += sstats[x].info.chr[0] + "\t" + to_string(sstats[x].info.min_pos) + "-" + to_string(sstats[x].info.max_pos) + "\t" + gene_name + "\t" + sstats[x].anno_class + "\t" + sstats[x].anno_subclass + "\t" +  to_string(sstats[x].nvar) + "\t" + to_string(sstats[x].qstat) + "\t" + scien(pvals[x])  + "\t" + sstats[x].flag + "\n";
		}
	}
	return out;
}


// NEEDS_MOD
void snpdata::skat(){
	if( z.size() < 1 ) {
		return;
	}
	ns = 0.0;
	
	w_adj = 1;
	
	Eigen::VectorXd zv = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(z.data(), z.size());
	Eigen::VectorXd b = zv;
	nvar = n.size();

	//outstr << nvar << ":" << nvar-LDREF.nmiss(info) << "\t";

	double zzmax = 0.0;
	qstat = 0.0;
	for( int i = 0; i < nvar; i++) {
		ns += n[i];
		b(i) /= sqrt(n[i]);
		qstat += zv(i)*zv(i);
		if( zv(i)*zv(i) > zzmax ) {
			zzmax = zv(i)*zv(i);
		}
	}
	ns /= (double) nvar;
	Eigen::MatrixXd sigma = LDREF.LD(info);

	flag = "OK";
	
	if(DEBUG) {
		cout << "\nSNP info:\n";
		info.print();
		cout << "\nLD matrix: \n";
		LDREF.print_LD(info);
	}

	pval = liu_pval(z, sigma, logdet);
	
	Eigen::MatrixXd iden = Eigen::MatrixXd::Identity(nvar, nvar);

	double pval_no_ld = liu_pval(z, iden, logdet);

	if( pval_no_ld > pval && pval < 0.001 ) {
		flag = "LD_MISMATCH";
		if(DEBUG) cout << "\nLD flagged as problematic ... \n";
		pval = pval_no_ld;
	}

	/*logdet = log(sigma.determinant());
	Eigen::VectorXd bEi = sigma.colPivHouseholderQr().solve(b);
	qstat = zEi.transpose() * z;
	theta = bEi.transpose() * b;
	// qstat = ns *theta;
	*/

	theta = 0.0;

	double denom = 0.0;

	for( int i = 0; i < nvar; i++) {
		theta += b(i)*b(i)*exp(0.5*zv(i)*zv(i)-0.5*zzmax-0.5*log(n[i]) - log(nvar));
		denom += exp(0.5*zv(i)*zv(i)-0.5*zzmax-0.5*log(n[i])  - log(nvar));
	}
	theta /= denom;

	bf_partial = log(denom) + 0.5*zzmax;
	bf = log(denom) + 0.5*zzmax - 0.5*log(tau_0);

	nvar_total = nvar;
	hsq = -1;
	
	qstat = sqrt(qstat);
	
	// outstr << theta << "\t" << qstat  << "\t" << scien(pval) << "\t" << flag;
}

vector<double> expWeight(vector<double>& ds, double& alp){
	vector<double> out;
	for( const double& d_i : ds ){
		out.push_back(
			exp( -abs(d_i)*alp )
		);
	}
	return out;
}

string gsub(const string pattern, const string replacement, string str) {
    size_t start_pos = str.find(pattern);
    if(start_pos != string::npos){
		str.replace(start_pos, pattern.length(), replacement);
    }
	return str;
}

void snpdata::cauchy_w(bool weight_by_dtss){
	if( z.size() < 1 ) {
		return;
	}
	ns = 0.0;
	
	//double pval_thresh_c = 5e-8;

	w_adj = 1;
	qstat = 0.0;
	theta = 0.0;
	
	nvar = z.size();
	
	vector<double> pvals;
	//vector<double> wei;
	vector<double>& wei = w;
	
	/*if( PRUNE_CAUCHY ){
		map<double,double> z_w_map;
		for( int i = 0; i < nvar; i++) {
			double z_i = abs(z[i]);
			if( z_w_map.find(z_i) == z_w_map.end() ){ 
				z_w_map[z_i] = w[i];
			}else{
				z_w_map[z_i] = min(w[i], z_w_map[z_i]);
			}
		}
		for (const auto& z_w_k : z_w_map) {
			double z_k = z_w_k.first;
			double w_k = z_w_k.second;
			
			double zstat = z_k*z_k;
			
			if( zstat <= 0.0 ){
				zstat = simChiNZ();
			}
			
			double pval = pchisq(zstat, 1);
			
			pvals.push_back(pval);
			pval_thresh_c = min(pval_thresh_c, pval);
			wei.push_back( w_k );
		} 
		nvar = pvals.size();
	}else{*/

	for( int i = 0; i < nvar; i++) {
		double zstat = z[i]*z[i];
		if( zstat <= 0.0 ){
			zstat = simChiNZ();
		}
		pvals.push_back(
			pchisq(zstat, 1)
		);
		//pval_thresh_c = min(pval_thresh_c, pvals[i]);
	}

	
	flag = "OK";
	
	vector<double> pval_alpha;
	
	bool verbose = false;
	
	string alpha_flag = "OK|";
	
	if( weight_by_dtss ){
		
		for(double& alpha : TSS_ALPHA){
			vector<double> wei_alpha = expWeight(wei, alpha);
			double p_i = cauchy_test_weighted(pvals, wei_alpha);
			if( p_i <= TSS_VERBOSITY_PVAL ){
				verbose = true;
			}
			alpha_flag += gsub("e-0", "e-", scien(alpha,0)) + ":" + scien(p_i,2) + ",";
			pval_alpha.push_back(p_i);
		}
		
		pval = cauchy_minP(pval_alpha);
			
		if( pval <= TSS_VERBOSITY_PVAL ){
			verbose = true;
		}
		if(verbose){
			flag = alpha_flag;
		}
		
	}else{
		
		pval = cauchy_test_weighted(pvals, wei);
		
	}
	
	//cerr << pval << "\n";
	//abort();
	
	nvar_total = nvar;

}


void snpdata::cauchy(){
	if( z.size() < 1 ) {
		return;
	}
	ns = 0.0;
	
	w_adj = 1;
	qstat = 0.0;
	theta = 0.0;
	
	nvar = z.size();
	
	vector<double> pvals;

	for( int i = 0; i < nvar; i++) {
		double zstat = z[i]*z[i];
		if( zstat <= 0.0 ){
			zstat = simChiNZ();
		}
		pvals.push_back(
			pchisq(zstat, 1)
		);
	}

	flag = "OK";
	
	pval = cauchy_minP(pvals);
	
	nvar_total = nvar;

}


/*
void unitdata::burden(string group) {
	//if( sstats.find(group) == sstats.end() ) {
	//	return;
	//}else{
		sstats[group].burden();
	//}
}*/

void unitdata::multiBurden() {
	
	double z_max = 0;
	vector<double> z_stats;
	
	for(const string& ts : Lform_groups){
		z_stats.push_back( sstats[ts].qstat );
		if ( abs( sstats[ts].qstat ) > z_max ){
			z_max = abs( sstats[ts].qstat );
		}
	}
	
	Eigen::MatrixXd bcor = get_Lform_covar(z_stats);
	
	double sum_z_sq = 0.0;
	
	for( const double& zs : z_stats ){
		sum_z_sq += zs*zs;
	}
	
	double nv = bcor.rows();
	
	double pval_bonf = min(1.00, nv*pchisq( z_max*z_max, 1));
	double pval_mmvn = min(pval_bonf, MVN_minP(z_max, bcor));
	//double pval_skat = liu_pval(sum_z_sq, bcor);
	
	if( pval_mmvn <= 0 ){
		//cerr << "\nWARNING: MVN_minP <= 0; using Bonferroni p-val\n";
		pval_mmvn = pval_bonf;
	}else if( pval_mmvn > pval_bonf){
		//cerr << "\nWARNING: MVN_minP > Bonferroni p-val ... \n";
		pval_mmvn = pval_bonf;
	}
	
	/*if( pval_mmvn < 1e-4 ){
		cout << "\n\nMax Z-score = " << z_max << ", " << tissues.size() << " tissues, " << nv << " unique stats\n";
		cout << "\n\nBonferroni Pval = " << pval_bonf << "\n";
		cout << "\n\nSOCS Pval = " << pval_skat << "\n";
		cout << "Adjusted Minimum P-value = " << pval_mmvn << "\n\n";
	}*/
	//cout << "TEST_STATS\t" << z_max << "\t" << tissues.size() << "\t" << nv << "\t" << pval_bonf << "\t" << pval_skat << "\t" << pval_mmvn << "\n";
	
	Lform_MMVN_pval = pval_mmvn;
	
	return;
}


void unitdata::globalPval(){
	global_pval = 1.000;
	min_pval = 1.000;
	naive_pval = 1.000;
	n_uniq_tests = 0.00;
	vector<double> pval_vec;
	double n_total_tests = 0.00;
	top_subclass = "";
	top_class = "";
	if( tissues.size() > 0 ){
		for( const string& ts : tissues ){
			if( sstats[ts].pval <= min_pval ){
				if( sstats[ts].pval == min_pval  ){
					if( top_subclass != "" ){
						top_subclass += "," + ts;
					}
					if( top_class != sstats[ts].anno_class ){
						top_class += "," + sstats[ts].anno_class;
					}
				}else{
					top_subclass = ts;
					top_class = sstats[ts].anno_class;
				}
				min_pval = sstats[ts].pval;
			}
			n_total_tests++;
		}
	}
	Lform_groups = tissues;
	if( groups.size() > 0 ){
		for( const string& x : groups ){
			if( sstats[x].pval <= min_pval ){
				if( sstats[x].pval == min_pval  ){
					if( top_subclass != "" ){
						top_subclass += "," + x;
					}
					if( top_class != sstats[x].anno_class ){
						top_class += "," + sstats[x].anno_class;
					}
				}else{
					top_subclass = x;
					top_class = sstats[x].anno_class;
				}
				min_pval = sstats[x].pval;
			}
			if( sstats[x].z.size() > 1 ){
				pval_vec.push_back(sstats[x].pval);
				global_pval = min(global_pval, sstats[x].pval);
				n_uniq_tests++;
			}else if( sstats[x].z.size() == 1 ){
				Lform_groups.push_back(x);
			}
			n_total_tests++;
		}
	}
	if( regels.size() > 0 ){
		for( const string& x : regels ){
			if( sstats[x].info.chr[0] == chrom ){
				if( sstats[x].pval <= min_pval ){
					if( sstats[x].pval == min_pval  ){
						if( top_subclass != "" ){
							top_subclass += "," + sstats[x].anno_class;
						}else{
							top_subclass = sstats[x].anno_subclass;
						}
						if( top_class != sstats[x].anno_class ){
							top_class += "," + sstats[x].anno_class;
						}
					}else{
						top_subclass = sstats[x].anno_subclass;
						top_class = sstats[x].anno_class;
					}
					min_pval = sstats[x].pval;
				}
				if( sstats[x].z.size() > 1 ){
					pval_vec.push_back(sstats[x].pval);
					global_pval = min(global_pval, sstats[x].pval);
					n_uniq_tests++;
				}else if( sstats[x].z.size() == 1 ){
					Lform_groups.push_back(x);
				}
			}
			n_total_tests++;
		}
	}
	Lform_MMVN_pval = 1.00;
	if(Lform_groups.size() > 0 ){
		n_uniq_tests++;
		multiBurden();
		pval_vec.push_back(Lform_MMVN_pval);
	}
	if( n_uniq_tests < 1.00 ){
		n_uniq_tests = 1.00;
	}
	if( n_total_tests <= 1.00 ){
		global_pval = min_pval;
		naive_pval = min_pval;
		return;
	}
	/*for( const string& x : regels ){
		if( sstats[x].pval < min_pval ){
			min_pval = sstats[x].pval;
			top_subclass = x;
			top_class = "Enhancer";
		}
		n_uniq_tests++;
	}*/
	//global_pval = n_uniq_tests*min(Lform_MMVN_pval, global_pval);
	//global_pval = min(global_pval, 1.00);

	global_pval = cauchy_minP(pval_vec);

	naive_pval = min(n_total_tests*min_pval, 1.00);

	//if( global_pval < min_pval ){
	//	global_pval = min(min_pval*n_uniq_tests, 1.00);
	//}
	return;
}



// NEEDS_MOD
void snpdata::burden(){
	if( z.size() < 1.0 ) {
		return;
	}
	Eigen::MatrixXd sigma;

	if(DEBUG) {
		cout << "\nSNP info:\n";
		info.print();
		cout << "\nLD matrix: \n";
		LDREF.print_LD(info);
	}

	sigma = LDREF.LD(info);
	
	int m = z.size();
	
	// outstr << m << ":" << m-LDREF.nmiss(info) << "\t";

	if(DEBUG) {
		cout << "\nweights: \n";
		for( int i = 0; i < m; i++) {
			cout << "\t" << w[i];
		}
		cout << "\n\nz-stats: \n";
		for( int i = 0; i < m; i++) {
			cout << "\t" << z[i];
		}
		cout << "\n";
	}

	ns = 0.0;
	theta = 0.0;

	nvar = z.size();

	double denom = 0.0;
	double denom_no_ld = 0.0;
	double stat = 0.0;
	double offd = 0.0;

	flag = "OK";

	for( int i = 0; i < m; i++) {
		ns += n[i];
		stat += z[i] * w[i];
		denom += w[i] * w[i];
		denom_no_ld  += w[i] * w[i];
		theta += z[i] * w[i] / sqrt(n[i]);
		if(i > 0) {
			for( int j = 0; j < i; j++) {
				denom += 2*sigma(i,j)*w[i]*w[j];
			}
		}
	}
	
	w_adj = denom;

	if(DEBUG) cout << "\ntest score = " << stat << "\n\ntest variance = " << denom << "\n\n";

	if( stat*stat/denom > stat*stat/denom_no_ld && stat*stat/denom > 8 ) {
		flag = "LD_MISMATCH";
		denom = denom_no_ld;
		if(DEBUG) cout << "\nLD flagged as problematic ... \n";
	}
	ns /= (double) m;
	theta /= sqrt(denom);

	double zscore = stat/sqrt(denom);
	if( m == 1 ) {
		zscore = z[0];
		if( w[0] < 0 ) {
			zscore *= -1;
		}
	}

	qstat = zscore;
	bf_partial = 0.5 * zscore * zscore;
	bf = bf_partial - 0.5 * log(1 + ns*tau_0);

	pval = pchisq( zscore*zscore, 1);

	//outstr << theta << "\t" << qstat  << "\t" << scien(pval) << "\t" << flag;
}


void gwasdata::pushElement(string& name, string& elid, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, int& index) {
	if( data.find(name) == data.end() ) {
		unitdata temp;
		temp.gene_name = name;
		data[name] = temp;
		elements.push_back(name);
	}
	if( omega.find(name) == omega.end() ) {
		dparam temp;
		omega[name] = temp;
	}
	data[name].push(elid, chr, pos, ref, alt, rsid, n, z, 1, false, index);
}


void gwasdata::push(string& gene, string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, int& index) {
	if( data.find(gene) == data.end() ) {
		unitdata temp;
		temp.gene_name = gene;
		data[gene] = temp;
		genes.push_back(gene);
	}
	if( omega.find(group) == omega.end() ) {
		dparam temp;
		omega[group] = temp;
		groups.push_back(group);
	}
	data[gene].push(group, chr, pos, ref, alt, rsid, n, z, 1, false, index);
}


void gwasdata::push(string& gene, string& group, string& chr, int& pos, string& ref, string& alt, string& rsid, double& n, double& z, double w, int& index, bool is_tissue) {
	if( data.find(gene) == data.end() ) {
		unitdata temp;
		temp.gene_name = gene;
		data[gene] = temp;
		genes.push_back(gene);
	}
	if( omega.find(group) == omega.end() ) {
		dparam temp;
		omega[group] = temp;
		if(is_tissue){
			tissues.push_back(group);
		}else{
			groups.push_back(group);
		}
	}
	data[gene].push(group, chr, pos, ref, alt, rsid, n, z, w, is_tissue, index);
}


double gwasdata::EM_round(int iteration) {
	dparam_l omega_new = omega;
	double p_new = 0.0;
	double delta = 0.0;
	bool bprint = false;
	if( iteration % 10 == 0.0 ) {
		bprint = true;
	}
	else {
		cerr << ". . . ";
	}
	// cout << "ITERATION" << iteration << "...\n";
	for( string& gene : genes ) {
		// cout << gene << " ............... \n";
		data[gene].update( omega, groups );
	}
	// if( bprint ){
	// cerr << "\nstarting groups... \n";
	// }
	for(string& group : groups) {
		// if( bprint ){
		// cerr << "group " << *group << "...\n";
		// }
		omega_new[group].tau = 0.0;
		omega_new[group].p = 0.0;
		for(string& gene : genes) {
			if( data[gene].par[group].p > 0.0 && data[gene].par[group].tau > 0.0 && data[gene].sstats[group].bf != 0.0 ) {
				omega_new[group].p += data[gene].par[group].p;
				omega_new[group].tau += data[gene].par[group].p*data[gene].par[group].tau;
			}
			else {
				// cerr << "\nERROR IN " << gene << "\t" << group << "...\n";
				// cerr << "\tp = " << data[gene].par[group].p << "\ttau = " << data[gene].par[group].tau << "...\n\n";
			}
		}
		omega_new[group].tau /= omega_new[group].p;
		omega_new[group].p /= genes.size();
		if( bprint ) {
			cerr << "tau " << group << " :: " << omega_new[group].tau << "\n";
			cerr << "ppr " << group << " :: " << omega_new[group].p << "\n";
		}
								 // || 1.0 > isnormal(omega_new[group].p) ){
		if( omega_new[group].p < p_lower ) {
			// cerr << "\n\nERROR P : " << group << " :: " << omega_new[group].p << "  NEW: " << p_lower << "\n\n";
			omega_new[group].p = p_lower;
		}
		p_new += omega_new[group].p;
	}
	if( p_new > p_upper ) {
		// cerr << "\n\nTOTAL_ERROR P_TOTAL: " << p_new << "  NEW: " << p_upper << "\n\n";
		for(string& group : groups) {
			omega_new[group].p *= p_upper/p_new;
		}
		p_new = p_upper;
	}
	else if( p_new < p_lower ) {
		// cerr << "\n\nTOTAL_ERROR P_TOTAL: " << p_upper << " :: " << p_new << "  NEW: " << p_lower << "\n\n";
		for(string& group : groups) {
			omega_new[group].p *= p_lower/p_new;
		}
		p_new = p_lower;
	}
	if( bprint ) {
		cerr << "\n=========== EM ALGORITHM ===========\n";
		cerr << "ROUND " << iteration << " ............... \n\n";
	}
	for(string& group : groups) {
								 // || 1.0 > isnormal(omega_new[group].tau) ){
		if( omega_new[group].tau < tau_lower) {
			// cerr << "\n\nERROR TAU : " << group << " :: " << omega_new[group].tau << "  NEW: " << tau_lower << "\n\n";
			omega_new[group].tau = tau_lower;
		}
		omega_new[group].p /= p_new;
		delta += abs(omega_new[group].p - omega[group].p);
		delta += abs(omega_new[group].tau - omega[group].tau);
		if( bprint ) {
			cerr << group << ":\tconditional pr = " << omega_new[group].p << "\ttau = " << omega_new[group].tau << "\n";
		}
	}
	if(bprint) {
		cerr << "\t\t delta = " << delta << ",  pr causal = " << p_new << "\n\n";
	}
	p_0 = p_new;
	omega = omega_new;
	return delta;
}




