#include "processLD.hpp"
#include "tabixpp/tabix.hpp"
#include "common.hpp"

using namespace std;

TabixList tfiles;
string chr_global = "0";
string vcf_path = "";

int JUMP_DIST = 200000;

int ld_panel_size = 0;

string chr_pfx = "";
bool checked_prefix = false;

bool memoize_LD = true;
bool PRELOAD_PANEL = false;

unordered_map<string, int> iCHR;
vector<string> sCHR;

ldref LDREF;

void setMemoizeLD(bool memLD){
	memoize_LD = memLD;
}

void setPreload(bool pl){
	PRELOAD_PANEL = pl;
}

void setJump(int jd){
	JUMP_DIST = jd;
}

void initCHR() {
	sCHR.resize(26);
	sCHR[0] = "0";
	for( int i = 1; i < 23; i++) {
		iCHR[to_string(i)] = i;
		iCHR["chr"+to_string(i)] = i;
		sCHR[i] = to_string(i);
	}
	iCHR["X"] = 23;
	iCHR["chrX"] = 23;
	sCHR[23] = "X";
	iCHR["Y"] = 24;
	iCHR["chrY"] = 24;
	sCHR[23] = "X";
	iCHR["M"] = 25;
	iCHR["chrM"] = 25;
	sCHR[25] = "M";
	for( int i = 1; i < 26; i++) {
		tfiles.opened[i-1] = false;
		ld_data tmp(i);
		LDREF.chr_data.push_back(tmp);
	}
}


void setPath(string& path) {
	vcf_path = path;
}



string gsubstr(string inp, string pat, string rep) {
	std::string::size_type n = 0;
	while ( ( n = inp.find( pat, n ) ) != std::string::npos ) {
		inp.replace( n, pat.size(), rep );
		n += rep.size();
	}
	return inp;
}


void snpinfo::push(string& ch, int& po, string& rf, string& al, string& rs) {
	chr.push_back(ch);
	pos.push_back(po);
	rsid.push_back(rs);
	ref.push_back(rf);
	alt.push_back(al);
	
	if( chr.size() == 1 ){
		min_pos = po;
		max_pos = po;
	}else{
		min_pos = min(po, min_pos);
		max_pos = max(po, max_pos);
	}
}


void snpinfo::push(string& ch, int& po, string& rf, string& al, string& rs, int& index) {
	chr.push_back(ch);
	pos.push_back(po);
	rsid.push_back(rs);
	ref.push_back(rf);
	alt.push_back(al);

	if( chr.size() == 1 ){
		min_pos = po;
		max_pos = po;
	}else{
		min_pos = min(po, min_pos);
		max_pos = max(po, max_pos);
	}
	LDREF.chr_data[iCHR[ch]-1].fetch(po, rf, al, index);

	ld_index.push_back(index);
}


void snpinfo::print() {
	for( int i = 0; i < pos.size(); i++) {
		cout <<  chr[i] << "\t";
		cout <<  pos[i] << "\t";
		cout << rsid[i] << "\t";
		cout <<  ref[i] << "\t";
		cout <<  alt[i] << "\n";
	}
}


int snpinfo::size() {
	return pos.size();
}


void snpinfo::pop(vector<int>& rm){
	sort(rm.begin(), rm.end(), greater<int>());
	for( const int& i : rm){
		chr.erase(chr.begin() + i);
		pos.erase(pos.begin() + i);
		rsid.erase(rsid.begin() + i);
		ref.erase(ref.begin() + i);
		alt.erase(alt.begin() + i);
		ld_index.erase(ld_index.begin() + i);
	}
	min_pos = 300000000;
	max_pos = 0;
	for( const int& p : pos ){
		min_pos = min(min_pos, p);
		max_pos = max(max_pos, p);
	}
}


ld_datum::ld_datum() {
	dir = 0;
	ld_count = 0;
	missing = true;
	fetched = false;
}


ld_datum::ld_datum(string& ch, int& po, string& re, string& al) {
	dir = 0;
	ld_count = 1;
	missing = true;
	fetched = false;
	chr = ch;
	pos = po;
	ref = re;
	alt = al;
}


ld_data::ld_data(int i) {
	ichr = i;
	schr = sCHR[i];
	nhaps = 0;
	nsnps = 0;
}


void ld_data::fetch(int& pos, string& ref, string& alt, int& index) {
	if( index == 0 ) {
		string key = to_string(pos) + ":" + ref + ":" + alt;
		if( snp_index.find(key) == snp_index.end() ) {
			ld_datum tmp(schr, pos, ref, alt);
			if( PRELOAD_PANEL ){
				tmp.fetchGeno();
			}
			pos_data.push_back(tmp);
			
			snp_index[key] = nsnps;
			index = nsnps;
			
			nsnps++;
		}else{
			index = snp_index[key];
			pos_data[index].ld_count++;
		}
	}
}

double ld_data::gcorr(int i, int j) {
	if(i == j) {
		return 1;
	}else if(i > j){
		swap(i,j);
	}
	ld_datum& pd_i = pos_data[i];
	ld_datum& pd_j = pos_data[j];
	if( memoize_LD ){
		if( pd_i.ld_count > 1 && pd_j.ld_count > 1 ){
			
			if( ld_vals.find(pair<int,int>(i,j)) != ld_vals.end() ){
				
				double& Rv = ld_vals[pair<int,int>(i,j)];
				
				if( abs(Rv) > 1 ){
					cout << "LD val = " << Rv << "\n\n";
					cout << " i = " << i << " , j = " << j << "\n\n";
					abort();
				}
				return Rv;
				
			}
		}
	}
	if( !pd_i.fetched ){
		pd_i.fetchGeno();
	}
	if( !pd_j.fetched ){
		pd_j.fetchGeno();
	}
	if(pd_i.dir == 0.0 || pd_j.dir == 0.0 ) {
		return 0;
	}
	if( pd_i.nhaps != pd_j.nhaps ){
		cerr << "WARNING: " << schr << ":" << pd_i.pos << " has n=" << pd_i.nhaps << ", while " << schr << ":" << pd_j.pos << " has n=" << pd_j.nhaps << "...\n";
		return 0;
	}else if( nhaps != pd_i.nhaps){
		nhaps = pd_i.nhaps;
	}
	if(min(pd_i.mac,pd_j.mac)==0.0 || max(pd_i.mac,pd_j.mac)==nhaps) {
		return 0;
	}

	double EG_i = (double) pd_i.mac/ (double) nhaps;
	double EG_j = (double) pd_j.mac/ (double) nhaps;
	double EGij;

	if( EG_i < EG_j ) {
		EGij = sizeIS(pd_i.carriers, pd_i.mac, pd_j.genotypes, nhaps);
	}
	else {
		EGij = sizeIS(pd_j.carriers, pd_j.mac, pd_i.genotypes, nhaps);
	}

	double SD_i = sqrt(EG_i)*sqrt(1.0 - EG_i);
	double SD_j = sqrt(EG_j)*sqrt(1.0 - EG_j);
	double C_ij = EGij - EG_i*EG_j;

	double D = pd_i.dir * pd_j.dir * C_ij;
	double R;

	if( SD_i * SD_j > 0.0 ) {
		R  = D / ( SD_i * SD_j );
	}
	if( abs(R) > 1.0 ) {
		R = (R > 0) ? 1.0 : -1.0;
	}
	
	if( memoize_LD ){
		if( pd_i.ld_count > 1 && pd_j.ld_count > 1 ){
			ld_vals[pair<int,int>(i,j)] = R;
		}
	}
	
	return R;
}

double ldref::gcorr(string& ch, int& i, int& j) {
	return chr_data[ iCHR[ch]-1 ].gcorr(i,j);
}

void ld_data::print_LD(vector<int>& idx) {
	for( const int& i : idx ) {
		for( const int& j : idx ) {
			if(i == j) {
				cout << "1.000" << "\t";
			}
			else {
				cout << gcorr(i, j) << "\t";
			}
		}
		cout << "\n";
	}
}


Eigen::MatrixXd ld_data::LD( vector<int>& idx ) {
	int M = idx.size();
	Eigen::MatrixXd sigma(M, M);
	for( int i = 0; i < M-1; i++) {
		sigma(i,i) = 1;
		for( int j = i + 1; j < M; j++) {
			sigma(i, j) = gcorr(idx[i], idx[j]);
			sigma(j, i) = sigma(i, j);
		}
	}
	sigma(M-1, M-1) = 1;
	return sigma;
}


int ld_data::nmiss() {
	int out = 0;
	for(int i = 0; i < nsnps; i++) {
		if( !pos_data[i].fetched ){
			pos_data[i].fetchGeno();
		}
		if( pos_data[i].missing ) {
			out++;
		}
	}
	return out;
}


int ld_data::nmiss( vector<int>& idx ) {
	int out = 0;
	if( idx.size() <= 1 ){
		return out;
	}
	for(vector<int>::iterator i = idx.begin(); i != idx.end(); ++i) {
		if( !pos_data[*i].fetched ){
			pos_data[*i].fetchGeno();
		}
		if( pos_data[*i].missing ) {
			out++;
		}
	}
	return out;
}


vector<int> ld_data::missing( vector<int>& idx ) {
	vector<int> out;
	if( idx.size() <= 1){
		return out;
	}
	for(int i = idx.size()-1; i >= 0; i--){
		if( !pos_data[idx[i]].fetched ){
			pos_data[idx[i]].fetchGeno();
		}
		if( pos_data[idx[i]].missing ) {
			out.push_back(i);
		}
	}
	return out;
}


void ldref::fetch(string& chr, int& pos, string& ref, string& alt, int& index) {
	if( index == 0 ) {
		chr_data[iCHR[chr]-1].fetch(pos, ref, alt, index);
	}
}

void ldref::print_LD(string& chr, vector<int>& idx ) {
	chr_data[iCHR[chr]-1].print_LD(idx);
}

void ldref::print_LD(snpinfo& sd) {
	print_LD( sd.chr[0], sd.ld_index );
}

Eigen::MatrixXd ldref::LD(string& chr, vector<int>& idx ) {
	return chr_data[iCHR[chr]-1].LD(idx);
}


Eigen::MatrixXd ldref::LD(snpinfo& sd) {
	return LD( sd.chr[0], sd.ld_index );
}


int ldref::nmiss(string& chr, vector<int>& idx ) {
	return chr_data[iCHR[chr]-1].nmiss(idx);
}


int ldref::nmiss(snpinfo& sd) {
	return nmiss( sd.chr[0], sd.ld_index );
}


vector<int> ldref::missing(string& chr, vector<int>& idx ){
	return chr_data[iCHR[chr]-1].missing(idx);
}


vector<int> ldref::missing(snpinfo& sd) {
	return missing(sd.chr[0], sd.ld_index );
}


void hdata::push_block (vector<int> ii, vector<int> hc, int nh) {
	iids.push_back(ii);
	hcts.push_back(hc);
	nhap.push_back(nh);
}


void hdata::push_map (vector<int> ma) {
	map.push_back(ma);
}


double sizeIS(vector<int> &iid1, int &mac1, vector<bool> &genotypes2, int &n ) {
	double out = 0.0;
	for( int k = 0; k < mac1; k++) {
		if( genotypes2[ iid1[k] ] ) {
			out++;
		}
	}
	
	return out/ ((double) n);
}


vector<int> getRegion (string str) {
	vector<int> v ;
	size_t prev_pos = 0, pos;
	while ((pos = str.find_first_of(":-,", prev_pos)) != std::string::npos) {
		if (pos > prev_pos) {
			v.push_back(stoi(str.substr(prev_pos, pos-prev_pos)));
		}
		prev_pos= pos+1;
	}
	if (prev_pos< str.length()) {
		v.push_back(stoi(str.substr(prev_pos, std::string::npos)));
	}
	return v;
}


string asRegion (string chr, int pos, int end) {
	string out = chr_pfx + gsubstr(chr, "chr", "") + ":" + to_string(pos) + "-" + to_string(end);
	return out;
}


string asRegion (int chr, int pos, int end) {
	string out = chr_pfx + to_string(chr) + ":" + to_string(pos) + "-" + to_string(end);
	return out;
}

string asRegion (string chr, int pos) {
	string out = chr_pfx + gsubstr(chr, "chr", "") + ":" + to_string(pos);
	return out;
}


string asRegion (int chr, int pos) {
	string out = chr_pfx + to_string(chr) + ":" + to_string(pos);
	return out;
}


string aflip (string &in) {
	if( in == "A" ) {
		return "T";
	}
	else if( in == "C" ) {
		return "G";
	}
	else if( in == "T" ) {
		return "A";
	}
	else if( in == "G" ) {
		return "C";
	}
	else {
		return "";
	}
}


int allele_check( string &ref, string &alt, string &ref0, string &alt0 ) {
	if( ref == ref0 && alt == alt0 ) {
		return 1;
	}
	else if( ref == alt0 && alt == ref0 ) {
		return -1;
	}
	else if( ref == aflip(ref0) && alt == aflip(alt0) ) {
		return 1;
	}
	else if( ref == aflip(alt0) && alt == aflip(ref0) ) {
		return -1;
	}
	else {
		return 0;
	}
}


void checkup() {
	for(int i = 0; i < 23; i++) {
		int jj = LDREF.chr_data[i].pos_data.size();
		if( jj > 0 ) cout << i+1 << "\t" << LDREF.chr_data[i].pos_data.size() << "\t" <<  LDREF.chr_data[i].nmiss() << "\n";
	}
}


string subchr (string inp, string chr) {
	chr = gsubstr(chr, "chr", "");
	if( inp.find("*") != string::npos ) {
		return gsubstr ( inp,  "*", chr);
	}else if( inp.find("$chr") != string::npos ) {
		return gsubstr ( inp,  "$chr", chr);
	}else if( inp.find("$") != string::npos ) {
		return gsubstr ( inp,  "$", chr);
	}
	return inp;
}

void addGeno(char& x, vector<bool>& geno, vector<int>& carr, int& n){
	switch (x){
		case '1' : 
			geno.push_back(1);
			carr.push_back(n);
			n++;
			break;
		case '0' : 
			geno.push_back(0);
			n++;
			break;
		default : 
			cerr << "FATAL ERROR: Complete, phased genotypes required in LD reference panel \n\n";
			abort();
	}
}

void setGeno(char& x, vector<bool>& geno, vector<int>& carr, int& n){
	switch (x){
		case '1' : 
			geno[n] = 1;
			carr.push_back(n);
			n++;
			break;
		case '0' :
			n++;
			break;
		default : 
			cerr << "FATAL ERROR: Complete, phased genotypes required in LD reference panel \n\n";
			abort();
	}
}

void ld_datum::fetchGeno() {
	if( fetched ){
		return;
	}
	int cdx = iCHR[chr]-1;
	if( !tfiles.opened[cdx] ) {
		//chr_global = chr;
		string chr_path = subchr(vcf_path, chr);
		struct stat buffer;
		if( stat (chr_path.c_str(), &buffer) == 0 ) {
			tfiles.ts[cdx].open( chr_path );
			tfiles.opened[cdx] = true;
			
			if( !checked_prefix ){
				string line; 
				while ( tfiles.ts[cdx].getNextLine(line) && !checked_prefix ){
					if( line[0] != '#' && line.length() > 0 ) {
							if( line.substr(0,3) == "chr" ){
								chr_pfx = "chr";
							}
							checked_prefix = true;
					}
				}
			}
			
		}
		else {
			cerr << "\nFATAL ERROR: " << vcf_path << " does not exist ... !\n";
			abort();
			// missing = true;
			// fetched = true;
			// return;
		}
	}
	
	bool passed = false;
	bool jumped = false;
	bool contin = true;
	
	missing = true;
	fetched = true;
	string region = asRegion( chr, pos );
	
	// cout << region << "\n";
	string line;
	while ( tfiles.ts[cdx].getNextLine(line) && contin ) {
		if( line[0] != '#' && line.length() > 0 ) {
			string ichr, iref, ialt, rsid, qual, filter, info, format, gts;
			int ipos;
			istringstream iss(line);
			iss >> ichr >> ipos >> rsid >> iref >> ialt >> qual >> filter >> info >> format;
			
			if( format.substr(0,2) != "GT" ){
				cerr << "\nFATAL ERROR: GT field is required for all variants in LD reference files\n\n";
				abort();
			}
			// cout << ichr << ":" << ipos << "\n";
			if( ipos == pos && iCHR[ichr] == iCHR[chr] ) {
				dir = allele_check(ref, alt, iref, ialt);
				if( abs(dir) > 0 ) {
					missing = false;
					contin = false;

					// if(DEBUG) cout << "SUCCESS\n";
					
					mac = 0;
					nhaps = 0;
					
					if( ld_panel_size == 0 ){
						
						while( iss >> gts ) {
							if( gts[1] != '|' ){
								cerr << "FATAL ERROR: phased genotypes required in LD reference panel\n\n";
								abort();
							}
							addGeno(gts[0], genotypes, carriers, nhaps);
							addGeno(gts[2], genotypes, carriers, nhaps);
						}
						mac = carriers.size();
						
						ld_panel_size = nhaps;
					}else{
						
						genotypes.resize(ld_panel_size, 0);
						while( iss >> gts ) {
							if( gts[1] != '|' ){
								cerr << "FATAL ERROR: phased genotypes required in LD reference panel\n\n";
								abort();
							}
							setGeno(gts[0], genotypes, carriers, nhaps);
							setGeno(gts[2], genotypes, carriers, nhaps);
						}
						if( nhaps != ld_panel_size ){
								cerr << "FATAL ERROR: inconsistent sample size in LD reference data\n\n";
								abort();
						}
						
						mac = carriers.size();
						
					}
					
					
					if( mac > 0.5*nhaps ){
						mac = nhaps - mac;
						dir *= (-1);
						genotypes.flip();
						carriers.clear();
						int i = 0;
						for( const bool gt : genotypes ){
							if(gt){
								carriers.push_back(i);
							}
							i++;
						}
					}
					// cout << chr << ":" << pos << ", N = " << nhaps << ", MAC = " << mac << "\n";
					// for( int j: carriers ){
					//	cout << j << ",";
					// }
					// cout << "\n";

				}
				//else if( DEBUG ) {
					// cout << "FAILED:" << ref << ":" << alt << "\n";
				//}
			}else{
				if( ipos < pos ){
					if( !jumped && pos - ipos >= JUMP_DIST ){
						tfiles.ts[cdx].setRegion( region );
						jumped = true;
					}
				}else{
					if( !passed ){
						passed = true;
						if( !jumped ){
							tfiles.ts[cdx].setRegion( region );
							jumped = true;
						}else{
							contin = false;
						}
					}else{
						contin = false;
					}
				}
				
			}
			//else if(DEBUG) {
				// cout << "FAILED:" << pos << "\n";
			//}
		}
	}
}
