#include "processGWAS.hpp"
#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <random>

using namespace std;

void print_usage() {
	cerr << "tools for analysis of GWAS summary statistics\n\n";
	cerr << "command line options \n";
	cerr << "    input\n";
	cerr << "      --gwas STR : annotated GWAS summary statistics file\n";
	cerr << "      --tss-bed STR : bed file specifying TSS locations \n";
	cerr << "      --anno-bed [data/anno.bed.gz] : annotation bed file\n";
	cerr << "      --anno-defs [defs.txt] : annotation hierarchy definitions\n";
	cerr << "      --betas [data/betas.txt.gz] : eSNP weight file\n";
	cerr << "      --ldref [data/chr*.vcf.gz] : LD reference panel (\"*\" as wildcard when split by chr)\n";
	cerr << "    variant filtering\n";
	cerr << "      --ldref-only : only retain variants with complete LD information\n";
	cerr << "      --tissues STR : only use eSNPs from tissues listed in file\n";
	cerr << "    output\n";
	cerr << "      --region STR : restrict analysis to specified region \n";
	cerr << "      --stdout : print to stdout rather than file \n";
	cerr << "      --prefix STR : prefix for output files \n";
	//cerr << "      --merge-tissues : use global rather than tissue-specific eQTL weights \n";
	cerr << "    other\n";
	cerr << "      --help : print this message and exit\n\n";
}


int main (int argc, char *argv[]) {
	cout.precision(5);

	string dfile = "";
	string afile = "";
	string sfile = "";
	string ldfile = "/net/snowwhite/home/corbinq/panels/G1K_EUR_3V5/chr$.vcf.gz";
	string outfile = "";
	string tfile = "";
	string tssfile = "";
	string bfile = "";
	string prefix = "";
	string target_gene = "";

	gwasdata gwinfo;
	annodef adefs;

	int region_mode = 0;
	string region = "";

	int extra = 0;
	int help = 0;
	int bayes = 0;
	int tmerge = 0;
	int debug_mode = 0;
	int exclude_missing = 0;
	int print_screen = 0;
	int twas = 0;
	
	char default_test = 'Q'; // default is 'Q' (SKAT)
	
	int cauchy_no_skat = 0;
	
	string TSS_VERBOSITY_PVAL_STR = "";
	double TSS_VERBOSITY_PVAL = 1.00;
	
	string TSS_ALPHA_STR = "";
	vector<double> TSS_ALPHA{1e-4, 5e-5, 1e-5, 5e-6};
	
	string TSS_WINDOW_STR = "";
	int TSS_WINDOW = 500000;

	int opt = 0;

	static struct option long_options[] = {
		{"help",   no_argument,  &help,  1},
		{"twas", no_argument, &twas, 1},
		{"ldref",      required_argument,  0,  'l'},
		{"gwas", required_argument, 0, 'g'},
		{"anno-defs", required_argument, 0, 'a'},
		{"defs", required_argument, 0, 'a'},
		{"anno-bed", required_argument, 0, 'f'},
		{"tss-bed", required_argument, 0, 's'},
		{"tss-alpha", required_argument, 0, 'x'},
		{"tss-window", required_argument, 0, 'y'},
		{"tss-verbosity", required_argument, 0, 'e'},
		{"tissues", required_argument, 0, 't'},
		{"betas", required_argument, 0, 'b'},
		{"merge-tissues", no_argument, &tmerge, 1},
		{"debug", no_argument, &debug_mode, 1},
		{"acat", no_argument, &cauchy_no_skat, 1},
		{"stdout", no_argument, &print_screen, 1},
		{"bayes", no_argument, &bayes, 1},
		{"region",    required_argument, 0,  'r' },
		{"prefix", required_argument, 0, 'p'},
		{"extra",    no_argument, &extra,  1},
		{"region",    required_argument, 0,  'r'},
		{"gene", required_argument, 0, 'v'},
		{"ldref-only", no_argument, &exclude_missing, 1},
		{0, 0, 0, 0}
	};
	int long_index =0;
	while ((opt = getopt_long(argc,argv,"p:l:g:a:t:s:b:d:",long_options,&long_index)) != -1) {
		switch (opt) {
			case 'f' : afile = optarg;
				break;
			case 'p' : prefix = optarg;
				break;
			case 'l' : ldfile = optarg;
				break;
			case 'd' : target_gene = optarg;
				break;
			case 'g' : sfile = optarg;
				break;
			case 'a' : dfile = optarg;
				break;
			case 't' : tfile = optarg;
				break;
			case 's' : tssfile = optarg;
				break;
			case 'x' : TSS_ALPHA_STR = optarg;
				break;
			case 'y' : TSS_WINDOW_STR = optarg;
				break;
			case 'e' : TSS_VERBOSITY_PVAL_STR = optarg;
				break;
			case 'b' : bfile = optarg;
				break;
			case 'r' : region = optarg; region_mode = 1;
				break;
			case '?': 
				/*if (isprint (optopt)){
					fprintf (stderr, "Unknown option `-%c'.\n", optopt); 
				}else{
					cerr << "Unknown option " << argv[optind - 1] << "\n";
				}*/
				cerr << "Try --help to see valid options.\n\n";
				return 1;
			//default:
			//	cerr << "Try --help to see options.\n\n";
			//	abort ();
		}
	}

	// if(fatal_opt) {
	// cerr << "hint: try --help to see all options \n\n";
	// return 1;
	// }
	if( print_screen ){
		cout << "## GAMBIT v0.2 (c) 2019 Corbin Quick (corbinq@gmail.com)\n";
	}else{
		cerr << "========================================================\n";
		cerr << " GAMBIT v0.2 (c) 2019 Corbin Quick (corbinq@gmail.com) \n";
		cerr << "========================================================\n\n";
	}
	
	if( help ) {
		print_usage();
		return 0;
	}
	if( ldfile == "") {
		cerr << "fatal error...  no LD reference files!\n\n";
		cerr << "hint: specify --ldref \"chr*.my_ref.vcf.gz\" ... \n";
		cerr << "      try --help to see more options\n\n";
		return 1;
	}
	else {
		setPath( ldfile );
	}
	if( debug_mode ) {
		debug();
	}
	if( sfile == "" ) {
		cerr << "fatal error...  no gwas input file!\n\n";
		cerr << "hint: specify --gwas my_summary_stats.txt.gz ... \n";
		cerr << "      try --help to see more options\n\n";
		return 1;
	}
	if( prefix == "" && !print_screen ) {
		//cerr << "FYI: no output --prefix specified!\n";
		//cerr << "     using \"output.funcwas\"  ... \n\n";
		prefix = sfile.substr(sfile.find("/")+1);
		prefix = prefix.substr(0, prefix.find(".gz"));
		prefix = prefix.substr(0, prefix.find(".txt"));
		prefix = prefix.substr(0, prefix.find(".az"));
		//prefix += ".output";
		cerr << "writing to " << prefix << ".*_out.txt ...\n\n";
	}
	if( exclude_missing ) {
		excludeMissing();
	}

	if( print_screen ){
		cerr.rdbuf(nullptr);
	}
	if( dfile == "" && afile == "" && bfile != "" ){
		twas = 1;
		cerr << "no definitions file or bed file specified...\n";
		cerr << "running in TWAS mode.\n\n";
	}

	initCHR();
	
	if( cauchy_no_skat ){
		default_test = 'C'; // rather than SKAT, using Cauchy/ACAT test (no LD)
	}
	
	if( TSS_VERBOSITY_PVAL_STR != "" ){
		TSS_VERBOSITY_PVAL = stod(TSS_VERBOSITY_PVAL_STR);
	}
	
	if( TSS_WINDOW_STR != "" ){
		TSS_WINDOW = stoi(TSS_WINDOW_STR);
	}
	
	if( TSS_ALPHA_STR != "" ){
		cerr << "specified TSS distance weight parameter = " << TSS_ALPHA_STR << " ... \n\n";
		TSS_ALPHA.clear();
		
		istringstream  iss(TSS_ALPHA_STR);
		string alpha_s;
		
		while(getline(iss, alpha_s, ',')) {
			TSS_ALPHA.push_back(stod(alpha_s));
		}
	}
	
	setVerbosityTSS(TSS_VERBOSITY_PVAL);
	setAlphaTSS(TSS_ALPHA);
	setWindowSizeTSS(TSS_WINDOW);

	if( tssfile != "" ) {
		cerr << "started processing TSS bed file...\n\n";
		gwinfo.tss_data.readTSS(tssfile, region);
	}

	if( bfile != "" ) {
		cerr << "started processing eSNP weights...\n\n";
		if( tmerge ) {
			cerr << "computing global eSNP weights only...\n\n";
		}
		if( tfile == "" ) {
			gwinfo.eweight_data.readBetas(bfile, region, tmerge);
		}else {
			gwinfo.eweight_data.readBetas(bfile, region, tfile, tmerge);
		}
	}

	if( afile != "" ) {
		cerr << "started processing annotation bed file...\n\n";
		gwinfo.regel_data.readElements(afile, region);
	}

	int n_haps;

	if( read_defs(dfile, adefs) ) {
		if( adefs.nodef ){
			//cerr << "no definitions file specified...\n\n";
		}else{
			cerr << "processed functional annotation definitions...\n\n";
		}
		int nstats = 0;
		int ngenes = 0;
		string fname = prefix + ".stratified_out.txt";
		string gname = prefix + ".summary_out.txt";
		string bname = prefix + ".bayes_out.txt";
		string hname = prefix + ".hyper_out.txt";
		if( read_gwas(sfile, region, gwinfo, adefs, twas) ) {
			if( target_gene != "" ) {
				ofstream genes_out(gname,ios::out);
				ofstream strat_out(fname,ios::out);
				
				genes_out << "#CHR\tPOS1\tPOS2\tGENE\tTOP_CLASS\tTOP_SUBCLASS\tMIN_UNADJ_PVAL\tPVAL\n";
				strat_out << "#CHR\tPOS\tGENE\tCLASS\tSUBCLASS\tNSNPS\tSTAT\tPVAL\tINFO\n";
				
				for(string& group : gwinfo[target_gene].groups) {
					if( exclude_missing ) gwinfo.data[target_gene][group].pop_missing();
					gwinfo[target_gene].runTest(group,default_test);
					nstats++;
				}
				for(string& tissue : gwinfo[target_gene].tissues) {
					if( exclude_missing ) gwinfo.data[target_gene][tissue].pop_missing();
					gwinfo[target_gene].runTest(tissue,'L');
				}
				gwinfo[target_gene].globalPval();
				strat_out << gwinfo[target_gene].print_all_groups();
				genes_out << gwinfo[target_gene].print_summary();
				
				strat_out.close();
				genes_out.close();
			}
			else {
				
				std::streambuf * buf;
				std::ofstream of;
				
				if(!print_screen){
					of.open(fname,ios::out);
					buf = of.rdbuf();
				}else{
					buf = cout.rdbuf();
				}
				
				ofstream genes_out(gname,ios::out);
				ostream strat_out(buf);
				
				strat_out << "#CHR\tPOS\tGENE\tCLASS\tSUBCLASS\tNSNPS\tSTAT\tPVAL\tINFO\n";
				
				genes_out << "#CHR\tPOS\tGENE\tN_SNPS\tN_CLASSES\tTOP_CLASS\tTOP_SUBCLASS\tMIN_UNADJ_PVAL\tNAIVE_PVAL\tPVAL\n";
				
				if ( gwinfo.elements.size() > 0 ){
					cerr << "starting regulatory elements ... \n";
					int nelem = 0;
					for(string& name : gwinfo.elements) {
						for(string& group : gwinfo[name].groups) {
							if( exclude_missing ) gwinfo.data[name][group].pop_missing();
							gwinfo.runTest(name,group,true,default_test);
							nelem++;
							if( nelem % 25 == 0 ) {
								cerr << "\rcalculated statistics for " << pretty(nelem) << " regulatory elements  ";
							}
						}
						//strat_out << gwinfo[name].print_all_groups();
					}
					cerr << "\rcalculated statistics for " << pretty(nelem) << " regulatory elements! \n\n";
				}
				
				
				for(string& gene : gwinfo.genes ) {
					// cout << "\nBEGIN GENE " << *gene << "\n";
					for(string& group : gwinfo[gene].groups) {
						if( exclude_missing ) gwinfo.data[gene][group].pop_missing();
						gwinfo.data[gene].runTest(group,default_test);
						nstats++;
						if( nstats % 25 == 0 ) {
							cerr << "\rcalculated " << pretty(nstats) << " statistics for " << pretty(ngenes) << " genes  ";
						}
					}
					for( string& tissue : gwinfo[gene].tissues ) {
						if( exclude_missing ) gwinfo.data[gene][tissue].pop_missing();
						gwinfo.data[gene].runTest(tissue,'L');
						nstats++;
						if( nstats % 25 == 0 ) {
							cerr << "\rcalculated " << pretty(nstats) << " statistics for " << pretty(ngenes) << " genes  ";
						}
					}
					gwinfo[gene].globalPval();
					strat_out << gwinfo[gene].print_all_groups();
					genes_out << gwinfo[gene].print_summary();
					ngenes++;
					// cout << "\nEND GENE " << *gene << "\n";
				}
				cerr << "\rcalculated " << pretty(nstats) << " statistics for " << pretty(ngenes) << " genes!  \n\n";
				

				//strat_out.close();
				genes_out.close();

				if(bayes) {
					ofstream bayesout(bname,ios::out);
					/*bayesout << "#chr\tpos\tgene\tlogBF\ttop_group\tpp_top_group\tpp_total\ttheta_top\n";
					double del = 500;
					int rounds = 0;
					while( ( rounds < 2000 && del > 0.0001 ) || rounds < 50 ) {
						del = gwinfo.EM_round(rounds);
						rounds++;
					}
					cerr << "\n";
					// cout << "## bayes stuff...\n" ;
					for(vector<string>::iterator gene = gwinfo.genes.begin(); gene != gwinfo.genes.end(); ++gene) {
						string wmax;
						double vmax = -100000000000;
						double mbf = -100000000000;
						for(vector<string>::iterator group = gwinfo.groups.begin(); group != gwinfo.groups.end(); ++group) {
							// cout << "\t" << *group << " BF = " << gwinfo.data[*gene][*group].bf;
							// cout << ", pp = " << gwinfo.data[*gene].par[*group].p << "\n";
							if( gwinfo.data[*gene].exists(*group) ) {
								// cout << *gene << "\t" << *group << "\t" << gwinfo.data[*gene].par[*group].p << "\t" << gwinfo.data[*gene][*group].bf << "\t" << gwinfo.data[*gene][*group].bf_partial << "\n" ;
								if( gwinfo.data[*gene][*group].bf > mbf ) {
									vmax = gwinfo.data[*gene].par[*group].p;
									wmax = *group;
									mbf = gwinfo.data[*gene][*group].bf;
								}
							}
						}
						gwinfo[*gene].print_range(wmax, bayesout);
						bayesout << "\t" << *gene << "\t";
						double tbf = 0;
						double pp = 0;
						for(vector<string>::iterator group = gwinfo.groups.begin(); group != gwinfo.groups.end(); ++group) {
							tbf += gwinfo.omega[*group].p*exp(gwinfo.data[*gene][*group].bf - mbf);
							pp += gwinfo.data[*gene].par[*group].p;
						}
						tbf = log(tbf) + mbf;
						bayesout << tbf << "\t" << wmax << "\t" << vmax << "\t" << pp << "\t" << gwinfo.data[*gene].par[wmax].tau << "\n";
					}
					bayesout.close();
					cerr << "\nbayes analysis complete ... ";
					cerr << "\nwriting hyperparameter estimates to " << hname << " ... \n";
					ofstream hyperout(hname,ios::out);
					hyperout << "GROUP\tPI\tTAU\tP_GENES\tN_GENES\tPI_FOLD\n";
					for(vector<string>::iterator group = gwinfo.groups.begin(); group != gwinfo.groups.end(); ++group) {
						hyperout << *group << "\t" << gwinfo.omega[*group].p << "\t" << gwinfo.omega[*group].tau << "\t";
						double n_genes = 0;
						for(vector<string>::iterator gene = gwinfo.genes.begin(); gene != gwinfo.genes.end(); ++gene) {
							if( gwinfo.data[*gene].exists(*group) ) {
								n_genes++;
							}
						}
						hyperout << n_genes/gwinfo.genes.size() << "\t" << n_genes << "\t" << gwinfo.omega[*group].p/(n_genes/gwinfo.genes.size()) << "\n";
					}*/
					
					bayesout << "#CHR\tSTART\tEND\tID\tGENE\tCLASS\tSUBCLASS\tBF_PARTIAL\tMLE\tN\tNVAR\n";
					for( string& gene : gwinfo.genes ) {
						for( string group : gwinfo[gene].groups ) {
							bayesout << gwinfo[gene].print_range(group) << "\t"
							 << gene << ":" << group << "\t"
							 << gene << "\t"
							 << adefs.meta[group] << "\t"
							 << group << "\t"
							 << gwinfo.data[gene][group].bf_partial << "\t"
							 << gwinfo.data[gene][group].theta << "\t"
							 << gwinfo.data[gene][group].ns << "\t"
							 << gwinfo.data[gene][group].nvar << "\n";
						}
						for(string& tissue : gwinfo[gene].tissues) {
							bayesout << gwinfo[gene].print_range(tissue) << "\t"
							 << gene << ":" << tissue << "\t"
							 << gene << "\t"
							 << "eQTL" << "\t"
							 << tissue << "\t"
							 << gwinfo.data[gene][tissue].bf_partial << "\t"
							 << gwinfo.data[gene][tissue].theta << "\t"
							 << gwinfo.data[gene][tissue].ns << "\t"
							 << gwinfo.data[gene][tissue].nvar << "\n";
						}
					}
					int neid = 0;
					for(string& name : gwinfo.elements) {
						for(string& group : gwinfo[name].groups) {
							bayesout << gwinfo[name].print_range(group) << "\t"
							 << "ELEM" << neid << "\t"
							 << gwinfo.regel_data.gene_maps[group] << "\t"
							 << name << "\t"
							 << gwinfo.regel_data.group_maps[group] << "\t"				
							 << gwinfo.data[name][group].bf_partial << "\t"
							 << gwinfo.data[name][group].theta << "\t"
							 << gwinfo.data[name][group].ns << "\t"
							 << gwinfo.data[name][group].nvar << "\n";
							 neid++;
						}
					}
					bayesout.close();
				}
			}
			cerr << "\ndone! thanks for using GAMBIT\n\n";
			return 0;
		}
		else {
			cerr << "fatal error... invalid input " << sfile << "!\n";
			return 1;
		}
	}
	else {
		cerr << "fatal error... invalid definitions " << dfile << "!\n";
		return 1;
	}
//	closeAll();
	return 0;
}
