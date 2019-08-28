#include "geneBasedTests.hpp"
#include "libMvtnorm/mvtnorm.h"
#include "cdflib/cdflib.hpp"
#include "distributionFunctions.hpp"

using namespace std;

double pval_HMP(vector<double>& p, vector<double>& w){
	
	// HMP references: Wilson (2019), PNAS
	
	double L = p.size();
	double HMP = 0;
	double sum_w = 0;
	
	for( int i = 0; i < L; i++){
		HMP += w[i]/p[i];
		sum_w += w[i];
	}
	HMP = sum_w/HMP;
	
	// Parameters from Wilson (2019):
	// LOC_L1 = 1 + digamma(1) + std::log(M_PI/2)
	// the location parameter when L==1, and
	// SCALE = M_PI/2 
	
	static double LOC_L1 = 0.874367040387922;
	static double SCALE = 1.5707963267949;
	
	return landau_ccdf(1/HMP, log(L) + LOC_L1, SCALE);
}

double pval_HMP(vector<double>& p){

	// HMP references: Wilson (2019), PNAS
	
	double L = p.size();
	double HMP = 0;

	for( int i = 0; i < L; i++){
			HMP += 1/p[i];
	}
	HMP = L/HMP;

	// Parameters from Wilson (2019):
	// LOC_L1 = 1 + digamma(1) + std::log(M_PI/2)
	// the location parameter when L==1, and
	// SCALE = M_PI/2

	static double LOC_L1 = 0.874367040387922;
	static double SCALE = 1.5707963267949;

	return landau_ccdf(1/HMP, log(L) + LOC_L1, SCALE);
}

double pval_ACAT(vector<double>& pvals){
	
	// ACAT references: Liu et al. (2019), AJHG; Liu and Xie (2018), arXiv
	
	if( pvals.size()==1 ){
		return pvals[0];
	}
	long double stat = 0.00;
	double pval_min = 1.00;
	for (double& p : pvals){
		stat += qcauchy( p );
		pval_min = min(pval_min, p); 
	}

	return pcauchy( stat/pvals.size() );
}

double pval_ACAT(vector<double>& pvals, vector<double>& weights){
	
	// ACAT references: Liu et al. (2019) AJHG; Liu and Xie (2018) arXiv
	
	if( pvals.size()==1 ){
		return pvals[0];
	}
	double stat = 0.0;
	double denom = 0.0;
	for (unsigned int i = 0; i < pvals.size(); i++){
		stat += weights[i]*qcauchy( pvals[i] );
		denom += weights[i];
	}
	return pcauchy( stat/denom );
}

double pval_LiuSKAT(vector<double>& Z, Eigen::MatrixXd& V, double& logdet) {
	
	// P-value approximation reference: Liu, Tang, Zhang (2009) Comput. Stat. Data Anal.
	// SKAT reference: Lee et al. (2012) AJHG
	
	double pval;
	int m = Z.size();
	logdet = 0.0;
	if( m != V.rows() ) {
			pval = -1;
	}
	else if( m < 1) {
			pval = 1;
	}
	else if( m == 1 ) {
			pval = pchisq(Z[0]*Z[0]/V.sum(), 1);
	}
	else {
		double stat = 0.0;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);
		es.compute(V);
		vector<double> SP(4, 0.0);
		double lambda;
		double lambda_p;
		for( int i = 0; i < m; i++) {
			stat += Z[i]*Z[i];
			lambda = es.eigenvalues()[i];
			lambda_p = lambda;
			logdet += log(lambda_p);
			for( int j = 0; j < 4; j++ ){
				SP[j] += lambda_p;
				lambda_p *= lambda;
			}
		}
		double s1 = SP[2]/pow(SP[1], 1.5);
		double s2 = SP[3]/(SP[1]*SP[1]);
		double a, d, l;
		if(s1*s1 > s2){
			a = 1/(s1 - sqrt(s1*s1 - s2));  
			d = (s1 *a - 1)*a*a;  
			l = a*a - 2*d;
		}else{
			l = 1/(s1*s1);  
			a = sqrt(l);  
			d = 0;
		}
		pval = pchisq( sqrt(2)*a*(stat - SP[0])/sqrt(2*SP[1]) + l + d , l, d);
	}
	return pval;
}

double pval_LiuSKAT(double sumZsq, Eigen::MatrixXd& V) {
	
	// P-value approximation reference: Liu, Tang, Zhang (2009) Comput. Stat. Data Anal.
	// SKAT reference: Lee et al. (2012) AJHG
	
	double pval;
	int m = V.cols();
	if( m != V.rows() ) {
			pval = -1;
	}
	else if( m < 1) {
			pval = 1;
	}
	else if( m == 1 ) {
			pval = pchisq(sumZsq/V.sum(), 1);
	}
	else {
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);
		es.compute(V);
		vector<double> SP(4, 0.0);
		double lambda;
		double lambda_p;
		for( int i = 0; i < m; i++) {
			lambda = es.eigenvalues()[i];
			lambda_p = lambda;
			for( int j = 0; j < 4; j++ ){
				SP[j] += lambda_p;
				lambda_p *= lambda;
			}
		}
		double s1 = SP[2]/pow(SP[1], 1.5);
		double s2 = SP[3]/(SP[1]*SP[1]);
		double a, d, l;
		if(s1*s1 > s2){
			a = 1/(s1 - sqrt(s1*s1 - s2));  
			d = (s1 *a - 1)*a*a;  
			l = a*a - 2*d; 
		}else{
			l = 1/(s1*s1);  
			a = sqrt(l);  
			d = 0;
		}
		pval = pchisq( sqrt(2)*a*(sumZsq - SP[0])/sqrt(2*SP[1]) + l + d , l, d);
	}
	return pval;
}

double pval_maxZsq(double& Z_max, Eigen::MatrixXd& Sigma ){
	
	int m = Sigma.rows();
	if( m != Sigma.cols() ){
		cerr << " << FATAL:: non-square LD matrix ... \n\n";
		abort();
	}
	int lt_len = 0.5*m*(m-1);
	
	double* Sigma_lt = new double[lt_len];
	
	double* upper = new double[m];
	double* lower = new double[m];
	double* delta = new double[m];
	int* infin = new int[m];
	
	int k = 0;
	for(int i = 0; i < m; i++ ){
		lower[i] = -abs(Z_max);
		upper[i] = abs(Z_max);
		delta[i] = 0.0;
		infin[i] = 2;
		for(int j = i+1; j < m; j++ ){
			Sigma_lt[k] = Sigma(i,j);
			k++;
		}
	}
	
	int nu_ = 0;
	int maxpts_ = 50000;
	double abseps_ = 1e-9;
	double releps_ = 0;
	
	double error;
	
	double value_ = 0.0;
	int inform_ = 0.0;
	
	mvtdst_(&m, &nu_, lower, upper, infin, Sigma_lt, delta, &maxpts_, &abseps_, &releps_, &error, &value_, &inform_);
	
	delete[] (Sigma_lt);
	delete[] (lower);
	delete[] (upper);
	delete[] (infin);
	delete[] (delta);
	
	return (1-value_);
}


double logBF( vector<double>& Z_scores,  Eigen::MatrixXd& Lambda, Eigen::MatrixXd& R, double& N, double& sigma2){
	
	int m = Z_scores.size(); 
	
	Eigen::VectorXd Z = Eigen::Map<Eigen::VectorXd>(Z_scores.data(), m);
	
	Eigen::MatrixXd R_Lambda = R * Lambda;
	
	Eigen::MatrixXd S = Eigen::MatrixXd::Identity(m,m) * sigma2/N + R_Lambda; 
	
	double qstat = (Z.transpose() * Lambda * S.ldlt().solve(Z) )(1,1)*0.5/sigma2;
	
	double logdet = log((R_Lambda*N/sigma2 + Eigen::MatrixXd::Identity(m,m)).partialPivLu().determinant());
	
	return (qstat - 0.5*logdet);
}

double logBF( double& Z,  double& Lambda, double& N, double& sigma2){ 
	
	double qstat = (Z * (Lambda/(sigma2/N + Lambda)  )* Z )*0.5/sigma2;
	
	double logdet = log(N/sigma2 + 1);
	
	return (qstat - 0.5*logdet);
}



