#include "distributionFunctions.hpp"

using namespace std;

double pchisq( double stat, double df) {
	double cdf, ccdf;
	cumchi ( &stat, &df, &cdf, &ccdf );
	return ccdf;
}

double pchisq( double stat, double df, double ncp) {
	double cdf, ccdf;
	cumchn( &stat, &df, &ncp, &cdf, &ccdf);
	return ccdf;
}

double qchisq( double cdf, double df) {
	double stat, ccdf, bd;
	int which = 2;
	ccdf = 1 - cdf;
	int status;
	cdfchi ( &which, &cdf, &ccdf, &stat, &df, &status, &bd );
	return stat;
}

double pcauchy(double x){
	return 0.5 + atan( x )/M_PI;
}

double qcauchy(double q){
	return tan( M_PI*(q - 0.5) );
}
