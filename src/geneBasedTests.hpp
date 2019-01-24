#ifndef GENEBASEDTESTS_HPP
#define GENEBASEDTESTS_HPP

#include "eigenmvn.hpp"

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

using namespace std;

double logBF( vector<double>&,  Eigen::MatrixXd&, Eigen::MatrixXd&, double&, double&);

// double pcauchy(double);
// double qcauchy(double);

double cauchy_minP(vector<double>&);

double cauchy_test_weighted(vector<double>&, vector<double>&);

double pchisq(double, double);
double pchisq(double, double, double);
double qchisq(double, double);

double skat_pval(vector<double>& , Eigen::MatrixXd& , double& );

double liu_pval(vector<double>& , Eigen::MatrixXd& , double& );
double liu_pval(double , Eigen::MatrixXd& );

double MVN_minP(double& , Eigen::MatrixXd&  );

#endif
