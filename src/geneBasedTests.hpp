#ifndef GENEBASEDTESTS_HPP
#define GENEBASEDTESTS_HPP

#include "distributionFunctions.hpp"

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
#include <set>
#include <map>

using namespace std;

double logBF( vector<double>&,  Eigen::MatrixXd&, Eigen::MatrixXd&, double&, double&);
double pval_HMP(vector<double>&);
double pval_HMP(vector<double>&, vector<double>&);

double pval_ACAT(vector<double>&);

double pval_ACAT(vector<double>&, vector<double>&);

double skat_pval(vector<double>& , Eigen::MatrixXd& , double& );

double pval_LiuSKAT(vector<double>& , Eigen::MatrixXd& , double& );
double pval_LiuSKAT(double , Eigen::MatrixXd& );

double pval_maxZsq(double& , Eigen::MatrixXd&  );

#endif
