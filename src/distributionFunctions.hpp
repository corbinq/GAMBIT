#ifndef DISTRIBUTIONFUNCTIONS_HPP
#define DISTRIBUTIONFUNCTIONS_HPP

#include "cdflib/cdflib.hpp"
#include "eigenmvn/eigenmvn.hpp"
#include "ROOT_Math/Landau.hpp"

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

double pcauchy(double);
double qcauchy(double);

double pchisq(double, double);
double pchisq(double, double, double);
double qchisq(double, double);

#endif
