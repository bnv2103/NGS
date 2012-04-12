#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"

//from boost:
#include "boost/boost/math/distributions.hpp"
#include "boost/boost/math/distributions/poisson.hpp"
using namespace boost::math;

//from bamtools:
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamMultiReader.h>

//from fastahack:
#include "fastahack/Fasta.h"
#include "fastahack/disorder.h"
#include "fastahack/split.h"

using namespace std;
using namespace BamTools;       

class Tracking {
public:
	double prob;
	vector<int> v_path;
	double v_prob;
	Tracking();
	Tracking(double, vector<int> &, double);
};

void forward_viterbi(vector<int>&, vector<int>&, map<int, double>&, map<int, map<int, double> >&, vector<map<int, map<int, double> > >&);
void setVariables();
double retPowOverFact(int lam, int val);
int getLineCount(string);

#endif // MAIN_H
