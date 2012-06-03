#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

using namespace std;

#include "snp.h"

SNP::SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string> gl3, bool known_par)
{
	refcount = altcount = nrefcount = 0;
	count = -1;
	position = snp_pos;
	ref = snp_ref;
	alt = snp_alt;
	reads = new READ*[50];
	add_read(type, read);
	gl = new double[3];
	gl[0] = pow(10, -(atof(gl3[0].c_str()))/10);
	gl[1] = pow(10, -(atof(gl3[1].c_str()))/10);
	gl[2] = pow(10, -(atof(gl3[2].c_str()))/10);
	known = known_par;
}

long SNP::GetPos()
{
	return position;
}

void SNP::append(int type, READ *read)
{
	add_read(type, read);
}

void SNP::add_read(int type, READ *read)
{
	if(type==0)
		refcount++;
	else if(type==1)
		altcount++;
	else
		nrefcount++;
	reads[++count] = read;
}

READ* SNP::GetRead(int pos)
{
	return reads[pos];
}

char SNP::GetRef()
{
	return ref;
}

char SNP::GetAlt()
{
	return alt;
}

int SNP::GetReadCount()
{
	return count + 1;
}

int SNP::GetRefCount()
{
	return refcount;
}

double* SNP::GetGenLik()
{
	return gl;
}

bool SNP::GetKnown()
{
	return known;
}

