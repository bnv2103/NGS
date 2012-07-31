#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

using namespace std;

#include "read.h"

READ::READ(void)
{
	snp_count = -1;
	hap = 0;
	hprob = 1.0;
	snps = new SNP*[1000];
	alleles = new char[1000];
cout << "Invoking plain READ" << endl;
}

READ::~READ()
{
//	cout << "Read destructor called: " << start << endl;
	delete [] alleles;
	delete [] snps;
}

READ::READ(const READ &read)
{
	*this = read;
}

READ::READ (int i)
{
	cout << "This is an int initialization of READ. This should not be invoked" << endl;
}

READ::READ (double i)
{
	cout << "This is a double initialization of READ. This should not be invoked" << endl;
}

READ::READ(long st, int len)
{
	length = len;
	snp_count = -1;
	hap = 1;
	hprob = 0.5;
	start = st;
	snps = new SNP*[1000];
	alleles = new char[1000];
}

void READ::addsnp(SNP *snp, char allele)
{
	snps[++snp_count] = snp;
	alleles[snp_count] = allele;
}

void READ::assignHaplotype(int haplotype, double prob)
{
	if(haplotype!=1&&haplotype!=2) {
		cout << "Assigning incorrect haplotype value " << haplotype << endl;
	}
//cout << "Prob:: " << prob << endl;
	if(prob<0.0||prob>1.0) {
		cout << "Assigning incorrect haplotype probability " << prob << " to haplotype " << haplotype << endl;
	}
	hap = haplotype;
	hprob = prob;
}

int READ::GetHap(void)
{
	return hap;
}

double READ::GetHapProb(void)
{
	return hprob;
}

long READ::GetPos(void)
{
	return start;
}

int READ::GetLen(void)
{
	return length;
}

SNP* READ::GetSnp(int pos)
{
	return snps[pos];
}

SNP** READ::GetSnpList()
{
	return snps;
}

int READ::GetSnpCount(void)
{
	return snp_count + 1;
}

char READ::GetAllele(int pos)
{
	return alleles[pos];
}

READ* READ::operator=(int i)
{
	cout << "Assigning int to READ. Should not be invoked" << endl;
	return this;
}

READ* READ::operator+=(double i)
{
	cout << "Adding double to READ. Should not be invoked" << endl;
	return this;
}

READ* READ::operator-=(double i)
{
	cout << "Subtracting double from READ. Should not be invoked" << endl;
	return this;
}

READ* READ::operator*=(double i)
{
	cout << "Multiplying READ by double. Should not be invoked" << endl;
	return this;
}

READ* READ::operator/=(double i)
{
	cout << "Dividing READ by double. Should not be invoked" << endl;
	return this;
}

READ* READ::operator+=(READ* read)
{
	cout << "Adding two READs. Implement" << endl;
	return this;
}

READ* READ::operator-=(READ* read)
{
	cout << "Subtracting two reads. Implement" << endl;
	return this;
}

READ* READ::operator*=(READ* read)
{
	cout << "Multiplying two reads. Implement??" << endl;
	return this;
}

READ* READ::operator/=(READ* read)
{
	cout << "Dividing two reads. Implement??" << endl;
	return this;
}

READ::operator char()
{
	cout << "Typecasting READ to char. Implement??" << endl;
	return 'q';
}

READ::operator int()
{
	cout << "Typecasting READ to int. Implement??" << endl;
	return 0;
}

READ::operator double()
{
	cout << "Typecasting READ to double. Implement??" << endl;
	return 0.0;
}

double READ::operator*(READ* read)
{
	cout << "Returning double from two multiplied reads. Should not be invoked" << endl;
	return 0.0;
}

double READ::operator-(double i)
{
	cout << "Subtracting double from READ to return double. Should not be invoked" << endl;
	return 0.0;
}

double READ::operator+(double i)
{
	cout << "Adding double to read to return double. Should not be invoked" << endl;
	return 0.0;
}

