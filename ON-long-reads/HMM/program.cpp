#include<stdio.h>
#include <assert.h>
#include <math.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<string>
#include <iomanip>
#include <fstream>
#include<map>
#include<vector>
#include<sstream>

using namespace std;

#include "utils.h"
#include "read.h"
#include "snp.h"
#include "obs.h"
#include "obsSeq.h"
#include "obsProb.h"
#include "discreteObsProb.h"
#include "gaussianObsProb.h"
#include "vectorObsProb.h"
#include "flexibleObsProb.h"
#include "stateTrans.h"
#include "plainStateTrans.h"
#include "gammaProb.h"
#include "explicitDurationTrans.h"
#include "initStateProb.h"
#include "hmm.h"

#define ncells 1

map<long,string> snp_map;

vector<SNP*> snp_list;
vector<READ*> reads_list;
vector<string> known_snps;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim))
		elems.push_back(item);
	return elems;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

char seq_type(const char* seq, int *slen)
{
	int i = 0;
	int len = 0;

	while(seq[i]>='0'&&seq[i]<='9') {
		len = 10*len + seq[i]-'0';
		i++;
	}
	if(len==0) {
		cout << "Invalid cigar. Beginning with " << seq[0];
		exit(1);
	}
	*slen = len;
	return seq[i];
}

void read_known_snp_file(const char *file)
{
	string line;
	FILE *snp_file = fopen(file, "r");
	if(snp_file==NULL) {
		printf("Cannot open snp file %s\n",file);
		exit(1);
	}

	while(true) {
		char str[20];
		fgets(str,sizeof(str),snp_file);
		if(feof(snp_file))
			break;
		line = str;
		vector<string> snp_line = split(line, ' ');
		if(snp_line[1]!="N") {
			known_snps.insert(known_snps.end(),snp_line[0]);
		}
	}
}

void read_snp_file(const char *snpfile)
{
	string line;
	FILE *snp_file = fopen(snpfile, "r");
	if(snp_file==NULL) {
		printf("Cannot open snp file %s\n",snpfile);
		exit(1);
	}

	while(true) {
		char str[1000];
		fgets(str,sizeof(str),snp_file);
		if(feof(snp_file))
			break;
		if(str[0]=='#')
			continue;
		line = str;
		vector<string> vcf_line = split(line, '\t');
		if(vcf_line[4].find(',')==string::npos) {
			snp_map.insert(pair<long,string>(atol(vcf_line[1].c_str()),vcf_line[3]+vcf_line[4]+vcf_line[9]));
		}
	}
	fclose(snp_file);
cout << "SNP MAP size = " << snp_map.size() << endl;
}

void read_sam_files(char *fq_files)
{
	int i;
	string line;
	FILE *fq_file = fopen(fq_files, "r");
	if(fq_file==NULL) {
		printf("Cannot open fq file %s\n",fq_files);
		exit(1);
	}

	while(true) {
		char str[20000];
		fgets(str,sizeof(str),fq_file);
		if(feof(fq_file))
			break;
		char start = str[0];

		if(start != '@') { // for each new read in the file
			int slen;
			line = str;
			vector<string> read_line = split(line, '\t');
			long start = atol(read_line[3].c_str());
			if(seq_type(read_line[5].c_str(),&slen)=='S') {
				start -= slen;
			}
			int length = (int)(read_line[9].length());
			int nsnps = 0;
			READ *read = new READ(start,length); // create new read object

			map<long,string>::iterator low = snp_map.lower_bound(start);
			map<long,string>::iterator high = snp_map.upper_bound(start+length-1);
	   			if( (low == snp_map.begin() || low->first <= start+length-1) && (high == snp_map.end() || high->first >= start) ) {
				for(map<long,string>::iterator snp_it1=low; snp_it1!=high; snp_it1++) {
					int snp_exists = 0;
					int type = 2;
					bool known = false;
					char snp_ref = snp_it1->second.c_str()[0];
					char alt_ref = snp_it1->second.c_str()[1];
					char allele = read_line[9][snp_it1->first-start];
					for(vector<string>::iterator known_snpit = known_snps.begin(); known_snpit != known_snps.end(); known_snpit++) {
						if(atol((*known_snpit).c_str())>snp_it1->first) {
							break;
						}
						if(atol((*known_snpit).c_str())==snp_it1->first) {
							known = true;
							break;
						}
					}
					string info = snp_it1->second;
					vector<string> format = split(info, ':');
					vector<string> gl3 = split(format[1], ',');

					if(allele==snp_ref)
						type = 0;
					else if(allele==alt_ref)
						type = 1;
					for(vector<SNP*>::iterator snp_it2 = snp_list.begin(); snp_it2 != snp_list.end(); snp_it2++) { // check for snps in vector
						if((*snp_it2)->GetPos() == snp_it1->first) {
							(*snp_it2)->append(type,read);
							snp_exists = 1;
							read->addsnp(*snp_it2, allele);
							break;
						}
					}
					if(!snp_exists) {
						SNP *snp = new SNP(snp_it1->first,snp_it1->second.c_str()[0],snp_it1->second.c_str()[1],type,read,gl3,known);
						snp_list.insert(snp_list.end(), snp);
						read->addsnp(snp, allele);
					}
				}
			}
			vector<READ*>::iterator rend = reads_list.end();
			reads_list.insert(rend, read); // add it to the vector
		}
	}
	fclose(fq_file);

	cout << "number of reads = " << reads_list.size() << endl;
	cout << "number of snps = " << snp_list.size() << endl;
}

void runHMM()
{
	int i;
	int nbDimensions = 1, nbSymbols = 250, nbStates = 2, T = 200;
	int isGaussian = 0;
	int *listNbSymbols = new int[(nbDimensions)+1];;
	double **transitionMatrix, **emissionMatrix;
	double transProb = 0.5, emissionProb;

	for (i=1;i<=nbDimensions;i++){
		listNbSymbols[i] = nbSymbols;
	}

	const char *outputName = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output";
	char hmmOutputName[256];
	char stateOutputName[256];
	char expectedObsOutputName[256];
	char distanceOutputName[256];

	sprintf(hmmOutputName,"%s%s",outputName, ".hmm");
	sprintf(stateOutputName,"%s%s",outputName, ".sta");
	sprintf(expectedObsOutputName,"%s%s",outputName, ".obs");
	sprintf(distanceOutputName,"%s%s",outputName, ".obs");

	CStateTrans *a;
	CObsProb *b;
	CInitStateProb *pi;
	CHMM *learnedHMM;
	CObs *obsType;

	obsType = new CFlexibleObs<READ>(nbDimensions);

	transitionMatrix = SetMatrix(nbStates, nbStates);
	AssignMatrix(transitionMatrix, transProb, nbStates, nbStates);
	a = new CPlainStateTrans(nbStates, transitionMatrix);

	emissionMatrix = SetMatrix(nbStates, nbSymbols);
	b = new CFlexibleObsProb(listNbSymbols,nbStates,nbDimensions,isGaussian);
	b->InitStateProb();

	pi = new CInitStateProb(nbStates);

	learnedHMM = new CHMM(a, b, pi);
 
	// 2. Read observation sequence into data structure
 
	CObsSeq *obsSeq = new CObsSeq(obsType, &snp_list, &reads_list);

	ofstream stateOutput(stateOutputName);
	ofstream expectedObsOutput(expectedObsOutputName);
	ofstream distanceOutput(distanceOutputName);
	stateOutput.close();
	expectedObsOutput.close();
	double normalizedLogProb = learnedHMM->FindViterbiDistance(obsSeq, distanceOutput, &reads_list);
	distanceOutput.close();

	delete learnedHMM;
	delete pi;
	delete b;
	delete a;
	//delete obsSeq;// check with purify
}

int main(int argc, char **argv)
{
	int i;
	const char *ext = ".sam";
	const char *snpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/simulated_reads_4.vcf";
	const char *file_base = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/test_4";
	const char *knownsnpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/common_snps_21_short.txt";

	char file_name[100];
	sprintf(file_name, "%s%s", file_base, ext);

	read_snp_file(snpfile);
	read_known_snp_file(knownsnpfile);
	read_sam_files(file_name);

	runHMM();

	return 0;
}

