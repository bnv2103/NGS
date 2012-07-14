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
#define FULL
//#define DEBUG

int region;
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

vector<string> split(const string &s, char delim,int flag) {
if(flag) cout << s << endl;
	vector<string> elems;
	return split(s, delim, elems);
}

int rev_s_len(const char *seq)
{
	int i = strlen(seq)-1;
	int len = 0;
	int j = 1;

	if(seq[i] != 'S')
		return 0;
	else {
		i--;
		while(seq[i]>='0'&&seq[i]<='9') {
			len = 10*j*(seq[i]-'0') + len;
			i--;
			j *= 10;
		}
		if(len==0) {
			cout << "Invalid CIGAR. Ends with " << seq[strlen(seq)-1] << endl;
			exit(1);
		}
		return len/10;
	}
}

int s_len(const char* seq)
{
	int i = 0;
	int len = 0;

	while(seq[i]>='0'&&seq[i]<='9') {
		len = 10*len + seq[i]-'0';
		i++;
	}
	if(len==0) {
		cout << "Invalid cigar. Beginning with " << seq[0] << endl;
		exit(1);
	}
	if(seq[i]=='S')
		return len;
	else
		return 0;
}

int get_deletion_count(string cigar, int limit)
{
	char *c_cigar = new char[cigar.length()+1];
	c_cigar = (char*)cigar.c_str();
	c_cigar[cigar.length()] = '\0';
	int it = 0;
	int len = 0;
	int dlen = 0;
	int mlen = 0;
	int ilen = 0;
	int slen;

	if(limit==-1) limit = 10000;
	while(it<strlen(c_cigar)) {
		if(c_cigar[it]>='0'&&c_cigar[it]<='9') {
			len = len*10 + c_cigar[it] - '0';
		} else if(c_cigar[it]=='S') {
			slen += len;
			len = 0;
		} else if (c_cigar[it]=='M') {
			mlen += len;
			len = 0;
/*
 * REVISIT: For future implementation, we need to account for deleted snps themselves
	if(mlen+dlen==limit) {
		cout << "Mlen = limit" << endl;
	}
*/
			if(mlen+dlen==limit) {
// Experimental:
//				return -1;
				break;
			} else if(mlen+dlen > limit) {
				break;
			}
		} else if(c_cigar[it]=='D') {
			dlen += len;
			len = 0;
		} else if(c_cigar[it]=='I') {
			ilen += len;
			len = 0;
		} else {
			cout << "Invalid character in cigar string " << c_cigar[it] << endl;
		}
		it++;
	}
	return dlen - ilen;
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
		vector<string> snp_line = split(line, ' ',0);
		if(snp_line[1]!="N") {
			known_snps.insert(known_snps.end(),snp_line[0]);
		}
	}
}

void read_snp_file(const char *snpfile)
{
	string line;
	int hash_ctr = 0;
	int known_hash = 10;
	int homon_hash = 100;
	FILE *snp_file = fopen(snpfile, "r");
	if(snp_file==NULL) {
		printf("Cannot open snp file %s\n",snpfile);
		exit(1);
	}

	vector<string>::iterator vec_start = known_snps.begin();
	while(true) {
		char str[1000];
		fgets(str,sizeof(str),snp_file);
		if(feof(snp_file))
			break;
		if(str[0]=='#')
			continue;
		line = str;
		vector<string> vcf_line = split(line, '\t',0);
		if(vcf_line[4].find(',')==string::npos) {
			long pos = atol(vcf_line[1].c_str());
			string info = vcf_line[9];
			vector<string> format = split(info, ':',0);
			vector<string> gl3 = split(format[1], ',',0);
			int known = 0;
			for(vector<string>::iterator known_snpit = vec_start; known_snpit != known_snps.end(); known_snpit++) {
				if(atol((*known_snpit).c_str()) > pos) {
					break;
				}
				if(atol((*known_snpit).c_str()) == pos) {
					hash_ctr = distance(known_snps.begin(), known_snpit);
					if((hash_ctr%known_hash)+1!=known_hash) {
						known = 1;
					} else if((hash_ctr%homon_hash)+1!=homon_hash) {
						known = 2;
					} else {
						known = 3;
					}
					vec_start = known_snpit;
					break;
				}
			}
			SNP *snp = new SNP(pos, vcf_line[3].c_str()[0], vcf_line[4].c_str()[0], gl3, known, atof(vcf_line[5].c_str()));
			snp_list.insert(snp_list.end(), snp);
//cout << "Insert.." << sizeof(*snp) << endl;
		}
	}
	fclose(snp_file);
	cout << "SNP MAP size = " << snp_list.size() << endl;
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

	vector<SNP*>::iterator snp_begin = snp_list.begin();
	while(true) {
		char str[20000];
		fgets(str,sizeof(str),fq_file);
		if(feof(fq_file))
			break;
		char start = str[0];

		if(start != '@') { // for each new read in the file
			line = str;
			vector<string> read_line = split(line, '\t',0);
			long start = atol(read_line[3].c_str());
#ifdef DEBUG
cout << start << endl;
#endif
			int slen = s_len(read_line[5].c_str());
			int rev_slen = rev_s_len(read_line[5].c_str());
			string cigar = read_line[5];
			int Ndels = get_deletion_count(cigar,-1);
if(Ndels==-1) {
	cout <<"Unexpected -1 for get_deletion_count" << endl;
}
			int length = (int)(read_line[9].length()) - slen - rev_slen + Ndels;
#ifdef DEBUG
cout << slen << "\t" << rev_slen << "\t" << Ndels << "\t" << length << endl;
#endif
			READ *read = new READ(start,length); // create new read object

//cout << start << "\t" << (*snp_begin)->GetPos() << endl;
			for(vector<SNP*>::iterator snp_it = snp_begin; snp_it != snp_list.end(); snp_it++) {
				if((*snp_it)->GetPos() >= start) {
					if((*snp_it)->GetPos() <= start+length-1) {
						int snp_exists = 0;
						int type = 2;
						char snp_ref = (*snp_it)->GetRef();
						char alt_ref = (*snp_it)->GetAlt();
						int ndels = get_deletion_count(cigar,(*snp_it)->GetPos()-start);
//if(ndels==-1) break;
						char allele = read_line[9][(*snp_it)->GetPos() - start + slen - ndels];

						if(allele==snp_ref)
							type = 0;
						else if(allele==alt_ref)
							type = 1;

						read->addsnp(*snp_it, allele);
						(*snp_it)->append(type,read);
					} else {
//cout << start+length-1 << "\t" << (*snp_begin)->GetPos() << "\t" << "missing/breaking.." << endl;
						break;
					}
				} else {
					snp_begin = snp_it;
//cout << "catch up" << "\t" << (*snp_begin)->GetPos() << endl;
				}
			}
//cout << "READ" << endl;
			vector<READ*>::iterator rend = reads_list.end();
			reads_list.insert(rend, read); // add it to the vector
		}
	}
	fclose(fq_file);

	cout << "number of reads = " << reads_list.size() << endl;
	cout << "number of snps = " << snp_list.size() << endl;
}

void runHMM(const char *outputName)
{
	int i;
//	int nbDimensions = 1, nbSymbols = 140000, nbStates = 2;
	int nbDimensions = 1, nbSymbols = reads_list.size(), nbStates = 2;
	int isGaussian = 0;
	int *listNbSymbols = new int[(nbDimensions)+1];;
	double **transitionMatrix, **emissionMatrix;
	double transProb = 0.5, emissionProb;

	for (i=1;i<=nbDimensions;i++){
		listNbSymbols[i] = nbSymbols;
	}

//	const char *outputName = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/";
	char stateOutputName[256];
	char expectedObsOutputName[256];
	char distanceOutputName[256];

	sprintf(stateOutputName,"%s%s",outputName, ".sta");
	sprintf(distanceOutputName,"%s%s",outputName, ".obs");

	CStateTrans *a;
	CObsProb *b;
	CInitStateProb *pi;
	CHMM *learnedHMM;
	CObs *obsType;

	obsType = new CFlexibleObs<READ*>(nbDimensions);

	transitionMatrix = SetMatrix(nbStates, nbStates);
	AssignMatrix(transitionMatrix, transProb, nbStates, nbStates);
	a = new CPlainStateTrans(nbStates, transitionMatrix);

	emissionMatrix = SetMatrix(nbStates, nbSymbols);
	b = new CFlexibleObsProb(listNbSymbols,nbStates,nbDimensions,isGaussian);
	b->InitStateProb();

	pi = new CInitStateProb(nbStates);

	learnedHMM = new CHMM(a, b, pi);
 
	// 2. Read observation sequence into data structure
	// GLOBAL
	CObsSeq *obsSeq = new CObsSeq(obsType);
	//CObsSeq *obsSeq = new CObsSeq(obsType, &snp_list, &reads_list);

	ofstream stateOutput(stateOutputName);
	ofstream distanceOutput(distanceOutputName);
	//GLOBAL
//	double normalizedLogProb = learnedHMM->FindViterbiDistance(obsSeq, distanceOutput, &reads_list, &snp_list);
	//double normalizedLogProb = learnedHMM->FindViterbiDistance(obsSeq, distanceOutput,stateOutput);
	//
	// Supposed to run forward and backward algorithms
	// Obtain P(x) and thus posteriors
	// Perform haplotype calling
	// Perform genotype calling
	// Perform EM
	learnedHMM->FindFBDistance(obsSeq, distanceOutput);
	distanceOutput.close();

	cout << endl << endl;
	for(vector<SNP*>::iterator snp_it = snp_list.begin(); snp_it != snp_list.end(); snp_it++) { // check for snps in vector
		double *gt = new double[3];
		gt = (*snp_it)->GetPosteriors();
		int gp = gt[0]>gt[1] ? (gt[0]>gt[2] ? 0:2) : (gt[1]>gt[2] ? 1:2);
		stateOutput << (*snp_it)->GetPos() << "\t" << (*snp_it)->GetKnown() << "\t" << gp << "\t" << gt[0] << "\t" << gt[1] << "\t" << gt[2] << endl;
	}
	stateOutput.close();
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

#ifdef FULL
	const char *snpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/vcfs/0.01_0.04_2_5000_2000_0_20.1.1.vcf";
	const char *knownsnps="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_20.list";
	const char *file_base = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reads/0.01_0.04_2_5000_2000_0_20.sorted";
	char knownsnpfile[200], *base;

	if(argc>1) {
		file_base=argv[1];
		snpfile=argv[2];
		region=atoi(argv[3]);
		sprintf(knownsnpfile,"/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_%d.list",region);
		base=argv[4];
	} else {
		strcpy(knownsnpfile,knownsnps);
		strcpy(base,"/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/output.1.1");
	}
#else
	const char *snpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/short_6.vcf";
	const char *file_base = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/short_6.sorted";
	const char *knownsnpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/short_snp.list";
	const char *base="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/output_short";
#endif

	char file_name[200];
	sprintf(file_name, "%s%s", file_base, ext);

cout << knownsnpfile << endl << snpfile << endl << file_name << endl << base << endl << endl;

	//read_known_snp_file(knownsnpfile);
	read_snp_file(snpfile);
	read_sam_files(file_name);

	runHMM(base);

	return 0;
}
