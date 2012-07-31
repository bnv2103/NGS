/* ************************************************************************ *
 * ************************************************************************ *

   File: hmm.h
   The class CHMM defines operations for HMM

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   Date:  2-18-99 

 * ************************************************************************ *

   Modification Log:
	

 * ************************************************************************ *
   Log for new ideas:
 * ************************************************************************ *
               Language and Media Processing
               Center for Automation Research
               University of Maryland
               College Park, MD  20742
 * ************************************************************************ *
 * ************************************************************************ */
 
//===============================================================================

//===============================================================================
//===============================================================================

//===============================================================================

const double DELTA = 5*1E-3;// threshold used for exiting Baum-Welch or K-Means loop
//===============================================================================

class CHMM{

public:

	CHMM(CStateTrans *a, CObsProb *b, CInitStateProb *pi);
//      CHMM(char* hmmFileName);
	~CHMM(void);
	double Forward(double **alpha, double *scale, CObs **obs, long T, boolean doLog);
	void Backward(double **beta, double *scale, CObs **obs, long T);
	double ForwardAlgo(double **alpha, double *scale, CObs **obs, long T, boolean doLog);
	void BackwardAlgo(double **beta, double *scale, CObs **obs, long T);
        double Viterbi(CObs **obs, long T, int *q);
        double ViterbiLog(CObs **obs, long T, int *q);
        double ViterbiLog(CObs **obs, long T, int *q, double *prob);
	double BaumWelchCore(CObs **obs, long T, double *gamma, double **xi, boolean doLog);
        double IterBaumWelch(CObsSeq *obsSeq, double *gamma, double **xi);
	void LearnBaumWelch(CObsSeq *obsSeq);
	void FindFBDistance(CObsSeq *obsSeq, ostream &outFile, long start, long end);
        double SegmentalKMeansCore(CObs **obs, long T);
        double IterSegmentalKMeans(CObsSeq *obsSeq);
        void LearnSegmentalKMeans(CObsSeq *obsSeq);
        void LearnHybridSKM_BW(CObsSeq *obsSeq);
        CObsSeq* GenerateSequences(long nbSequences, long nbObs, int seed);
        CObs** GenerateObservations(long T);
 #if 1
        CObsSeq* ReadSequences(ifstream &inputSeqFile);
 #endif
	void PrintStatesAndExpectedObs(CObsSeq *obsSeq, 
						 ostream& stateFile, ostream& bestObsFile);
        double FindDistance(CObsSeq *obsSeq, ostream &outFile);
	//GLOBAL
        //double FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile, vector<READ*> *reads_list, vector<SNP*> *snp_list);
        double FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile, ostream &gtFile);
        double FindCrossEntropyDistance(CObsSeq *obsSeq, ostream &outFile);
        double FindQToPiProb(long T, int *q);
        void Print(ostream &outFile);
        int GetN(void){return mN;};
	void GetCommonSnpList(CObs **obs, SNP** reads_snp_list, int *common_snp_count, int *index, int t);
	double compute_g_h(SNP **reads_snp_list, int count, CObs **obs, int t, int *index, int hap, int type);
	double compute_c_g(SNP **reads_snp_list, int count, CObs **obs, int t, int *index, int hap, int type);
	double compute_r_s(SNP **reads_snp_list, int common_snp_count, CObs **obs, int t, int *index, int hap);
	double binomial(int r, int n, double p);
	unsigned long long com (int n, int m);
	double compute_new_emission(int obs, int type, int hap);
//	double compute_new_emission(SNP **reads_snp_list, int count, CObs **obs, int t, int *index, int hap);
	double compute_new_emission(SNP**reads_snp_list, int count, CObs **obs, int t, int *index, int hap, double obslik[3], double genlik[3]);
	void haplotypeProbability(vector<SNP*>::iterator snp_it, double happ[3]);
	void genotypeProbability(vector<SNP*>::iterator snp_it, double genp[2]);
	//double* haplotypeProbability(SNP *snp_it);
	//GLOBAL
	//void UpdateGenotypes(vector<SNP*> *snp_list);
	void UpdateGenotypes(long start, long end);

protected:
        int mN;// nb of states
	CStateTrans *mA;
	CObsProb *mB;
	CInitStateProb *mPi;
};

//===============================================================================
//===============================================================================
