
#include<vector>

class READ;

class SNP {

public:

  SNP(void);
  //SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string>, bool known);
  SNP(long snp_pos, char snp_ref, char snp_alt, vector<string>, int known, double qual);

  ~SNP(void) {}

  long GetPos();
  void append(int type, READ *read);
  READ *GetRead(int pos);
  char GetRef();
  char GetAlt();
  int GetReadCount();
  int GetRefCount();
  int GetAltCount();
  int GetErrCount();
  double* GetGenLik();
  int GetKnown();
  double GetQualScore();
  void add_posteriors(double posterior[3]);
  void assign_genotype(int gt, double genp);
  void PrintPosterior();
  void PrintLR();
  double* GetPosteriors();

private:
  int known;
  char ref;
  char alt;
  int refcount;
  int altcount;
  int errcount;
  int count;
  int posterior_count;
  int genotype;
  long position;
  double qual;
  double likelihood_ratio;
  double genprob;
  READ **reads;
  double *gl;
  double **posteriors;
  double *posterior;

void add_read(int type, READ *read);
void CalculateLikelihoodRatio();

};

//===============================================================================
//===============================================================================
