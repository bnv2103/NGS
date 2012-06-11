
#include<vector>

class READ;

class SNP {

public:

  SNP(void);
  //SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string>, bool known);
  SNP(long snp_pos, char snp_ref, char snp_alt, vector<string>, bool known);

  ~SNP(void) {}

  long GetPos();
  void append(int type, READ *read);
  READ *GetRead(int pos);
  char GetRef();
  char GetAlt();
  int GetReadCount();
  int GetRefCount();
  double *GetGenLik();
  bool GetKnown();
  void add_posteriors(double posterior[3]);
  void PrintPosterior();
  void PrintLR();

private:
  bool known;
  char ref;
  char alt;
  int refcount;
  int altcount;
  int nrefcount;
  int count;
  int posterior_count;
  long position;
  double likelihood_ratio;
  READ **reads;
  double *gl;
  double **posteriors;
  double *posterior;

void add_read(int type, READ *read);
void CalculateLikelihoodRatio();

};

//===============================================================================
//===============================================================================
