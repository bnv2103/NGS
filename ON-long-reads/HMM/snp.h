
#include<vector>

class READ;

class SNP {

public:

  SNP(void);
  SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string>, bool known);

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

private:
  long position;
  char ref;
  char alt;
  int refcount;
  int altcount;
  int nrefcount;
  int count;
  READ **reads;
  double *gl;
  bool known;

void add_read(int type, READ *read);

};

//===============================================================================
//===============================================================================
