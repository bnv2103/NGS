
class SNP;

class READ {

public:
  READ(void);
  READ(const READ&);
  READ(long start, int length);
  READ (int i);
  READ (double i);

  ~READ(void);

  void addsnp(SNP *snp, char allele);
  void assignHaplotype(int haplotype, double prob);
  long GetPos(void);
  int GetLen(void);
  int GetHap(void);
  double GetHapProb(void);
  SNP *GetSnp(int pos);
  char GetAllele(int pos);
  int GetSnpCount(void);
  SNP **GetSnpList();

  READ* operator=(int i);
  READ* operator+=(double i);
  READ* operator-=(double i);
  READ* operator*=(double i);
  READ* operator/=(double i);
  READ* operator+=(READ* read);
  READ* operator-=(READ* read);
  READ* operator*=(READ* read);
  READ* operator/=(READ* read);
  operator char();
  operator int();
  operator double();
  double operator*(READ* read);
  double operator-(double i);
  double operator+(double i);

private:
  int length;
  int snp_count;
  int hap;
  long start;
  double hprob;
  SNP **snps;
  char *alleles;
};

