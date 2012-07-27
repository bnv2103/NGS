// Usage information
// > ./a.out -snp_rate=<double>, -err_rate=<double>, -coverage=<int>, -max_read_len=<int>, -min_read_len=<int>, -ref_file=<string> -snp_file=<string> -indel_file=<string> -ncells=<int> -out_base=<string> -region=<int>
//
// All parameters are optional. Check is performed for duplicate parametes, but not for inconsistency of values.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define FULL
#define DEBUG

int reads_hap[111] = {1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1};
//long int nbp[23] = {0, 24724971, 24295114, 19950182, 19127306, 18085786, 17089999, 15882142, 14627482, 14027325, 13537473, 13445238, 13234953, 11414298, 10636858, 10033891, 8882725, 7877474, 7611715, 6381165, 6243596, 4694432, 4969143};
long int nbp[23] = {0, 247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432};

struct snp_struct {
	int position;
	char all[2];
};

struct indel_struct {
	int position;
	char reflen;
	char all[2][100];
};

struct read_struct {
	char *str;
	long start;
	int len;
	double del_errors;
	double snp_errors;
};

// By default, it will run only for chr21 if FULL is defined, or a small portion of it if FULL is not defined.
#ifdef FULL
struct snp_struct snps[300000];
struct indel_struct indels[35000];
long int N_BP = 46944323; // chr21
int reads_max = 140000;
#else
struct snp_struct snps[100];
struct indel_struct indels[50];
long int N_BP = 20000;
int reads_max = 1000;
#endif

double def_snp_rate = 0.01;
double def_err_rate = 0.04;
double snp_rate;
double err_rate;

int CLT = 30;
char A[3] = "TCG";
char T[3] = "CGA";
char C[3] = "GAT";
char G[3] = "ATC";

int def_coverage = 10;
int coverage;
int max_read_len = 5000;
int min_read_len = 2000;
int max_ref_len = 250000000;
int ncells = 0;
int region = 21;

char *snp_file, *indel_file;
FILE *fq_file, *log_file;
const char *file_base = "simulated_reads";
const char *ref_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reference/bcm_hg18.fasta";

int read_cmp(const void *r1, const void *r2)
{
	return (((struct read_struct*)r1)->start < ((struct read_struct*)r2)->start ? -1 : (((struct read_struct*)r1)->start == ((struct read_struct*)r2)->start ? 0 : 1) );
}

void simulate(int chrnum);

int main(int argc, char **argv)
{
	snp_rate = def_snp_rate;
	err_rate = def_err_rate;
	coverage = def_coverage;
	snp_file = (char*)malloc(sizeof(char)*200);
	indel_file = (char*)malloc(sizeof(char)*200);

	// override parameters
	if(argc>1) {
		int param = 1;
		int snp_ck = 0, err_ck = 0, max_ck = 0, min_ck = 0, cov_ck = 0, ref_ck = 0, snf_ck = 0, ind_ck = 0, cel_ck = 0, out_ck = 0, reg_ck = 0;
		for(param=1;param<argc;param++) {
			char *token[2];
			token[0] = strtok(argv[param], "=");
			if(token[0]==NULL) {
				printf("Ignoring incorrect parameter #%d\n", param);
			} else if((token[1] = strtok(NULL, "=")) == NULL) {
				printf("Ignoring incorrect parameter #%d\n", param);
			} else {
				if(strcmp(token[0],"-snp_rate")==0) {
					if(snp_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						snp_ck++;
						snp_rate = atof((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-err_rate")==0) {
					if(err_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						err_ck++;
						err_rate = atof((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-min_read_len")==0) {
					if(min_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						min_ck++;
						min_read_len = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-max_read_len")==0) {
					if(max_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						max_ck++;
						max_read_len = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-coverage")==0) {
					if(cov_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						cov_ck++;
						coverage = atoi((const char *)token[1]);
						if(coverage<2) {
							printf ("Insufficient coverage %dX provided. Ignoring..\n",coverage);
							coverage = def_coverage;
						}
					}
				} else if(strcmp(token[0],"-ref_file")==0) {
					if(ref_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						ref_ck++;
						ref_file = token[1];
					}
				} else if(strcmp(token[0],"-snp_file")==0) {
					if(snf_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						snf_ck++;
						snp_file = token[1];
					}
				} else if(strcmp(token[0],"-indel_file")==0) {
					if(ind_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						ind_ck++;
						indel_file = token[1];
					}
				} else if(strcmp(token[0],"-ncells")==0) {
					if(cel_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						cel_ck++;
						ncells = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-out_base")==0) {
					if(out_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						out_ck++;
						file_base = token[1];
					}
				} else if(strcmp(token[0],"-region")==0) {
					if(reg_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						reg_ck++;
						region = atoi((const char *)token[1]);
					}
				} else {
					printf("Invalid parameter %s=%s. Ignoring..\n", token[0],token[1]);
				}
			}
		}
		if(snp_rate>err_rate) {
			printf("SNP rate (default 0.01) should be lower than error rate (default 0.04). Using default values now...\n");
				snp_rate = def_snp_rate;
				err_rate = def_err_rate;
		}
	}

	int chrnum;
	char file_name[200];
	char log_file_name[200];
	char *ext = ".fq";
	char *log = ".log";

	sprintf(file_name, "%s%s", file_base, ext);
	sprintf(log_file_name, "%s%s", file_base, log);
	fq_file = fopen(file_name, "w");
	//log_file = fopen(log_file_name, "w");

	if (fq_file == NULL)
		printf("Error opening output file [%s]\n", (const char *) file_name);
	if (log_file == NULL)
		printf("Error opening output file [%s]\n", (const char *) log_file_name);

	simulate(region);
	fclose(fq_file);
}

void simulate(int chrnum)
{
	int c=0;
	long int i=0;
	N_BP = nbp[chrnum];
	sprintf(snp_file, "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_%d.list",chrnum);
	sprintf(indel_file, "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/indel_%d.list",chrnum);
	long int n_reads = 2 * (nbp[chrnum] * coverage)/(max_read_len + min_read_len);
	struct read_struct *reads = (struct read_struct*)malloc(sizeof(struct read_struct)*n_reads);

	for(c=0;c<n_reads;c++) {
		reads[c].str = NULL;
	}

	int ctr=0;
	char line[100];
	char *ref_s = (char*) malloc(sizeof(char)*max_ref_len);
	
	// Read ref fasta here
	FILE *ref = fopen(ref_file,"rt");
	if(ref == NULL) { printf("Can't open ref file \"%s\"\n", ref_file); exit(1); }

	int read_flag = 0;
	while(1) {
		fgets(line,sizeof(line),ref);
		if(feof(ref))
			break;
		int c_ptr = 0;
		int line_len = strlen(line) - 1;
		if(line[0]=='>') {
			if(read_flag==1) {
				break;
			} else {
				if((line_len==5&&line[4]-48==chrnum)||(line_len==6&&(line[4]-48)*10+line[5]-48==chrnum))
					read_flag = 1;
					continue;
			}
		}
		if(read_flag==0)
			continue;
		for(c_ptr = 0;c_ptr<line_len;c_ptr++) {
			ref_s[ctr+c_ptr] = line[c_ptr];
		}
		ctr += line_len;
	}
	ref_s[ctr] = '\0';
	fclose(ref);

	// Read snp file here
	FILE *snp_f = fopen(snp_file,"rt");
	if(snp_f == NULL) { printf("Can't open snp file \"%s\"\n", snp_f); exit(1); }

	int sctr = 0;
	while(1) {
		int pos=0,hap=0;
		char ref,alt;

		fgets(line,sizeof(line),snp_f);
		if(feof(snp_f))
			break;
		strtok(line," ");
		snps[sctr].position = atoi(strtok(NULL," "));
		strtok(NULL," ");
		strtok(NULL," ");
		strtok(NULL," ");
		strtok(NULL," ");
		snps[sctr].all[0] = strtok(NULL," ")[0];
		snps[sctr].all[1] = strtok(NULL," ")[0];
		sctr++;
	}
	fclose(snp_f);

	// Read indel file here
	FILE *indel_f = fopen(indel_file,"rt");
	if(indel_f == NULL) { printf("Can't open snp file \"%s\"\n", indel_f); exit(1); }

	int ictr = 0;
	while(1) {
		int pos=0,hap=0,sit=0;
		char sptr[100];

		fgets(line,sizeof(line),indel_f);
		line[strlen(line)-1] = '\0';
		if(feof(indel_f))
			break;
		strtok(line," ");
		indels[ictr].position = atoi(strtok(NULL," "));
		memset(sptr,'\0',100);
		strcpy(sptr,strtok(NULL," "));
		indels[ictr].reflen = strlen(sptr);
		strtok(NULL," ");
		strtok(NULL," ");
		strtok(NULL," ");
		memset(indels[ictr].all[0],'\0',100);
		memset(indels[ictr].all[1],'\0',100);
		strcpy(indels[ictr].all[0],strtok(NULL," ")); // ref and alt are misnomers here. they're actually the two alleles
		strcpy(indels[ictr].all[1],strtok(NULL," "));
		ictr++;
	}
	fclose(indel_f);

	srand ( (unsigned)time ( NULL ) );
	for(i=0;i<n_reads;i++) {
#ifdef FULL
		int rnd = rand();
		reads[i].start = (int)((rnd*(N_BP-max_read_len))/RAND_MAX) + 1;
#else
		int read_position_start = 9880000;
		int rnd = rand();
		reads[i].start = (int)((rnd*N_BP)/RAND_MAX) + read_position_start;
#endif
	}
	for(i=0;i<n_reads;i++) {
		int rnd = rand();
		reads[i].len = (int)((rnd*(long int)(max_read_len-min_read_len))/RAND_MAX) + min_read_len;
		reads[i].str = (char*)malloc(sizeof(char)*(reads[i].len+1));
	}
	for(i=0;i<n_reads;i++) {
		int clt;
		double sum=0, ssum=0;
		for(clt=0;clt<CLT;clt++) {
			sum += (((rand()*(err_rate*1000))/RAND_MAX) + (err_rate*1000/2));
			ssum += (((rand()*(snp_rate*1000))/RAND_MAX) + (snp_rate*1000/2));
		}
		reads[i].del_errors = sum/(CLT*10);
		reads[i].snp_errors = ssum/(CLT*10);
		sum = 0; ssum = 0;
	}
#ifdef DEBUG
printf("Before sorting\n");
for(i=0; i<n_reads; i++) {
//	printf("Read %d starts at %d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif
	qsort(reads, n_reads, sizeof(struct read_struct), read_cmp);

#ifdef DEBUG
printf("After sorting\n");
for(i=0; i<n_reads; i++) {
//	printf("Read %d starts at %d,%d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].start+reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif

	int snp_ct_start = 0, indel_ct_start = 0;
	char *qualstr = (char *)malloc(sizeof(char) * max_read_len);

	for(i=0;i<n_reads;i++) {
		int j, k = 0;

		int j_st = 0, k_st = 0;
		int snp_flag = 0, indel_flag = 0;
		for(j=0;j<reads[i].len;j++) {
			int homon_hash = 100;
			int snp_ct = 0, common_snp = 0;
			int ctj_start = snp_ct_start > j_st ? snp_ct_start : j_st;
			for(snp_ct = ctj_start; snp_ct <= ctr; snp_ct++) {
				common_snp = 0;
				if(snps[snp_ct].position < reads[i].start) {
					if(snp_flag==0) {
						snp_ct_start = snp_ct+1;
						snp_flag = 1;
					}
				} else if(snps[snp_ct].position == reads[i].start+j+1) {
					j_st = snp_ct+1;
					reads[i].str[k] = snps[snp_ct].all[reads_hap[i%100]-1];
					common_snp = 1;
#ifdef DEBUG
//printf("Replacing common snp %c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),snps[snp_ct].all[reads_hap[i%100]-1],i,k+1,j+1,reads[i].start+j+1);
#endif
					break;
				} else {
					break;
				}
			}

			int indel_ct = 0, common_indel = 0;
			int ctk_start = indel_ct_start > k_st ? indel_ct_start : k_st;
			for(indel_ct = ctk_start; indel_ct <= ctr; indel_ct++) {
				common_indel = 0;
				if(indels[indel_ct].position < reads[i].start) {
					if(indel_flag==0) {
						indel_ct_start = indel_ct+1;
						indel_flag = 1;
					}
				} else if(indels[indel_ct].position == reads[i].start+j+1) {
					k_st = indel_ct+1;
					int inct = 0;
					while(reads[i].str[k++] = (indels[indel_ct].all[reads_hap[i%100]-1])[inct++]);
					k--;
					common_indel = 1;
					//j += indels[indel_ct].reflen - strlen(indels[indel_ct].all[reads_hap[i%100]-1]);
					j += (indels[indel_ct].reflen - 1);
#ifdef DEBUG
printf("Indeling site %c with %s on read %d up to position %d,%d,%d\n",*(ref_s+reads[i].start+j),indels[indel_ct].all[reads_hap[i%100]-1],i,k,j+1,reads[i].start+j+1);
#endif
					break;
				} else {
					break;
				}
			}
			if(common_indel==1) continue;

			double check = (100*(double)rand())/RAND_MAX;
			if(check > reads[i].del_errors+1.0) {
				if(common_snp==0) {
					reads[i].str[k] = *(ref_s+reads[i].start+j);
#ifdef DEBUG
//printf("inserting %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),i,k+1,j+1,reads[i].start+j+1);
#endif
				}
				k++;
			} else if(check <= reads[i].snp_errors) {
				if(*(ref_s+reads[i].start+j) == 'A'||*(ref_s+reads[i].start+j) == 'a') {
					reads[i].str[k] = A[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'T'||*(ref_s+reads[i].start+j) == 't') {
					reads[i].str[k] = T[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'C'||*(ref_s+reads[i].start+j) == 'c') {
					reads[i].str[k] = C[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'G'||*(ref_s+reads[i].start+j) == 'g') {
					reads[i].str[k] = G[(int)((long int)3*rand()/RAND_MAX)];
				} else { reads[i].str[k] = *(ref_s+reads[i].start+j); }
				k++;
#ifdef DEBUG
//printf("replacing %c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),reads[i].str[k-1],i,k,j+1,reads[i].start+j+1);
#endif
			} else {
#ifdef DEBUG
//printf("Deleting %c on read %d\t%d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),i,reads[i].start,k+1,j+1,reads[i].start+j+1);
#endif
			}
		}
		reads[i].str[k] = '\0';
		int l = 0;
		for(l=0;l<strlen(reads[i].str);l++) {
			qualstr[l] = '>';
		}
		qualstr[l] = '\0';

		fprintf(fq_file, "@Start_POS:%d:%d\n", reads[i].start+1,reads_hap[i%100]);
		fprintf(fq_file, "%s\n", reads[i].str);
		fprintf(fq_file, "+\n");
		fprintf(fq_file, "%s\n", qualstr);
//		fprintf(log_file,"%d\t%d\n",startsi,reads_hap[i%100]);
	} // end n_reads
	for(i=0;i<n_reads;i++) {
		free(reads[i].str);
	}
	free(qualstr);
	free(ref_s);
}

