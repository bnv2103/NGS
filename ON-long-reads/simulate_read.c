// Usage information
// > ./a.out -snp_rate=<double>, -err_rate=<double>, -coverage=<int>, -max_read_len=<int>, -min_read_len=<int>, -ref_file=<string> -snp_file=<string> -ncells=<int> -out_base=<string>
//
// All parameters are optional. Check is performed for duplicate parametes, but not for inconsistency of values.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define FULL
//#define SIMULATION
#define DEBUG
int reads_hap[250] = {1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1};

struct snp_struct {
	int position;
	int hap;
	char ref;
	char alt;
};

struct read_struct {
	char *str;
	long start;
	int len;
	double del_errors;
	double snp_errors;
};

#ifdef FULL
struct snp_struct snps[81000];
#else
struct snp_struct snps[100];
#endif

//int read_cmp(const struct read_struct *r1, const struct read_struct *r2)
int read_cmp(const void *r1, const void *r2)
{
	return (((struct read_struct*)r1)->start < ((struct read_struct*)r2)->start ? -1 : (((struct read_struct*)r1)->start == ((struct read_struct*)r2)->start ? 0 : 1) );
}

int main(int argc, char **argv)
{
	double snp_rate = 0.01;
	double err_rate = 0.04;
	double def_snp_rate = snp_rate;
	double def_err_rate = err_rate;

	int CLT = 30;
	char A[3] = "TCG";
	char T[3] = "CGA";
	char C[3] = "GAT";
	char G[3] = "ATC";

#ifdef FULL
	long int N_BP = 46944323;
	int reads_max = 140000;
	const char *snp_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/common_snps_21.txt";
#else
	long int N_BP = 20000;
	int reads_max = 1000;
	const char *snp_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/common_snps_21_short.txt";
#endif
	int coverage = 10;
	int def_coverage = coverage;
	int max_read_len = 5000;
	int min_read_len = 2000;
	int max_ref_len = 50000000;
	int ncells = 0;
	struct read_struct reads[reads_max];

	int c=0;
	for(c=0;c<reads_max;c++) {
		reads[c].str = NULL;
	}

	//REVISIT: Phased snps needed here
	//const char *snp_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/phased_indi_snps.txt";
	const char *ref_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/bcm_hg18_chr21.fasta";
	const char *file_base = "simulated_reads";

	// override parameters
	if(argc>1) {
		int param = 1;
		int snp_ck = 0, err_ck = 0, max_ck = 0, min_ck = 0, cov_ck = 0, ref_ck = 0, snf_ck = 0, cel_ck = 0, out_ck = 0;
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

	long int i=0;
	long int n_reads = (2 * N_BP * coverage)/(max_read_len+min_read_len);
	
	// Read ref fasta here
	FILE *ref = fopen(ref_file,"rt");
	if(ref == NULL) { printf("Can't open ref file \"%s\"\n", ref_file); exit(1); }

	// Read snp file here
	// cut -d' ' -f 4,181 genotypes_chr21_CEU_r27_nr.b36_fwd.txt | awk -F "" '{if(NR!=1) {for(i=1;i<=NF;i++) {if(i!=(NF-1)) printf $i;}printf "\n";}}' > phased_indi_snps.txt 
	FILE *snp_f = fopen(snp_file,"rt");
	if(snp_f == NULL) { printf("Can't open snp file \"%s\"\n", snp_f); exit(1); }

	int ctr=0;
	char line[100];
	char *ref_s = (char*) malloc(sizeof(char)*max_ref_len);
	
	while(fgets(line,sizeof(line),ref)) {
		int c_ptr = 0;
		int line_len = strlen(line) - 1;
		for(c_ptr = 0;c_ptr<line_len;c_ptr++) {
			ref_s[ctr+c_ptr] = line[c_ptr];
		}
		ctr += line_len;
	}
	ref_s[ctr] = '\0';
	fclose(ref);

	// Add snps, indels to ref file
#ifdef OLD
	while(1) {
		int pos = 0, r = 0;
		char allele;
		if(fgets(line,sizeof(line),snp_f)) {
			pos = atoi(line);
			allele = line[strlen(line)-2];
		} else {
			break;
		}
		ref_s[pos-1] = allele;
	}
#endif
	ctr = 0;
	while(1) {
		int pos=0,hap=0;
		char ref,alt;

		// 10% of these snps will be made novel snps with low prior values in the HMM.
		// How we infer these 10% for now is basically every 10th such snp.
		if(fgets(line,sizeof(line),snp_f)) {
			pos = atoi(line);
			ref = line[strlen(line)-4];
			alt = line[strlen(line)-2];
			hap = ctr%2+1;
			snps[ctr].position = pos;
			snps[ctr].hap = hap;
			snps[ctr].ref = ref;
			snps[ctr].alt = alt;
		} else
			break;
		ctr++;
	}

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
#ifdef SIMULATION
		reads[i].start = i*((N_BP-max_read_len)/n_reads) + 1;
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
		//sum += 15;
		reads[i].del_errors = sum/(CLT*10);
		reads[i].snp_errors = ssum/(CLT*10);
//printf("%f\t%f\n",ssum,reads[i].snp_errors);
		sum = 0; ssum = 0;
	}
#ifdef DEBUG
printf("Before sorting\n");
for(i=0; i<n_reads; i++) {
	printf("Read %d starts at %d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif
	qsort(&reads, n_reads, sizeof(struct read_struct), read_cmp);

#ifdef DEBUG
printf("After sorting\n");
for(i=0; i<n_reads; i++) {
	printf("Read %d starts at %d,%d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].start+reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif
	// RAND_MAX = 2147483647
	char file_name[30];
	char log_file_name[30];
	char *ext = ".fq";
	char *log = ".log";
	sprintf(file_name, "%s_%d%s", file_base, ncells, ext);
	sprintf(log_file_name, "%s_%d%s", file_base, ncells, log);
	FILE *fq_file = fopen(file_name, "w");
	FILE *log_file = fopen(log_file_name, "w");

	if (fq_file == NULL)
		printf("Error opening output file [%s]\n", (const char *) file_name);
	if (log_file == NULL)
		printf("Error opening output file [%s]\n", (const char *) log_file_name);

	char *qualstr = (char *)malloc(sizeof(char) * max_read_len);

	int ct_start = 0;
//	srand ((unsigned)time(NULL));
	for(i=0;i<n_reads;i++) {
		int j, k = 0;

/* REVISIT: This piece deletes any common snps as well. Skipping for prototype.
		for(j=0;j<lensi;j++) {
			int ct = 0, common_snp = 0;
			for(ct=0;ct<=ctr;ct++) {
				common_snp = 0;
				if(snps[ct].position==startsi+j+1) {
//printf("Replace\n");
					if(reads_hap[i]==snps[ct].hap) {
						readsi[k] = snps[ct].alt;
						common_snp = 1;
printf("Replacing common snp %c/%c with %c on read %d at position %d,%d,%d\n",*(ref_s+startsi+j),snps[ct].ref,snps[ct].alt,i,k+1,j+1,startsi+j+1);
						break;
					}
				}
			}

			double check = (100*(double)rand())/RAND_MAX;
			if(check > del_errorsi) {
				if(common_snp==0) {
					readsi[k] = *(ref_s+startsi+j);
printf("inserting %c on read %d at position %d,%d,%d\n",*(ref_s+startsi+j),i,k,j+1,startsi+j+1);
				}
				k++;
			} else if(check <= snp_errorsi) {
				if(*(ref_s+startsi+j) == 'A'||*(ref_s+startsi+j) == 'a') {
					readsi[k] = A[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+startsi+j) == 'T'||*(ref_s+startsi+j) == 't') {
					readsi[k] = T[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+startsi+j) == 'C'||*(ref_s+startsi+j) == 'c') {
					readsi[k] = C[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+startsi+j) == 'G'||*(ref_s+startsi+j) == 'g') {
					readsi[k] = G[(int)((long int)3*rand()/RAND_MAX)];
				} else { readsi[k] = *(ref_s+startsi+j); }
				k++;
printf("replacing %c with %c on read %d at position %d,%d,%d\n",*(ref_s+startsi+j),readsi[k-1],i,k,j+1,startsi+j+1);
			} else {
printf("Deleting %c on read %d at position %d,%d,%d\n",*(ref_s+startsi+j),i,k+1,j+1,startsi+j+1);
			}
		}
*/


// REVISIT: This piece does not delete common snps. Using for prototype only
		int j_st = 0;
		int flag = 0;
		for(j=0;j<reads[i].len;j++) {
			int homon_hash = 100;
			int ct = 0, common_snp = 0;
			int ctj_start = ct_start > j_st ? ct_start : j_st;
			for(ct = ctj_start; ct <= ctr; ct++) {
				common_snp = 0;
				if(snps[ct].position < reads[i].start) {
					if(flag==0) {
						ct_start = ct+1;
						flag = 1;
					}
				} else if(snps[ct].position == reads[i].start+j+1) {
					j_st = ct+1;
					// Every 100th common snp is made homozygous null reference
					if(reads_hap[i%100] == snps[ct].hap || (ct%homon_hash)+1==homon_hash) {
						reads[i].str[k] = snps[ct].alt;
						common_snp = 1;
						k++;
#ifdef DEBUG
printf("Replacing common snp %c/%c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),snps[ct].ref,snps[ct].alt,i,k+1,j+1,reads[i].start+j+1);
#endif
					}
					break;
				// Experimental: Looks good
				//} else if(snps[ct].position > reads[i].start + reads[i].len) {
				} else {
					break;
				}
			}
			if(common_snp==1) continue;

			double check = (100*(double)rand())/RAND_MAX;
			if(check > reads[i].del_errors) {
				reads[i].str[k] = *(ref_s+reads[i].start+j);
#ifdef DEBUG
printf("inserting %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),i,k,j+1,reads[i].start+j+1);
#endif
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
printf("replacing %c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),reads[i].str[k-1],i,k,j+1,reads[i].start+j+1);
#endif
			} else {
#ifdef DEBUG
printf("Deleting %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),i,k+1,j+1,reads[i].start+j+1);
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
	fclose(fq_file);
	fclose(log_file);
	free(qualstr);
	free(ref_s);
}

