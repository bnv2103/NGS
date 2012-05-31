
/* ************************************************************************ *
 * ************************************************************************ *

   File: hmm.C
   The class CHMM defines operations for HMM

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   building up on original C code by Tapas Kanungo
   Date:  2-18-99 

 * ************************************************************************ *

   Modification Log:
	4-14-99: Compute log(A) and log(pi) in respective classes
	4-16-99: Compute ViterbiLog with state duration probabilities

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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

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
#include "stateTrans.h"
#include "plainStateTrans.h"
#include "gammaProb.h"
#include "explicitDurationTrans.h"
#include "initStateProb.h"
#include "hmm.h"

//===============================================================================

#define DEBUG

//===============================================================================

//===============================================================================

CHMM::CHMM(CStateTrans *a, CObsProb *b, CInitStateProb *pi)
{
	mA = a;
	mB = b;
	mPi = pi;

	mN = mA->GetN();
}

//===============================================================================

CHMM::~CHMM(void)
{
// Nothing for now
}

//===============================================================================

double CHMM::Forward(double **alpha, double *scale, CObs **obs, long T, boolean doLog)
     // Scaling is used to prevent roundoff errors
     // Same scaling is used for backward and forward procedures
     // so that the scales cancel out in the Baum Welch formula
     // Quantity returned is actually - log(P(O | model),
     // i.e. exponential of this quantity is 1 / P(O | model)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double logProb;
	double bi1, aij, bjt1;

// 1. Initialization
	for (i = 1; i <= mN; i++) {
	  bi1 = mB->at(i, obs[1]);
	  alpha[1][i] = mPi->at(i) *  bi1;
	}
	scale[1] = Normalize(alpha[1], mN);

// 2. Induction
	for (t = 1; t <= T - 1; t++) {
	  for (j = 1; j <= mN; j++) {
	    sum = 0.0;
	    for (i = 1; i <= mN; i++){
	      aij = mA->at(i, j);
	      sum += alpha[t][i] * aij;
	    }
	    bjt1 =  mB->at(j, obs[t+1]);
	    alpha[t+1][j] = sum * bjt1;
	  }
	  scale[t+1] = Normalize(alpha[t+1], mN);
	}
	logProb = 0.0;

// 3. Termination
	if(doLog){
	  for (t = 1; t <= T; t++){
	    logProb += log(scale[t]);
	  }
	}// endif

	return logProb;// zero returned if doLog is false
}

//===============================================================================

void CHMM::Backward(double **beta, double *scale, CObs **obs, long T)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
// 1. Initialization
		for (i = 1; i <= mN; i++){
			beta[T][i] = 1.0/scale[T];
		}
 
// 2. Induction
     for (t = T - 1; t >= 1; t--){
	  for (i = 1; i <= mN; i++){
	    sum = 0.0;
	    for (j = 1; j <= mN; j++){
	      sum += mA->at(i,j) * mB->at(j, obs[t+1]) * beta[t+1][j];
	    }
	    beta[t][i] = sum/scale[t];
	    }
     }
}

//===============================================================================

double CHMM::Viterbi(CObs **obs, long T, int *q)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
{
	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double	maxval, val;
	double pprob;
	double *prevDelta, *delta, *tmp;
	int **psi;
	
	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

//1. Initialization
	for (i = 1; i <= mN; i++) {
		prevDelta[i] = mPi->at(i) * mB->at(i, obs[1]);
		psi[1][i] = 0;
	}	

// 2. Recursion
	for (t = 2; t <= T; t++) {
		for (j = 1; j <= mN; j++) {
			maxval = 0.0;
			argmaxval = 1;	
			for (i = 1; i <= mN; i++) {
				val = prevDelta[i] * mA->at(i, j);
				if (val > maxval) {
					maxval = val;	
					argmaxval = i;	
				}
			}
			delta[j] = maxval * mB->at(j, obs[t]);
			psi[t][j] = argmaxval; // reverse order of indices?
		}
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}

// 3. Termination
	pprob = 0.0;
	q[T] = 1;
	for (i = 1; i <= mN; i++) {
	  if (prevDelta[i] > pprob) {
			pprob = prevDelta[i];	
			q[T] = i;
	  }
	}

// 4. Path (state sequence) backtracking

	for (t = T - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete delta;
	
	return pprob;
}

//===============================================================================

double CHMM::ViterbiLog(CObs **obs, long T, int *q)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
// This implementation uses logarithms to avoid underflows.
{

	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double	maxval, val, bVal, logVal;
	double logProb;
	double *prevDelta, *delta, *tmp;
	int **psi;
	double **logBiOt;
	int zeroProbCount;

	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

// We do not preprocess the logs for B
// because the data have variable length T
	logBiOt =  SetMatrix(mN, T);
	for (t = 1; t <= 1; t++){
            zeroProbCount = 0;
            for (i = 1; i <= mN; i++){ 
              bVal = mB->at(i, t); // obs[t] is an entire seq of snps in a read
                if(bVal<=0.0){
                    logVal = -1000.0;
                    zeroProbCount++;
                }
                else{
                    logVal = log(bVal);
                }
                logBiOt[i][t] = logVal;
            }// for i
            if(zeroProbCount == mN){// unseen obs, Viterbi decides which seen obs to use
                cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                for (i = 1; i <= mN; i++){		  
                    logBiOt[i][t] = log(1.0/mN);
                }
            }
 	}// for t
	// May have to move this closing brace to the end or rather
	// perform this log conversion in the recursion step.
// 1. Initialization
	for (i = 1; i <= mN; i++){
		prevDelta[i] = mPi->logAt(i) + logBiOt[i][1];
		psi[1][i] = 0; // What's this for?
		mA->InitViterbiDurations(i); // NO-OP
	}
 
// 2. Recursion
	for (t = 2; t <= T; t++) { // successive reads
		zeroProbCount = 0;
		int common_snp_count = 0;
		for (j = 1; j <= mN; j++) {
			maxval = prevDelta[1] + mA->logAt(1, j);
			argmaxval = 1;
			for (i = 1; i <= mN; i++) {// previous state
				val = prevDelta[i] + mA->logAt(i, j);
				if (val > maxval) {
					maxval = val;
					argmaxval = i;
				}
			}

			SNP *reads_snp_list[250];
			int *index = new int[500];
			common_snp_count = 0;
			int hap = (logBiOt[j][t-1] > logBiOt[(j%mN)+1][t-1]) ? 1 : 2; // 1 -> H1==H2, 2 -> H1!=H2

			GetCommonSnpList(obs, reads_snp_list, &common_snp_count, index, t);

			// I end up not assigning the actual emission matrix, but only the log matrix for the computations

			logBiOt[j][t] = 1.0;
			for(int count=0; count<common_snp_count; count++) {
				double emission = compute_new_emission(reads_snp_list, count, obs, t, index, hap);
				if(emission <= 0.0) {
					logVal = -1000.0;
                    			zeroProbCount++;
                		} else{
                    			logVal = log(emission);
                		}
                		logBiOt[j][t] *= logVal;
			}
/*
			logBiOt[j][t] = 0.0;
			for(int count=0; count<common_snp_count; count++) {
				double g_h = 0.0, c_g = 0.0, c_h = 0.0;
				for(int type=0; type<3; type++) {
					g_h += compute_g_h(reads_snp_list, count, obs, t, index, hap, type+1);
					c_g += compute_c_g(reads_snp_list, count, obs, t, index, hap, type+1);
					c_h += g_h * c_g;
				double emission = compute_emission(reads_snp_list, count, obs, t, index, hap);
				}
				if(c_h <= 0.0) {
					logVal = -1000.0;
                    			zeroProbCount++;
                		} else{
                    			logVal = log(c_h);
                		}
                		logBiOt[j][t] += logVal;
			}
			 logBiOt[j][t] /= common_snp_count;
			 logBiOt[j][t] *= compute_r_s(reads_snp_list, common_snp_count, obs, t, index, hap);
*/

			delta[j] = maxval + logBiOt[j][t]; 
			psi[t][j] = argmaxval; // What's this for?
			mA->UpdateViterbiDurations(argmaxval, j);// ***dfd 4-16-99
		}
		if(zeroProbCount > common_snp_count) { // need to check per snp
                	cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                	for (i = 1; i <= mN; i++){		  
                    		logBiOt[i][t] = log(1.0/mN);
                	}
		}
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}
 
// 3. Termination
	logProb = prevDelta[1];
	q[T] = 1;
	for (i = 1; i <= mN; i++) {
		if (prevDelta[i] > logProb) {
		      logProb = prevDelta[i];
		      q[T] = i;
		}
	} 
 
// 4. Path (state sequence) backtracking
	for (t = T - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete [] delta;
	delete [] logBiOt[1];
	delete [] logBiOt;
#if 0
	delete [] logA[1];
	delete [] logA;
	delete [] logPi;
#endif	
	return logProb;
}

double CHMM::ViterbiLog(CObs **obs, long T, int *q, double *probarray)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
// This implementation uses logarithms to avoid underflows.
{

	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double  prob;
	double	maxval, val, bVal, logVal;
	double logProb;
	double *prevDelta, *delta, *tmp;
	int **psi;
	double **logBiOt;
	int zeroProbCount;

	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

// We do not preprocess the logs for B
// because the data have variable length T
	logBiOt =  SetMatrix(mN, T);
	for (t = 1; t <= 1; t++){
            zeroProbCount = 0;
            for (i = 1; i <= mN; i++){ 
              bVal = mB->at(i, t); // obs[t] is an entire seq of snps in a read
                if(bVal<=0.0){
                    logVal = -100.0;
                    zeroProbCount++;
                }
                else{
                    logVal = log(bVal);
                }
                logBiOt[i][t] = logVal;
            }// for i
            if(zeroProbCount == mN){// unseen obs, Viterbi decides which seen obs to use
                cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                for (i = 1; i <= mN; i++){		  
                    logBiOt[i][t] = log(1.0/mN);
                }
            }
probarray[t] = -10000;
 	}// for t

// 1. Initialization
// Initialization is performed for prevDelta and psi only.
// PrevDelta stores probability to stage k-1 for all states

	for (i = 1; i <= mN; i++){
		prevDelta[i] = mPi->logAt(i) + logBiOt[i][1];
		psi[1][i] = 0; // What's this for?
		if(probarray[1] < logBiOt[i][1]) probarray[1] = logBiOt[i][1];
		mA->InitViterbiDurations(i); // NO-OP
	}
 
// 2. Recursion
	int nread = 2;
	for (t = 2; t <= T; t++) { // successive reads
		zeroProbCount = 0;
		int common_snp_count = 0;
		SNP *reads_snp_list[250];
		int *index = new int[500];

		GetCommonSnpList(obs, reads_snp_list, &common_snp_count, index, t);

		//REVISIT: There is a bug here. In case two consecutive reads with 0 overlapping snps are encountered, it may not
		//continue to work as expected. We might have to reset haplotype assumptions, or compare the last read with overlapping
		//snps to the next such one (again to which there might be very little chance)
		if(common_snp_count==0) {
			cout << "Read " << t-1 << " and " << t << " have no overlapping snps. Skipping to the next pair.." << endl;
			continue;
		}
		for (j = 1; j <= mN; j++) {
			// prevDelta[i] = v(k)(i)
			maxval = prevDelta[1] + mA->logAt(1, j);
			argmaxval = 1;
			prob = maxval;
			for (i = 1; i <= mN; i++) {// previous state
				val = prevDelta[i] + mA->logAt(i, j);
				if (val > maxval) {
					maxval = val;
					argmaxval = i;
					prob = maxval;
				}
			}
			// maxval = max(k) (v(k)(i)a(k)(l))
			// need to add emission here and find max state for next read

			// I end up not assigning the actual emission matrix, but only the log matrix for the computations

			logBiOt[j][nread] = 0.0;
			// REVISIT: Change this to find hap based on psi/q arrays?
			int hap = prevDelta[j]>prevDelta[(j%mN)+1] ? 1 : 2;
			int abs_hap = prevDelta[j] > prevDelta[(j%mN)+1] ? j : (j%mN)+1;
if(j==9) cout << "Read " << t << " : " << abs_hap << " : " << -prevDelta[j] << " : " << -logBiOt[j][nread-1] << " : " << -prevDelta[j%mN+1] << " : " << -logBiOt[j%mN+1][nread-1] << " : " << endl;
			for(int count=0; count<common_snp_count; count++) {
				double emission = compute_new_emission(reads_snp_list, count, obs, t, index, hap);
//cout << reads_snp_list[count]->GetKnown() << " : " << reads_snp_list[count]->GetPos() << " : " << emission << endl;
				if(emission <= 0.0) {
					logVal = -100.0;
                    			zeroProbCount++;
                		} else{
                    			logVal = log(emission);
                		}
                		logBiOt[j][nread] += logVal;
			}
			// logBiOt[j][nread] now contains the mN different emissions to be added to maxval to determine the max state at this stage

			delta[j] = maxval + logBiOt[j][nread]; 
			psi[nread][j] = argmaxval; // What's this for?
			if(probarray[nread] < logBiOt[j][nread]) probarray[nread] = logBiOt[j][nread];
			mA->UpdateViterbiDurations(argmaxval, j);// ***dfd 4-16-99
		}
		// By now all v(l)(i+1)s are computed from which the max has to be found
		if(zeroProbCount > common_snp_count) { // need to check per snp
                	cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                	for (i = 1; i <= mN; i++){		  
                    		logBiOt[i][nread] = log(1.0/mN);
                	}
		}
		nread++;
		// prevDelta updated to find the max from in the next read
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}
	nread--;
 
// 3. Termination
	logProb = prevDelta[1];
	//q[T] = 1;
	q[nread] = 1;
	for (i = 1; i <= mN; i++) {
		if (prevDelta[i] > logProb) {
		      logProb = prevDelta[i];
		      q[nread] = i;
		      //q[T] = i;
		}
	} 
 
// 4. Path (state sequence) backtracking
	//for (t = T - 1; t >= 1; t--){
	for (t = nread - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete [] delta;
	delete [] logBiOt[1];
	delete [] logBiOt;
#if 0
	delete [] logA[1];
	delete [] logA;
	delete [] logPi;
#endif	
	return logProb;
}

void CHMM::GetCommonSnpList(CObs**obs, SNP** reads_snp_list, int *common_snp_count, int *index, int t)
{
	READ prev_read = ((CFlexibleObs<READ>*)(obs[t-1]))->Get(1);
	int prev_read_snp_count = prev_read.GetSnpCount();
	SNP **prev_snp_list = prev_read.GetSnpList();
	READ curr_read = ((CFlexibleObs<READ>*)(obs[t]))->Get(1);
	int curr_read_snp_count = curr_read.GetSnpCount();
	SNP **curr_snp_list = curr_read.GetSnpList();
	int it1 = 0, it2 = 0;

	while(it1<prev_read_snp_count && it2<curr_read_snp_count) {
		SNP* snp1 = prev_snp_list[it1];
		SNP* snp2 = curr_snp_list[it2];

///*
		if(snp1->GetPos() == snp2->GetPos()) {
			reads_snp_list[*common_snp_count] = snp1;
			index[2*(*common_snp_count)] = it1;
			index[(2*(*common_snp_count))+1] = it2;
			(*common_snp_count)++;it1++;it2++;
		} else if(snp1->GetPos() < snp2->GetPos()) {
			it1++;
//cout << prev_read.GetPos() << ", " << curr_read.GetPos() << endl;
//cout << "It1++:" << snp1->GetPos() << ", " << snp2->GetPos() << endl;
		} else if(snp1->GetPos() > snp2->GetPos()) {
			it2++;
//cout << "It2++" << endl;
		}
//*/
	}

}

//===============================================================================


double CHMM::compute_new_emission(SNP **reads_snp_list, int count, CObs **obs, int t, int *index, int hap)
{
	double err_rate = 0.01;
	double known_snp_rate = 0.1;
	double novel_snp_rate = 0.001;
	double gen_prior[3];
	double gen_priorn[3];
	double snp_rate = novel_snp_rate;
	double gen_lik[3] = {0.9, 0.01, 0.09};
	double obs_lik[3] = {0.0, 0.0, 1.0};
	double prob = 0.0;
	double lik = 0.0;

	char ref = reads_snp_list[count]->GetRef();
	char all1 = ((CFlexibleObs<READ>*)(obs[t-1]))->Get(1).GetAllele(index[2*count]);
	char all2 = ((CFlexibleObs<READ>*)(obs[t]))->Get(1).GetAllele(index[(2*count)+1]);
	int obt = (ref==all1&&ref==all2) ? 1 : (ref!=all1&&ref!=all2) ? 3 : 2; // ref,ref -> 1, ref,nref -> 2, nref,nref -> 3

	if(reads_snp_list[count]->GetKnown())
		snp_rate = known_snp_rate;
	gen_prior[0] = 1 - snp_rate - snp_rate*snp_rate;
	gen_prior[1] = snp_rate;
	gen_prior[2] = snp_rate*snp_rate;
	gen_priorn[0] = 1;
	gen_priorn[1] = 1;
	gen_priorn[2] = 1;

	for(int j=0; j<3; j++) {
		lik += (reads_snp_list[count]->GetGenLik())[j] * gen_prior[j];
//cout << reads_snp_list[count]->GetPos() << " : " << (reads_snp_list[count]->GetGenLik())[j] << endl;
	}

//cout << reads_snp_list[count]->GetKnown() << " : " << reads_snp_list[count]->GetPos();
	for(int i=0; i<3; i++) {
		gen_lik[i] = ((reads_snp_list[count]->GetGenLik())[i] * gen_prior[i])/lik;
//cout << " :\t" << gen_lik[i];

		switch(i) {
		case 0:
			if(obt==1&&hap==1)
				obs_lik[i]=
			else if(obt==1&&hap==2)
				obs_lik[i]=
			else if(obt==2&&hap==1)
				obs_lik[i]=
			else if(obt==2&&hap==2)
				obs_lik[i]=
			else if(obt==3&&hap==1)
				obs_lik[i]=
			else if(obt==3&&hap==2)
				obs_lik[i]=
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		case 1:
			if(obt==1&&hap==1)
				obs_lik[i]=
			else if(obt==1&&hap==2)
				obs_lik[i]=
			else if(obt==2&&hap==1)
				obs_lik[i]=
			else if(obt==2&&hap==2)
				obs_lik[i]=
			else if(obt==3&&hap==1)
				obs_lik[i]=
			else if(obt==3&&hap==2)
				obs_lik[i]=
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		case 2:
			if(obt==1&&hap==1)
				obs_lik[i]=
			else if(obt==1&&hap==2)
				obs_lik[i]=
			else if(obt==2&&hap==1)
				obs_lik[i]=
			else if(obt==2&&hap==2)
				obs_lik[i]=
			else if(obt==3&&hap==1)
				obs_lik[i]=
			else if(obt==3&&hap==2)
				obs_lik[i]=
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		}

		//prob += gen_lik[i] * obs_lik[i];
		prob += gen_lik[i] * obs_lik[i] * gen_priorn[i];
	}
//cout << endl;
//cout << reads_snp_list[count]->GetKnown() << " : " << reads_snp_list[count]->GetPos() << " : " << all1 << " : " << all2 << " : " << prob << endl;
	return prob;
/*
		switch(i) {
			case 0:
		if(obt==1&&hap==1)
			obs_lik[i]=(1-err_rate)*(1-err_rate);
		else if(obt==1&&hap==2)
			obs_lik[i]=0;
		else if(obt==2&&hap==1)
			obs_lik[i]=err_rate*(1-err_rate);
		else if(obt==2&&hap==2)
			obs_lik[i]=0;
		else if(obt==3&&hap==1)
			obs_lik[i]=err_rate*err_rate;
		else if(obt==3&&hap==2)
			obs_lik[i]=0;
		else
			cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
			case 1:
		if(obt==1&&hap==1)
			obs_lik[i]=(1-err_rate)*(1-err_rate); //CONFIRM!!
		else if(obt==1&&hap==2)
			obs_lik[i]=err_rate*(1-err_rate);
		else if(obt==2&&hap==1)
			obs_lik[i]=(1-err_rate)*(1-err_rate); // CONFIRM!! 0/e*e ?
		else if(obt==2&&hap==2)
			obs_lik[i]=(1-err_rate)*(1-err_rate);
		else if(obt==3&&hap==1)
			obs_lik[i]=err_rate*err_rate; //CONFIRM!!
		else if(obt==3&&hap==2)
			obs_lik[i]=err_rate*(1-err_rate);
		else
			cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
			case 2:
		if(obt==1&&hap==1)
			obs_lik[i]=(1-err_rate)*(1-err_rate); //CONFIRM!!
		else if(obt==1&&hap==2)
			obs_lik[i]=err_rate*(1-err_rate);
		else if(obt==2&&hap==1)
			obs_lik[i]=err_rate*(1-err_rate);
		else if(obt==2&&hap==2)
			obs_lik[i]=0;
		else if(obt==3&&hap==1)
			obs_lik[i]=(1-err_rate)*(1-err_rate);
		else if(obt==3&&hap==2)
			obs_lik[i]=err_rate*(1-err_rate);
		else
			cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		}
*/
}

double CHMM::BaumWelchCore(CObs **obs, long T, double *gamma, double **xi,
				 boolean doLog)
// Operations on a single observation sequence in the Baum-Welch loop are grouped here
{
  int	i, j, t;
  double logProb;
  double bjBeta;
  double *scale;
  double **alpha, **beta;

  alpha = SetMatrix(T, mN);// different size every time
  beta = SetMatrix(T, mN);

  scale = SetVector(T);
  
  logProb = Forward(alpha, scale, obs, T, doLog);
  Backward(beta, scale, obs, T);

  for (t = 1; t <= T - 1; t++){
    for (j = 1; j <= mN; j++) {
      gamma[j] = alpha[t][j] * beta[t][j];
      bjBeta =  beta[t+1][j] * mB->at(j, obs[t+1]);
      for (i = 1; i <= mN; i++){ 
	 xi[i][j] = alpha[t][i] *  mA->at(i, j) * bjBeta;
      }
    }//end for j
    Normalize(xi, mN, mN);
    Normalize(gamma, mN);
    
    mA->BWSum(xi);
    mB->BWSum(gamma, obs[t]);

    if(t==1){
      mPi->BWSum(gamma);
    }
  }// end for t

  // Step for t = T
  for (j = 1; j <= mN; j++) {
    gamma[j] = alpha[T][j] * beta[T][j];
  }
  Normalize(gamma, mN);
  mB->BWSum(gamma, obs[T]);

  delete [] alpha[1];
  delete [] alpha;
  delete [] beta[1];
  delete [] beta;
  delete [] scale;
  
  return logProb;
}

//===============================================================================

double CHMM::IterBaumWelch(CObsSeq *obsSeq, double *gamma, double **xi)
{
	  double deltaAB, deltaA, deltaB, deltaPi, delta;
	  double logProb;
	  int i;
	  const boolean NOLOG = FALSE;

	  mA->StartIter();// Zero sums used to cumulate results from each sequence
	  mB->StartIter();
	  mPi->StartIter();
  
	  for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:

	    logProb = BaumWelchCore(obsSeq->mObs[i], obsSeq->mNbObs[i], gamma, xi, NOLOG);

	  }// end loop over observation files

	  deltaA = mA->EndIter();
	  deltaB = mB->EndIter();
	  deltaPi = mPi->EndIter();
	  
	  deltaAB = deltaA > deltaB ? deltaA : deltaB;
	  delta = deltaAB > deltaPi ? deltaAB : deltaPi;

	  return delta;
}

//===============================================================================

void CHMM::RunFwdBwd(CObsSeq *obsSeq)
{
  double *scale;
  double **alpha, **beta;
  double logProb;
  int i;

    for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:
  long T = obsSeq->mNbObs[i];
  CObs **obs = obsSeq->mObs[i];
  alpha = SetMatrix(T, mN);// different size every time
  beta = SetMatrix(T, mN);
  scale = SetVector(T);

        // logProb = BaumWelchCore(obsSeq->mObs[i], obsSeq->mNbObs[i], gamma, xi, NOLOG);
        logProb = Forward(alpha, scale, obsSeq->mObs[i], obsSeq->mNbObs[i], FALSE);
        Backward(beta, scale, obs, T);
    }
}

void CHMM::LearnBaumWelch(CObsSeq *obsSeq)
// Follows (6.110) p. 369 of Rabiner-Huang
// The Baum-Welch loop is done over all the sequences
{
	int  iCount = 0;
	double delta;

	double *gamma = SetVector(mN);
	double **xi = SetMatrix(mN, mN);
	
	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

// Baum-Welch loop starts here
	do  {	
	  delta = IterBaumWelch(obsSeq, gamma, xi);

	  iCount++;
	  cout<< endl << "BW iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
	
	delete [] xi[1];
	delete [] xi;
	delete [] gamma;
}

//===============================================================================

double CHMM::SegmentalKMeansCore(CObs **obs, long T)
// Operations on a single observation sequence in the segmental K-means loop
{
  int	t;
  int thisQ, nextQ;
  double logProb;
  int *q;

  q = SetIntVector(T);// best state sequence
  
  logProb = ViterbiLog(obs, T, q);
  
  mPi->SKMSum(q[1]);

  nextQ = q[1];
  for (t = 1; t <= T - 1; t++){
    thisQ = nextQ;
    nextQ = q[t+1];
    mA->SKMSum(thisQ, nextQ);
    mB->SKMSum(thisQ, obs[t]);
  }// end for t

  // Step for t = T
  mB->SKMSum(nextQ, obs[T]);

  delete [] q;  
  return logProb;
}

//===============================================================================

double CHMM::IterSegmentalKMeans(CObsSeq *obsSeq)
{
	  double deltaAB, deltaA, deltaB, deltaPi, delta;
	  double logProb;
	  int i;

	  mA->StartIter();// Zero sums used to cumulate results from each sequence
	  mB->StartIter();
	  mPi->StartIter();
  
	  for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:

	    logProb = SegmentalKMeansCore(obsSeq->mObs[i], obsSeq->mNbObs[i]);

	  }// end loop over observation files

	  deltaA = mA->EndIter();
	  deltaB = mB->EndIter();
	  deltaPi = mPi->EndIter();
	  
	  deltaAB = deltaA > deltaB ? deltaA : deltaB;
	  delta = deltaAB > deltaPi ? deltaAB : deltaPi;

	  return delta;
}

//===============================================================================

void CHMM::LearnSegmentalKMeans(CObsSeq *obsSeq)
// Follows (6.15.2) p. 383 of Rabiner-Huang
// The Segmental K-Means loop is done over all the sequences
{
	int  iCount = 0;
	double delta;

	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

// Baum-Welch loop starts here
	do  {	
	  delta = IterSegmentalKMeans(obsSeq);

	  iCount++;

	  cout<< endl << "SKM iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
}

//===============================================================================

void CHMM::LearnHybridSKM_BW(CObsSeq *obsSeq)
// Follows (6.15.2) p. 383 of Rabiner-Huang
// Combine Segmental K-Means and Baum-Welch
{
	int  iCount = 0;
	double delta;
	double *gamma = SetVector(mN);
	double **xi = SetMatrix(mN, mN);

	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

	for(int i=0;i<10;i++){
	  delta = IterBaumWelch(obsSeq, gamma, xi);
	  iCount++;
	  cout<< endl << "BW iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;
	}
// Baum-Welch loop starts here
	do  {	
	  delta = IterSegmentalKMeans(obsSeq);

	  iCount++;

	  cout<< endl << "SKM iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
}

//===============================================================================


CObsSeq* CHMM::GenerateSequences(long nbSequences, long nbObs, int seed)
{
  int i;
  int anyState = 1;
  CObs* obsType = mB->PickObservation(anyState);// to pass observation type
  CObsSeq* obsSeq = new CObsSeq(obsType, nbSequences, nbObs);

  MyInitRand(seed);

  for(i=1;i<=nbSequences;i++){
	   obsSeq->mObs[i] = GenerateObservations(nbObs);
  }
  return obsSeq;
}

//===============================================================================


CObs** CHMM::GenerateObservations(long T)
// Generate an observation sequence of length T using A, B and Pi
{
   int t = 1;
	int currentState;
	CObs** obs;

	obs = new CObs*[T+1];
	assert(obs != NULL);

        currentState = mPi->PickInitialState();
        mA->InitDuration(currentState);
        obs[1] = mB->PickObservation(currentState);
 
        for (t = 2; t <= T; t++) {
                currentState =  mA->PickNextState(currentState);
                obs[t] =  mB->PickObservation(currentState);
        }
	return obs;
}

//===============================================================================

void CHMM::PrintStatesAndExpectedObs(CObsSeq *obsSeq,
				     ostream& stateFile, ostream& bestObsFile)
// Print sequences of most probable states 
// and sequences of expected observations 
// corresponding to those states
{
	long i, j, t, T, nbSequences;
	double logProb;
 	int *q;// most probable state sequence
	CObs **stateToObsMap;// list of expected observations for each state
	CObs *expectedObs;

	obsSeq->PrintHeader(stateFile);
	obsSeq->PrintHeader(bestObsFile);
	
	stateToObsMap = mB->MapStateToObs();// expected observation for each state

	for(i=1;i<=obsSeq->GetNbSequences();i++){ // Loop over observation files:
		T = obsSeq->mNbObs[i];
		stateFile <<"T= "<< T << endl;
		bestObsFile <<"T= "<< T << endl;
		
  		q = SetIntVector(T);// best state sequence
		logProb = ViterbiLog(obsSeq->mObs[i], T, q);

		for (t=1; t <= T; t++){
			stateFile << q[t] << " ";
			
			expectedObs = stateToObsMap[q[t]];
			expectedObs->Print(bestObsFile);
		}
		stateFile << endl;
		bestObsFile << endl;
		delete [] q;
	}
#if 1
	for (j=1; j<= mN; j++){
		delete stateToObsMap[j];
	}
#endif
	delete [] stateToObsMap;
}

//===============================================================================

#if 0

CObsSeq* CHMM::ReadSequences(ifstream &inputSeqFile)
// Read observation file into a data structure to make training faster
{
  long i, j, T, nbSequences;
  char magicID[32];
  CObs* obs;
  CObsSeq* obsSeq = new CObsSeq;

  mB->ReadFileHeader(inputSeqFile);// P5 or P6

  inputSeqFile >> magicID;
  assert(strcmp(magicID, "nbSequences=")==0);
  inputSeqFile >> nbSequences;
  obsSeq->mNbSequences = nbSequences;
  obsSeq->mObsCount = 0;
  obsSeq->mObs = new CObs**[nbSequences+1];
  obsSeq->mNbObs = new long[nbSequences+1];

  for(i=1;i<=nbSequences;i++){
    cout <<"Sequence "<< i <<endl;
    inputSeqFile >> magicID;
    assert(strcmp(magicID, "T=")==0);

    inputSeqFile >> T;// nb of observations for each sequence
    obsSeq->mNbObs[i] =  T;
  	 obsSeq->mObsCount += T;
    obsSeq->mObs[i] = new CObs*[T+1];// This array is from 1 to T
    assert(obsSeq->mObs[i] != NULL);

    for(j=1; j <=T; j++){
      obs = mB->ReadObsFrom(inputSeqFile);// obs type depends on type of mB
      obsSeq->mObs[i][j] = obs;
    }
  }

  return obsSeq;
}

#endif

//===============================================================================

double CHMM::FindDistance(CObsSeq *obsSeq, ostream &outFile)
// Find product of probabilities for each sequence
{
  long i, T, nbSequences;

  double logProb;
  double sumProb = 0.0;
  double *scale;
  double **alpha, **beta;
  const boolean YESLOG = TRUE;
  long stepCount = 0;
  double distance;


  nbSequences = obsSeq->mNbSequences;
  //outFile << "nbSequences= " << nbSequences << endl;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    T = obsSeq->mNbObs[i];
    stepCount += T;

    alpha = SetMatrix(T, mN);// different size every time
    beta = SetMatrix(T, mN);

    scale = SetVector(T);
  
    logProb = Forward(alpha, scale, obsSeq->mObs[i], T, YESLOG);

    sumProb += logProb;

//    cout<<i<<". T = "<<T<< ", logProb = " <<logProb<<", -prob/step = "<<-logProb/T << endl;
   // outFile << -logProb/T << endl;
    cout << -logProb/T << endl;

    delete [] alpha[1];
    delete [] alpha;
    delete [] beta[1];
    delete [] beta;
    delete [] scale;
  }

  cout << endl << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
  cout << "Average Forward distance: " << distance << endl;
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average Forward prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================

double CHMM::FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile, vector<READ*> *reads_list)
// Find for each sequence log prob corresponding to best state segmentation
{
  long i, T, nbSequences;

  double logProb;
  double sumProb = 0.0;
  int *q;// most probable state sequence
  double *prob;
  long stepCount = 0;
  double distance;


  nbSequences = obsSeq->mNbSequences;
 // outFile << "nbSequences= " << nbSequences << endl;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    T = obsSeq->mNbObs[i];
    stepCount += T;
  
    q = SetIntVector(T);// best state sequence
    prob = SetVector(T);
    logProb = ViterbiLog(obsSeq->mObs[i], T, q, prob);
    sumProb += logProb;

    //outFile << -logProb/T << endl;
    cout << -logProb/T << endl;
    int j;
    for(j=1;j<=T;j++) {
	//cout << q[j] << endl;
	outFile << q[j] << ":" << prob[j] << ": ";
	int k_max = ((*reads_list)[j-1])->GetSnpCount();
	for(int k=1;k<=k_max;k++)
		 outFile << " (" << (*reads_list)[j-1]->GetSnp(k-1)->GetPos() << ", " << (*reads_list)[j-1]->GetAllele(k-1) << ")";
		 //outFile << "(" << (*reads_list)[j-1]->GetSnp(k-1)->GetRef() << "," << (*reads_list)[j-1]->GetSnp(k-1)->GetAlt() << ")" << (*reads_list)[j-1]->GetAllele(k-1) << " ";
	outFile << endl;
    }
    delete [] q;
  }

 // outFile << endl << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
 // outFile << "Nb of steps: " << stepCount << endl;
  cout << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
//  outFile << "Average Viterbi distance: " << distance << endl;
//  cout << "Average Viterbi distance: " << distance << endl;
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average Viterbi prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================

double CHMM::FindCrossEntropyDistance(CObsSeq *obsSeq, ostream &outFile)
// Find for each sequence log prob corresponding to best state segmentation
{
  long i, T, nbSequences;

  const double kWeight = 1.0;
  double logProb, weightedLogProb;
  double sumProb = 0.0;
  int *q;// most probable state sequence
  long stepCount = 0;
  double distance;
  double logQToPiProb;


  nbSequences = obsSeq->mNbSequences;
  outFile << "nbSequences= " << nbSequences << endl;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    T = obsSeq->mNbObs[i];
    stepCount += T;
  
    q = SetIntVector(T);// best state sequence
    logProb = ViterbiLog(obsSeq->mObs[i], T, q);
    logQToPiProb = FindQToPiProb(T, q);
    weightedLogProb = (kWeight*logProb + logQToPiProb)/(1.0+kWeight);
    sumProb += weightedLogProb;

    outFile <<-weightedLogProb/T << endl;
    cout <<-logProb/T << ", "<< -logQToPiProb/T << endl;
    delete [] q;
  }

  cout << endl << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
  cout << "Average cross-entropy distance: " << distance << endl;
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average cross-entropy prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================

double CHMM::FindQToPiProb(long T, int *q)
// Find probability of seeing a Pi distribution given that we have a Qi distribution
{
  int i, t;
  double *probQ;
  double piProb, qProb;
  double qToPiProb = 0.0;

  probQ = SetVector(mN);// prob distribution of observed Qs
  SetToZero(probQ, mN);

  for (t = 1; t <= T; t++){
    probQ[q[t]]++;
  }
  NonZeroNormalizeRow(probQ, mN);

  for(i=1;i<=mN;i++){
    piProb = mPi->at(i);
    qProb = probQ[i];
    qToPiProb += piProb * log(qProb/piProb);
  }
  qToPiProb *= T;// each number of observations is T * mPi->at(i)

  delete [] probQ;

  return qToPiProb;
}

//===============================================================================

void CHMM::Print(ostream &outFile)
{
	mA->Print(outFile);
	mB->Print(outFile);
	mPi->Print(outFile);
}

//===============================================================================

//===============================================================================
//===============================================================================

