/* ************************************************************************ *
 * ************************************************************************ *

   File: flexibleObsProb.C
   The class CFlexibleObsProb defines operations 
   for the observation distribution B
   when observations are flexibles with independent components

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   Date:  2-27-99 

 * ************************************************************************ *

   Modification Log:
   January 2003: Added a switch for eaither Gaussian or discrete pdfs of flexible components
                 Also, discrete components read their nbSymbols from an array

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
#include "obs.h"
#include "obsProb.h"
#include "discreteObsProb.h"
#include "gaussianObsProb.h"
#include "flexibleObsProb.h"

//===============================================================================

CFlexibleObsProb::CFlexibleObsProb(int *listNbSymbols, int nbStates, int nbComponents, int isGaussian)
// A list of nb of symbols is passed, defining how many are used for each component
{
	int i;
	
	mN = nbStates;
	mDimension = nbComponents;
	mComponentProb = new CObsProb*[mDimension+1];

	for(i=1;i<= mDimension; i++){
                if (isGaussian){
                    mComponentProb[i] = new CGaussianObsProb(listNbSymbols[i], nbStates);
                }
                else{ // discrete (histograms)
                    mComponentProb[i] = new CDiscreteObsProb(listNbSymbols[i], nbStates);
                }
	}
}

//===============================================================================

CFlexibleObsProb::~CFlexibleObsProb(void)
{
	int i;
	
	for(i=1;i<= mDimension; i++){
		delete mComponentProb[i];
	}
	delete [] mComponentProb;
}

//===============================================================================

void CFlexibleObsProb::SetVal(int state,int symbo,double value)
{
	cout << "SetVal not to be called from CFlexibleObsProb" << endl;
}

void CFlexibleObsProb::InitStateProb(double initProb)
{
	int i;
	double otherProb;

	otherProb = 1.0 - initProb;
	for(i=1;i<=mDimension;i++) {
		mComponentProb[i]->SetVal(1,1,initProb);
		mComponentProb[i]->SetVal(2,1,otherProb);
	}
}
/*
int CFlexibleObsProb::GetCount(void)
{
	int val;
	CIntObs *intObs = new CIntObs;
	double componentProb;
	
		val = (((CFlexibleObs<READ>*)this)->Get(1)); // component of obs
		if(!val) break;
		intObs->Set(val);
		componentProb = mComponentProb[1]->GetSnpCount();// prob of component
	delete intObs;
}
*/
double CFlexibleObsProb::at(int state, CObs *obs)
// return prob of seeing this obs flexible in this state
// equal to the product of the probs of seeing the components of the obs flexible in state
{
	int i;
	int val;
	CIntObs *intObs = new CIntObs;
	double prob = 1.0;
	double componentProb;
	
	for(i=1;i<= mDimension; i++){
		val = *(((CFlexibleObs<READ*>*)obs)->Get(i)); // component of obs
		if(!val) break;
		intObs->Set(val);
		componentProb = mComponentProb[i]->at(state, intObs);// prob of component
		prob *= componentProb;// multiply individual probabilities
	}

	delete intObs;
	return prob;
}

double CFlexibleObsProb::at(int state, int symbol)
// return prob of seeing this obs flexible in this state
// equal to the product of the probs of seeing the components of the obs flexible in state
{
	int i;
	double prob = 1.0;
	double componentProb;
	
	for(i=1;i<= mDimension; i++){
		componentProb = mComponentProb[i]->at(state, symbol);// prob of component
		prob *= componentProb;// multiply individual probabilities
	}

	return prob;
}

void CFlexibleObsProb::addEmission(int state, int symbol, double emission)
{
	int i;

	for(i=1;i<= mDimension; i++){
		mComponentProb[i]->addEmission(state, symbol, emission);// prob of component
	}
}

//===============================================================================

void CFlexibleObsProb::Start(void)
{
	int i;
	
	for(i=1;i<= mDimension; i++){
		mComponentProb[i]->Start();
	}
}

//===============================================================================

void CFlexibleObsProb::StartIter(void)
{
	int i;
	for(i=1;i<= mDimension; i++){
		mComponentProb[i]->StartIter();
	}
}

//===============================================================================

void CFlexibleObsProb::BWSum(double *gamma, CObs *obs)
{
	int i;
	int val;
	CIntObs *intObs = new CIntObs;
	
	for(i=1; i<=mDimension; i++){
		val = (int)*(((CFlexibleObs<READ*>*)obs)->Get(i)); // component i of obs flexible
		intObs->Set(val);
		mComponentProb[i]->BWSum(gamma, intObs);
	}
	delete intObs;
}

//===============================================================================

void CFlexibleObsProb::SKMSum(int state, CObs *obs)
{
	int i;
	int val;
	CIntObs *intObs = new CIntObs;
	
	for(i=1;i<= mDimension; i++){
		val = (int)*(((CFlexibleObs<READ*>*)obs)->Get(i)); // component of obs flexible
		intObs->Set(val);
		mComponentProb[i]->SKMSum(state, intObs);
	}
	delete intObs;
}

//===============================================================================

double CFlexibleObsProb::EndIter()
// Finish Baum-Welch iteration; return magnitude of largest change
{
	int i;
	double maxDiff = 0.0;
	double componentMaxDiff = 0.0;
	
	for(i=1;i<= mDimension; i++){
		componentMaxDiff = mComponentProb[i]->EndIter();
		if(componentMaxDiff > maxDiff) maxDiff = componentMaxDiff;
	}
   return maxDiff;
}

//===============================================================================

void CFlexibleObsProb::End()
// Finish Baum-Welch session
{
	int i;
	for(i=1;i<= mDimension; i++){
		mComponentProb[i]->End();
	}
}

//===============================================================================

CObs* CFlexibleObsProb::PickObservation(int state)
// Return a random observation given a state using the observation prob. matrix B
{
	int i;
	CIntObs *intObs;
	CFlexibleObs<READ*> *vectObs = new CFlexibleObs<READ*>(mDimension);

	for(i=1;i<= mDimension; i++){
		intObs = (CIntObs*)(mComponentProb[i]->PickObservation(state));
cout << "Alert: Entering PickObservation. Following line will fail" << endl;
		//vectObs->Set(intObs->Get(), i);
	}
	return vectObs;
}

//===============================================================================

CObs** CFlexibleObsProb::MapStateToObs(void)
// return array of expected observation flexibles for all states
{
        int i, j;
	CIntObs ***expectedComponentObs;
	CFlexibleObs<READ*> **expectedObs;
	CFlexibleObs<READ*> *flexibleObs;
	int obsVal;

	expectedComponentObs = new CIntObs**[mN+1];// array of maps 
	expectedObs = new CFlexibleObs<READ*>*[mN+1];

	for(i=1;i<= mDimension; i++){
	  expectedComponentObs[i] =  (CIntObs**)mComponentProb[i]->MapStateToObs();
	}  

	for (j=1; j <= mN; j++){// states
	  expectedObs[j] = new CFlexibleObs<READ*>(mDimension);
	  flexibleObs = expectedObs[j];
	  for(i=1;i<= mDimension; i++){
	         obsVal = expectedComponentObs[i][j]->Get();
cout << "Alert: Entering MapStateToObs. Following line will fail" << endl;
		 //flexibleObs->Set(obsVal, i);
	  }
	}
	delete [] expectedComponentObs[1];
	delete [] expectedComponentObs;

	return (CObs**)expectedObs;
}

//===============================================================================

void CFlexibleObsProb::Print(ostream &outFile)
{
	int i;

	for(i=1;i<= mDimension; i++){
		mComponentProb[i]->Print(outFile);
	}
}

//===============================================================================

CObs* CFlexibleObsProb::ReadObsFrom(ifstream &inFile)
{
	CFlexibleObs<READ*> *obs = new CFlexibleObs<READ*>(mDimension);
	obs->ReadFrom(inFile);

	return obs;
}

//===============================================================================

#if 0

void CFlexibleObsProb::ReadFileHeader(ifstream &inFile)
{
  CFlexibleObs::ReadHeader(inFile);
}

//===============================================================================

void CFlexibleObsProb::PrintFileHeader(ostream &outFile)
{
  CFlexibleObs::PrintHeader(outFile);
}

//===============================================================================
#endif
//===============================================================================
//===============================================================================
