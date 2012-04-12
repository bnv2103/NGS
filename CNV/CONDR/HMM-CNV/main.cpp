/* 
 * File:   main.cpp
 * Author: Andy Martin
 *
 * Created on September 28, 2011, 10:26 
 * 
 * Input:
 * CLA 1 - Genome to analyze (.bam)
 * CLA 2 - List of control genomes (.txt, containing list of .bam's)
 * CLA 3 - Reference sequence (.fasta)
 * CLA 4 - Window size (integer)
 * CLA 5 - Output filename
 * 
 */

#include "main.h"

const int MAX_CHR = 50;
int WINDOW_SIZE = 500;
int MAX_DEPTH = WINDOW_SIZE;
int DIFF_WINDOW_SIZE;
vector<int> genomeDepth[MAX_CHR];
vector< vector<int> > controlDepth[MAX_CHR];
vector<int> states; 
vector<int> observations;
vector< map <int,map <int, double> > > emission_probability[MAX_CHR];
map<string,int> chrMapping;
map<int,double> start_probability;
map<int,map<int, double> > transition_probability;
ofstream oFile;

int main(int argc, char** argv){
    try{
        if (argc < 5) throw "Not enough arguments (5 req)";
		setVariables();
		srand( time(NULL) );     //seed clock for random number generation
		BamReader reader;
		
        string bam_file = argv[1]; string bam_index_file = bam_file+".bai";
        string control_file_list = argv[2];
        string fasta_file = argv[3]; string fasta_index_file = fasta_file+".fai";
	DIFF_WINDOW_SIZE = atoi(argv[4])/WINDOW_SIZE;
	string output_file = argv[5]; oFile.open(output_file.c_str());
	if (!oFile.is_open()) throw "Unable to open output file";
		
        vector<unsigned int> chrLengths;
		
        ifstream is; string str; is.open(fasta_index_file.c_str(), ifstream::in);
        for (int i = 0; (is); ++i){
            getline(is, str);
            if (str != ""){
				int firstTab = str.find('\t',0);
				chrMapping[string(str.substr(0, firstTab))] = i;
				chrLengths.push_back(atoi(str.substr(firstTab, str.find('\t', firstTab+1)).c_str())/WINDOW_SIZE);
            }
            else break;
        }
        is.close();
        int nChromosomes = chrLengths.size();
		oFile << "finished loading " << nChromosomes << " chromosomes from FASTA file." << endl;
		cout << "finished loading " << nChromosomes << " chromosomes from FASTA file." << endl;
		
        int nControls = getLineCount(control_file_list);
        oFile << "observed " << nControls << " samples in the control file." << endl;
	cout << "observed " << nControls << " samples in the control file." << endl;
		
		
        int n = 0;
        oFile << "attempting to resize depth vectors..\n";
        cout << "attempting to resize depth vectors..\n";
        for (vector<unsigned int>::iterator it = chrLengths.begin(); it < chrLengths.end(); ++it){
            oFile << "  working on chr" << n << " (size == " << *it << ")\n";
            cout << "  working on chr" << n << " (size == " << *it << ")\n";
            genomeDepth[n].resize(*it, 0);
            emission_probability[n].resize(*it);
            controlDepth[n].resize(nControls);
            for (unsigned int j = 0; j < controlDepth[n].size(); ++j){ controlDepth[n][j].resize(*it,0); }
            ++n;
        }
        oFile << "finished!\n\n";
	cout << "finished!\n\n";
		
        oFile << "attempting to load primary genome (" << bam_file << ")\n";
	cout << "attempting to load primary genome (" << bam_file << ")\n";
		if (!reader.Open(bam_file)){
			throw "failed to open main bam file";
		}
		
		BamAlignment al;
		while (reader.GetNextAlignment(al)){
			if (al.IsFirstMate() && al.IsMapped()){
				genomeDepth[al.RefID][al.Position/WINDOW_SIZE] += 1;	
			}
		}
		reader.Close();
        
	is.open(control_file_list.c_str(), ifstream::in);
        for (int lineNumber = 0; (is); ++lineNumber){
            getline(is, str);
            if (str != ""){
                oFile << "attempting to load background genome (" << str << ")\n";
				cout << "attempting to load background genome (" << str << ")\n";
				reader.Open(str);
				while (reader.GetNextAlignment(al)){
					if (al.IsFirstMate() && al.IsMapped()){
						controlDepth[al.RefID][lineNumber][al.Position/WINDOW_SIZE] += 1;
					}
				}	
            }
            else break;
        }
		
		//oFile << "printing information--" << endl;
		for (int chrNum = 0; chrNum < nChromosomes; ++chrNum){
			for (unsigned int window = 0; window < chrLengths[chrNum]; ++window){
				if (genomeDepth[chrNum][window] >= MAX_DEPTH) genomeDepth[chrNum][window] = MAX_DEPTH - 1;
				cout << "gen: " << genomeDepth[0][window] << "\tcontrol: " << controlDepth[0][0][window] << endl;
			}
		}
		
        oFile << "attempting to generate & set emission probabilities..\n";
        cout << "attempting to generate & set emission probabilities..\n";
        int runningTally, runningVarianceTally, sampleMean;
		float sampleVariance, term2, alpha, beta;
        for (int chrNum = 0; chrNum < nChromosomes; ++chrNum){
            oFile << "  working on chr" << chrNum << " (size == " << controlDepth[chrNum][0].size() << ")\n";
            cout << "  working on chr" << chrNum << " (size == " << controlDepth[chrNum][0].size() << ")\n";
            for (unsigned int j = 0; j < controlDepth[chrNum][0].size(); ++j){
                runningTally = 0; runningVarianceTally = 0;
                for (unsigned int lineNum = 0; lineNum < controlDepth[chrNum].size(); ++lineNum) 
					runningTally += controlDepth[chrNum][lineNum][j];
                sampleMean = float(runningTally)/nControls;
                for (unsigned int lineNum = 0; lineNum < controlDepth[chrNum].size(); ++lineNum) 
					runningVarianceTally += pow((float(controlDepth[chrNum][lineNum][j]) - sampleMean), 2);
                sampleVariance = float(runningVarianceTally)/nControls;
                term2 = (sampleMean*(1-sampleMean)/sampleVariance) - 1;
				
                alpha = sampleMean*term2;
                beta = (1-sampleMean)*term2;
                double lambda, temp_lambda, emission_prob_tally = 0;
                if (sampleVariance == 0) lambda = sampleMean;
                else if (alpha <= 0) lambda = 0;
                else if (beta <= 0) lambda = 1;
                else {
                    beta_distribution<> dist(alpha, beta);
                    lambda = quantile(dist, rand()); 
                }
				
                // HMM Emission Probabilities:
                for (int k = -2; k < 3; ++k){
                    for (int value = 0; value < MAX_DEPTH; ++value){
                        //pmf = (e^-lambda)*(lambda^k)/(k!)
						//vector< map <int,map <int, double> > > emission_probability[MAX_CHR];
						
						if (k == -2) temp_lambda = 0.05*lambda;
						else if (k == -1) temp_lambda = 0.5*lambda;
						else if (k == 0) temp_lambda = lambda;
						else if (k == 1) temp_lambda = 1.5*lambda;
						else if (k == 2) temp_lambda = 2*lambda;
						
						//cout << chrNum << " " << j << " " << k << " " << value << " lambda-" << lambda << " " << retPowOverFact(temp_lambda, value);
						//emission_prob_tally += retPowOverFact(temp_lambda,value);
						//cout << "(new e_p_t == " << emission_prob_tally << endl;
						emission_probability[chrNum][j][k][value] = retPowOverFact(temp_lambda,value); 
						//emission_prob_tally = 0;
						
					}
                }
				
            }
        }
		
        oFile << "finished!\n\n\n";
        cout << "finished!\n\n\n";
		
		for (int chr = 0; chr < nChromosomes; ++chr){
			oFile << "attempting to run viterbi algorithm on chr" << chr << endl;
			cout << "attempting to run viterbi algorithm on chr" << chr << endl;
			forward_viterbi(genomeDepth[chr], states, start_probability, transition_probability, emission_probability[chr]);
		}
		
        return 1;
    }
    catch(char* ERROR){
        cout << string(ERROR) << endl;
        cout << " * Input:\
		\n * CLA 1 - Genome to analyze (.bam)\
		\n * CLA 2 - List of control genomes (.txt, containing list of .bam's)\
		\n * CLA 3 - Reference sequence (.fasta)\
		\n * CLA 4 - Window size (integer)\
		\n * CLA 5 - Output filename\n\n";
		return 0;
    }
}

//this method compute total probability for observation, most likely viterbi path 
//and probability of such path
void forward_viterbi(vector<int> &obs, vector<int> &states, map<int, double> &start_p, map<int, map<int, double> > &trans_p, vector<map<int, map<int, double> > > &emit_p) {
	map<int, Tracking> T;
	
	for(vector<int>::iterator state=states.begin(); state!=states.end();state++) {
		vector<int> v_pth;
		v_pth.push_back(*state);
		
		T[*state] = Tracking(start_p[*state], v_pth, start_p[*state]);
	}
	
	int obsNo = 0;
	for(vector<int>::iterator output=obs.begin(); output!=obs.end();output++) {
		map<int, Tracking> U;
		//go through all possible states..
		for(vector<int>::iterator next_state=states.begin(); next_state!=states.end(); next_state++) {
			Tracking next_tracker;
			//go through all possible states again..
			for(vector<int>::iterator source_state=states.begin(); source_state!=states.end(); source_state++) {
				Tracking source_tracker = T[*source_state];
				double p, prob = ( emit_p[obsNo/*chr position*/][*source_state][*output]*trans_p[*source_state][*next_state] );
				if (prob != 0) p = -log(prob);
				else p = 999999;
				
				//	cout << prob << " --> p = " << p << endl;
				source_tracker.prob += p;
				source_tracker.v_prob += p;
				
				// wrong? v
				next_tracker.prob += source_tracker.prob;
				
				//	cout << "comparing .. " << source_tracker.v_prob << " vs " << next_tracker.v_prob << endl;
				if(source_tracker.v_prob < next_tracker.v_prob) {
					next_tracker.v_path = source_tracker.v_path;
					next_tracker.v_path.push_back(*next_state);
					next_tracker.v_prob = source_tracker.v_prob;
				}
			}
			
			U[*next_state] = next_tracker;
		}
		
		T = U;
		obsNo += 1;
	}
	
	// apply sum/max to the final states
	Tracking final_tracker;
	
	for(vector<int>::iterator state=states.begin(); state!=states.end(); state++) {
		Tracking tracker = T[*state];
		final_tracker.prob += tracker.prob;
		
		if(tracker.v_prob < final_tracker.v_prob) {
			final_tracker.v_path = tracker.v_path;
			final_tracker.v_prob = tracker.v_prob;
		}
	}
	
	int count[6] = {0,0,0,0,0,0};
	int position = 0;
	for(vector<int>::iterator state=final_tracker.v_path.begin(); state!=final_tracker.v_path.end(); state++) {
		if (position % DIFF_WINDOW_SIZE == 0) cout << endl << "detected: ";

		cout << *state << " ";
		/*
		if (*state < 0) cout << "-2 ";
		else if (*state == 0) cout << "0 ";
		else cout << "2 ";
		*/
		count[*state + 2] += 1;
		position++;
	}
	cout << "\nStateArray == " << endl;
	for (int i = -2; i < 3; ++i) cout << i << " == " << count[i + 2] << endl;
	return;
}

double retPowOverFact(int lam, int val){
	// returns  exp(-lambda)*pow(lambda, value))/fact(value)
	double temp = 1;
	for (int n_lambda=lam; n_lambda; n_lambda--) temp /= 2.7182818284;
	while (val) {
		temp *= lam;
		temp /= val;
		val --;
	}
	return temp;
}

int getLineCount(string filename){
    FILE *f=fopen(filename.c_str(),"rb");
    int c=0,b;while ((b=fgetc(f))!=EOF) c+=(b==10)?1:0;fseek(f,0,SEEK_SET);
    return c;
}

void setVariables(){
	for (int i = 0; i < MAX_DEPTH; ++i) observations.push_back(i);
	
	states.push_back(-2); 
	states.push_back(-1); 
	states.push_back(0);
	states.push_back(1); 
	states.push_back(2);
	
	start_probability[-2] = 0.05; 
	start_probability[-1] = 0.05; 
	start_probability[0] = 0.8;
	start_probability[1] = 0.05; 
	start_probability[2] = 0.05;
	
	for (int i = -2; i < 3; ++i){
		for (int j = -2; j < 3; ++j){
			if (i == j) {
				cout << "tp[" << i << "][" << j << "] == 0.4" << endl;
                                transition_probability[i][j] = 0.4;
                        }
			if (j == 0) {
				cout << "tp[" << i << "][" << j << "] == 0.6" << endl;
				transition_probability[i][j] = 0.6;
			}
			else {
				cout << "tp[" << i << "][" << j << "] == 0.1" << endl;
				transition_probability[i][j] = 0.1;
			}
			
		}
	}
	
	return;
}


Tracking::Tracking() {
	prob = 999999999.0;
	v_prob = 99999999.0;
}

Tracking::Tracking(double p, vector<int> & v_pth, double v_p) {
	prob = p;
	v_path = v_pth;
	v_prob = v_p;
}
