
/*
 * This class defines the states included in the Markov Model
 * They are based on underlying biology and as such we define six states
 * 
 * NORMAL
 * HOMOZYGOUS_DELETE
 * HETEROZYGOUS_DELETE
 * COPY_NEUTRAL_LOH
 * INSERTION
 * CHROMATIN_CHANGE
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import cern.jet.random.Gamma;
import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;

/**
 * Holds the state name and the associated ratios expected <p>
 * Eg: In a het deletion, we expect [observed coverage]/[average coverage in other samples]~1/2 and 
 * [observed heterozygosity]/[average heterozygosity in other samples]~0 <p>
 * Also contains methods to compute the emission and transition probabilities.
 * 
 * @author arthiramachandran
 *
 */
public class State
{
	int stateName;
	double rpkmRatio;
	double snpRatio;
	int E_LengthOfState;

	static double rateOfOccurenceOfCNV = 0;

	static HashMap<State, HashMap<State, Double>> transitionProbabilities = new HashMap<State, HashMap<State, Double>>();

	// Associating state names with constants
	/* for RNA-seq
	 * public static final int NORMAL = 0, HOMOZYGOUS_DELETE = 1, HETEROZYGOUS_DELETE = 2, COPY_NEUTRAL_LOH = 3, 
	 * INSERTION = 4, DOUBLE_INSERT = 5; CHROMATIN_CHANGE = 5, UNKNOWN_TYPE = 6;
	 */
	public static final int NORMAL = 0, HOMOZYGOUS_DELETE = 1, HETEROZYGOUS_DELETE = 2, 
	INSERTION = 3, DOUBLE_INSERT = 4; 
	// Decides whether or not we compute the transition probability at run time or can pre-compute
	public static final int computeAtRunTime = -1;
	public static HashMap<Integer, Double> LOG_FACTORIAL = new HashMap<Integer, Double>();

	/**
	 * Gets the transition probability from state s1 -> state s2 <p>
	 * For CNV state -> normal state, transition probability is dependent on length of intron
	 * 
	 * @param s1 State transitioning from 
	 * @param s2 State transitioning to
	 * @param lengthOfIntron length of the intron
	 * @param lengthOfExon length of the exon
	 * @return transition probability from s1 -> s2
	 */
	public static double getTransitionProbability(State s1, State s2, int lengthOfIntron, int lengthOfExon)
	{
		// TODO see if this an appropriate fix
		double intronLength = 0;
		if ( lengthOfIntron < 0 )
			intronLength = 0.0000001;
		else
			intronLength = (double)lengthOfIntron;
		//System.out.println("Rate of occurence: " + rateOfOccurenceOfCNV );
		Double transitionProb = 0.0;
		if (s1 == null || s2 == null)
		{
			System.err.println("Invalid states in transition Probabilities");
			System.exit(0);
		}

		if (s1.equals(s2) )
		{
			if (s1.stateName!=NORMAL) // self loop in Markov Model
			{
				double lambda = (double)1.0 / s1.E_LengthOfState;
				transitionProb = Math.exp(-lambda * intronLength);
			}	
			else if (s1.stateName == NORMAL) // normal -> normal
			{
				// see below part
				// Poisson process with lambda = rate of occurrence of CNV
				double lambda = rateOfOccurenceOfCNV * lengthOfExon;
				transitionProb = lambda*Math.exp(-lambda)/5;
				lambda = rateOfOccurenceOfCNV;
				transitionProb = 1 - HiddenMarkovModel.States.size()*Math.exp(-lambda * intronLength);
				transitionProb = Math.exp(-lambda * intronLength);
				
				/*if (transitionProb < 0) // shouldnt happen
					transitionProb = 0.0;
					*/
			}
		}
		else if (! s1.equals(s2))
		{
			if (s2.stateName == NORMAL) // abnormal -> normal
			{
				double lambda = (double)1.0 / s1.E_LengthOfState;
				transitionProb = (1-Math.exp(-lambda * intronLength));
			}
			else if (s1.stateName == NORMAL) // normal -> abnormal
			{
				// Poisson process with lambda = rate of occurrence of CNV
				double lambda = rateOfOccurenceOfCNV * lengthOfExon;
				transitionProb = lambda*Math.exp(-lambda)/5;
				lambda = rateOfOccurenceOfCNV;
				transitionProb = (1-Math.exp(-lambda * intronLength)) / (HiddenMarkovModel.States.size()-1);
			}
			else transitionProb = 0.0;
		}
		else
			transitionProb = 0.0;

		return transitionProb;
	}

	/**
	 * Initializes all the states with the values in the parameter file
	 * 
	 * @param parameterFile name of parameter file
	 * @return HashMap of states indexed by their name
	 */
	/*
	 * Initializes all the state information: expected RPKM, SNP values and the allowed deviations from expectation
	 * This information is based on prior observations
	 * Contains Hard coded assumptions about the data
	 */
	public static HashMap<String, State> initializeStates(String parameterFile)
	{
		HashMap<String, State> states = new HashMap<String, State>();

		// Read parameter file and input states
		String line = "";
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(parameterFile));
			while( (line = br.readLine()) != null)
			{
				if (line.startsWith("//") || line.equals("")) // ignore comments
					continue;
				else if (line.startsWith("rate"))
				{
					rateOfOccurenceOfCNV = Double.parseDouble(line.split("=")[1]);
				}
				else
				{
					String[] lineElements = line.split("\\."); 
					String stateName = lineElements[0];
					State s = null;
					if (states.containsKey(stateName))
						s = states.get(stateName);
					else
					{
						s = new State();
						if (stateName.equals("NORMAL"))
							s.stateName = State.NORMAL;
						else if (stateName.equals("HOMOZYGOUS_DELETE"))
							s.stateName = State.HOMOZYGOUS_DELETE;
						else if (stateName.equals("HETEROZYGOUS_DELETE"))
							s.stateName = State.HETEROZYGOUS_DELETE;
						else if (stateName.equals("INSERTION"))
							s.stateName = State.INSERTION;
						else if (stateName.equals("DOUBLE_INSERT"))
							s.stateName = State.DOUBLE_INSERT;
						/* if RNA-seq
						 * else if (stateName.equals("COPY_NEUTRAL_LOH"))
						 * s.stateName = State.COPY_NEUTRAL_LOH;
						 * else if (stateName.equals("CHROMATIN_CHANGE"))
						 * s.stateName = State.CHROMATIN_CHANGE;
						 */
					}
					String fieldName = lineElements[1].split("=")[0].trim();
					String fieldValue = line.split("=")[1].trim();
					if (fieldName.equals("rpkmRatio"))
						s.rpkmRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("snpRatio"))
						s.snpRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("E_LengthOfState"))
						s.E_LengthOfState = Integer.parseInt(fieldValue);
					states.put(stateName, s);
				}			
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process parameter file");
			e.printStackTrace();
			System.exit(0);
		}

		return states;
	}

	/**
	 * For pretty-printing and for serving coffee with the code
	 * @return String character describing the name of the state
	 */
	public String getStateName()
	{
		switch(this.stateName)
		{
			case State.NORMAL: 				return "Normal"; 
			case State.HOMOZYGOUS_DELETE: 	return "Homozygous Deletion"; 
			case State.HETEROZYGOUS_DELETE: return "Heterozygous Deletion"; 
			case State.INSERTION: 			return "Insertion";
			case State.DOUBLE_INSERT: 	return "Double Insertion";
			//case State.COPY_NEUTRAL_LOH: 	return "Copy Neutral Loss of Heterozygosity"; 
			//case State.CHROMATIN_CHANGE: 	return "Chromatin Change"; 
			//case State.UNKNOWN_TYPE: 		return "Unknown Type"; 
		}
		return "";
	}

	// Printing functions
	public String toString()
	{
		String stateName = this.getStateName();
		return stateName + "\t" + this.rpkmRatio + "\t" + this.snpRatio + "\t";
	}

	// Defining our equality
	@Override public boolean equals(Object aThat) {
		if ( this == aThat ) return true;

		if ( !(aThat instanceof State) ) return false;
		State that = (State)aThat;

		return	(this.stateName == that.stateName);
	}

	// TODO: remove.. most likely useless function
	/**
	 * To remove
	 */
	public static double getEmissionProbability(double observedFPKM, double observedSNPs, double expectedFPKM,
			double expectedSNPs, double stdDevFPKM, double stdDevSNPs, double normalizationFPKM, double normalizationSNPs, State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming normal. Hence for other states, we multiple the mean
		// by the ratios for that state. We assume the same std dev for all the states distributions

		// Prob(value x | expected, stdDev, state) -> assuming this is a POISSON DISTRUIBUTION

		observedFPKM = observedFPKM/normalizationFPKM;
		expectedFPKM = expectedFPKM/normalizationFPKM;
		observedSNPs = observedSNPs/normalizationSNPs;
		expectedSNPs = expectedSNPs/normalizationSNPs;

		int lambda = (int)(s.rpkmRatio*expectedSNPs);
		double probSNPs = Math.exp((double)lambda) * Math.pow(lambda, (int)observedSNPs) / cern.jet.math.Arithmetic.factorial((int)observedSNPs);
		lambda = (int)(s.rpkmRatio*expectedFPKM);
		double probFPKM = Math.exp((double)lambda) * Math.pow(lambda, (int)observedFPKM) / cern.jet.math.Arithmetic.factorial((int)observedFPKM);

		// handling case where all the baselines are equal
		if (stdDevFPKM == 0.0)
			if (observedFPKM == s.rpkmRatio*expectedFPKM)
				probFPKM = 1;
			else 
				probFPKM = 0;
		if (stdDevSNPs == 0.0)
			if (observedSNPs == s.snpRatio*expectedSNPs)
				probSNPs = 1;
			else
				probSNPs = 0;

		// TODO: how to combine the two?
		return (probFPKM * probSNPs);
	}

	/**
	 * To remove
	 */
	// TODO: remove.. most likely useless function
	public static Double getEmissionProbability(Exon exon, Exon expected, Exon stdDev, State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming poisson. Hence for other states, we multiple the mean
		// by the ratios for that state. 

		double observedSNPs = exon.SNPs;
		double expectedSNPs = expected.SNPs; 
		double observedFPKM = exon.FPKM;
		double expectedFPKM = expected.FPKM;

		Poisson distFPKM = new Poisson (s.rpkmRatio*expectedFPKM, RandomEngine.makeDefault());
		double probFPKM = distFPKM.pdf((int)observedFPKM);
		Poisson distSNPs = new Poisson (s.rpkmRatio*expectedSNPs, RandomEngine.makeDefault());
		double probSNPs = distSNPs.pdf((int)observedSNPs);

		// Since the package handles 0! = 0 instead of 0! = 1, we have to do that manually		
		if (observedSNPs == 0.0)
		{
			double lambda = s.rpkmRatio*expectedSNPs;
			probSNPs = Math.exp(-lambda) * Math.pow(lambda, (int)observedSNPs) / cern.jet.math.Arithmetic.factorial((int)observedSNPs);
		}
		if (observedFPKM == 0.0)
		{
			double lambda = s.rpkmRatio*expectedFPKM;
			probFPKM = Math.exp(-lambda) * Math.pow(lambda, (int)observedFPKM) / cern.jet.math.Arithmetic.factorial((int)observedFPKM);
		}	


		// handling case where all the baselines are equal
		if (s.rpkmRatio*expectedFPKM == 0.0)
			if (observedFPKM == s.rpkmRatio*expectedFPKM)
				probFPKM = 1;
			else 
				probFPKM = 0;
		if (s.snpRatio*expectedSNPs == 0.0)
			if (observedSNPs == s.snpRatio*expectedSNPs)
				probSNPs = 1;
			else
				probSNPs = 0;


		if (Double.isInfinite(probSNPs) || Double.isNaN(probSNPs)  || probSNPs == 0.0)
			probSNPs = 10e-30;
		if (Double.isInfinite(probFPKM) || Double.isNaN(probFPKM) || probFPKM == 0.0)
			probFPKM = 10e-30;


		// TODO: how to combine the two?
		return (probFPKM * probSNPs);
	}

	/**
	 * Calculated the log of the factorial of the number
	 * @param number
	 * @return log(number!)
	 */
	public static double logFactorial(int number)
	{
		// storing some numbers for quicker computation
		Double value = LOG_FACTORIAL.get(number);
		if (value == null)
		{
			value = 0.0;
			for (int i=1; i<=number; i++)
			{
				value += Math.log(i);
			}
			LOG_FACTORIAL.put(number, value);
		}
		return value;
	}

	/**
	 * Calculates the emission probability for the observation at a particular exon, in state s
	 * 
	 * @param exon exon for which emission probability is computed
	 * @param expected expected observations for that exon (based on other samples)
	 * @param stdDev std deviations for that exon (based on other samples)
	 * @param s state
	 * @param nSamples number of samples
	 * @param gamma to remove
	 * @param gammaK gamma parameter k of the coverage
	 * @param gammaTheta gamma parameter Theta of the coverage
	 * @param gammaSNPsK gamma parameter k of the heterozygosity
	 * @param gammaSNPsTheta gamma parameter Theta of the heterozygosity
	 * @return log(emission probability(expected | s))
	 */
	public static double getLogEmissionProbability(Exon exon, Exon expected, Exon stdDev,
			State s, int nSamples, Gamma gamma, double gammaK, double gammaTheta, double gammaSNPsK, double gammaSNPsTheta)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming poisson. Hence for other states, we multiple the mean
		// by the ratios for that state. 

		double observedSNPs = exon.SNPs;
		double expectedSNPs = expected.SNPs; 
		double observedFPKM = exon.FPKM;
		double expectedFPKM = expected.FPKM;

		// using gamma distribution as a conjugate prior.
		// we start off with a gamma prior with arbitrary parameters k, theta
		// the posterior has parameters k+\sum{x_i} / (\theta/(n \theta + 1))
		// the expected value of the lambda posterior is (n * y_mean + a)/(n+b) 
		// from Bayesian modeling using WinBUGS By Ioannis Ntzoufras
		double a = .05, b = 0; // arbitrary numbers. Only there for worst case. Usually won't matter
		double lambdaSNPs = (s.snpRatio * nSamples * expectedSNPs + a * exon.length())/(nSamples + b);
		double logProbSNPs = logPoisson(observedSNPs, lambdaSNPs);
		a=100; b=2;
		double lambdaFPKM = (s.rpkmRatio * nSamples * expectedFPKM + a * exon.length())/(nSamples + b);
		double logProbFPKM = logPoisson(observedFPKM, lambdaFPKM);

		// handling case where all the baselines are equal, since log(0) = Inf 
		if (lambdaFPKM == 0.0)
			if (observedFPKM == 0.0)
				logProbFPKM = 0; // log 1 = 0; expected = observed = 0 so it fits with the expectation
			else
				logProbFPKM = Double.NEGATIVE_INFINITY;
		if (lambdaSNPs == 0.0)
			if (observedSNPs == 0.0)
				logProbSNPs = 0;
			else
				logProbSNPs = Double.NEGATIVE_INFINITY;

		/*
		 * Instead of gamma, using a negative binomial as the poisson parameter prior
		 * The mean/std dev of the neg. bin is inferred from the genomewide average
		 * parameters of poisson: p ~ genomewide average: average coverage/number of sites * length of exon
		 * r calculated from p : k_dash * ((1-p)/p)
		 * Using gamma as a prior distribution based on the genome wide coverage
		 * Using the sample values to update the gamma (based on the conjugate prior formualtion)
		 * This results in a new gamma distribution with different parameters
		 */
		
		double newGammaK = updateGammaK(s.rpkmRatio, nSamples, gammaK, expectedFPKM);
		double newGammaTheta = updateGammaTheta(nSamples, gammaTheta);
		gamma = new Gamma(newGammaK, newGammaTheta, null);
		// over 10 quantiles:
		logProbFPKM = Double.NEGATIVE_INFINITY;
		jsc.distributions.Gamma gammaDist = new jsc.distributions.Gamma(newGammaK, newGammaTheta);
		double prevx = 0;
		for(double p = 0.1; p <= 0.9; p+=0.1)
		{
			double x = gammaDist.inverseCdf(p); //inverse cdf of gamma with p^{th} percentile;
			if (prevx==0)
				prevx = x;
			x = (x + prevx)/2;
			logProbFPKM = Probability.logSum(logProbFPKM, logPoisson(observedFPKM, x) + logGamma(x, newGammaK, newGammaTheta));
			prevx = x;
		}
		
		if (observedFPKM == 0.0)
			if (s.rpkmRatio*expectedFPKM == 0)
				logProbFPKM = 0; // log 1 = 0; expected = observed = 0 so it fits with the expectation

		/*
		 * no need to do the same for SNPs
		double newGammaSNPsK = updateGammaK(s.snpRatio, nSamples, gammaSNPsK, expectedSNPs);
		double newGammaSNPsTheta = updateGammaTheta(nSamples, gammaSNPsTheta);
		gamma = new Gamma(newGammaSNPsK, newGammaSNPsTheta, null);
		// over 10 quantiles:
		logProbSNPs = Double.NEGATIVE_INFINITY;
		jsc.distributions.Gamma gammaSNPsDist = new jsc.distributions.Gamma(newGammaSNPsK, newGammaSNPsTheta);
		for(double p = 0.05; p <= 0.9; p+=0.05)
		{
			double x = gammaSNPsDist.inverseCdf(p); //inverse cdf of gamma with p^{th} percentile;
			logProbSNPs = Probability.logSum(logProbSNPs, logPoisson(observedSNPs, x) + logGamma(x, newGammaSNPsK, newGammaSNPsTheta));
		}
		if (observedSNPs == 0.0)
			if (s.snpRatio*expectedSNPs == 0)
				logProbSNPs = 0; // log 1 = 0; expected = observed = 0 so it fits with the expectation
		 */
		
		return (logProbFPKM + logProbSNPs);
	}

	/**
	 * Calculates the log of the pdf of a gamma distribution
	 * 
	 * @param x - observation
	 * @param k - parameter of Gamma distribution
	 * @param theta - parameter of Gamma distribution
	 * @return log(pdf(Gamma(k, Theta)))
	 */
	private static double logGamma(double x, double k, double theta)
	{
		return (k-1) * Math.log(x) - (x / theta) 
		- k * Math.log(theta) - cern.jet.stat.Gamma.logGamma(k);
	}

	/**
	 * Calculates log of pdf of Poisson distribution
	 * 
	 * @param x - observed
	 * @param lambda - parameter of Poisson distribution
	 * @return log(pdf(Poisson(lambda)))
	 */
	private static double logPoisson(double x, double lambda)
	{
		return -lambda + ((int)x * Math.log(lambda)) - logFactorial((int)x);
	}

	/**
	 * Updates the theta for the Gamma prior (computes the posterior hyperparameter)
	 * @param nSamples number of samples
	 * @param gammaTheta original Theta value
	 * @return new Theta value
	 */
	private static double updateGammaTheta(int nSamples, double gammaTheta)
	{
		return gammaTheta / (nSamples * gammaTheta + 1);
	}

	/**
	 * Updates the k for the Gamma prior (computes the posterior hyperparameter)
	 * 
	 * @param ratio ratio of observed to expected in the particular state
	 * @param nSamples number of samples
	 * @param gammaK original k value
	 * @param expectedFPKM expected coverage
	 * @return new k value
	 */
	private static double updateGammaK(double ratio, int nSamples, double gammaK, double expectedFPKM)
	{
		return gammaK + ratio * expectedFPKM * nSamples;
	}

}
