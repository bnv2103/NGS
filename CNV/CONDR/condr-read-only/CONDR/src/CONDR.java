import java.io.*;
import java.util.*;
import cern.jet.random.*;

/**
 * Main calling function of program
 * Gets data; parses it; computes most likely states
 * @author arthiramachandran
 *
 */
public class CONDR
{
	static ArrayList<Exon> Exons = new ArrayList<Exon>();
	static ArrayList<Exon> ExpectedValues = new ArrayList<Exon>();
	static ArrayList<Exon> StdDeviations = new ArrayList<Exon>(); // TODO reconsider a format for these so they're not arrays of exons

	static String exonFileName = "";
	static String outputFileName = "";
	static ArrayList<String> baselineExonFileNames = new ArrayList<String>();
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static String parameterFileName = "";
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;
	static double threshold = 0;

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();

		/*
		 * -c "16-16" -e "WSimulated4.region16.out.integerValues.exon" -t -b "WZ186.region16.out.integerValues.exon,WZ740.region16.out.integerValues.exon,WZ3313.region16.out.integerValues.exon,WZ561389.region16.out.integerValues.exon" -p "SimulationParameterFile"
		 * 
		 * -c "16-16" -e "WZ1034.region16.out.exon" -t -b "WZ186.region16.out.exon,WZ740.region16.out.exon,WZ3313.region16.out.exon,WZ561389.region16.out.exon"
		 * -c "16-16" -e "WSimulated3.region16.out.exon.partial" -t -b "WZ186.region16.out.exon.partial,WZ740.region16.out.exon.partial,WZ3313.region16.out.exon.partial,WZ561389.region16.out.exon.partial" -p "ParameterFile"
		 * 
		 * -c "16-16" -e "./bin/WSimulated10000_1.partialChr16.09202010.exon" -t -b "WZ186.region16.out.integerValues.exon,WZ740.region16.out.integerValues.exon,WZ3313.region16.out.integerValues.exon,WZ561389.region16.out.integerValues.exon" -p "./bin/SimulationParameterFile.09202010.10000"
		 */

		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal order
		 */

		/*
		File dir1 = new File (".");
		File dir2 = new File ("..");
		try {
			System.out.println ("Current dir : " + dir1.getCanonicalPath());
			System.out.println ("Parent  dir : " + dir2.getCanonicalPath());
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		*/

		try
		{
			int nSamples = baselineExonFileNames.size();
			System.out.println("Number of baseline: " + nSamples);
			if (nSamples < 10)
				for(String s : baselineExonFileNames)
					System.out.println(s);
			for (int chromosome : chromosomes )
			{
				double currentTime = 0, totalExonReadTime = 0;

				//System.out.println("Chromosome " + chromosome);
				System.out.println("Calculating expected values from given files....");
				ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, chromosome);
				StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, chromosome, ExpectedValues);
				Exon.sortExons(ExpectedValues);
				Exon.sortExons(StdDeviations);

				System.out.println("Reading exon with measurements file....");
				currentTime = System.currentTimeMillis();
				Exons = Exon.readAndStoreExonFile(exonFileName, chromosome);
				totalExonReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);
				double gammaK = Exon.getGenomeWideGammaParametersK(Exons);
				double gammaTheta = Exon.getGenomeWideGammaParametersTheta(Exons, gammaK);
				double gammaSNPsK = Exon.getGenomeWideGammaParametersSNPsK(Exons);;
				double gammaSNPsTheta = Exon.getGenomeWideGammaParametersSNPsTheta(Exons, gammaK);;
				System.out.println("K:" + gammaK + "\tTheta:"+gammaTheta);
				System.out.println("K:" + gammaSNPsK + "\tTheta:"+gammaSNPsTheta);
				Gamma priorGamma = new Gamma(gammaK, gammaTheta, null);
				
				// TODO add chromosome checks
				System.out.println("Calculating States....");
				currentTime = System.currentTimeMillis();
				HiddenMarkovModel.getStates(Exons, ExpectedValues, StdDeviations, parameterFileName, threshold, nSamples, 
						priorGamma, gammaK, gammaTheta, gammaSNPsK, gammaSNPsTheta);
				System.out.println("Finished running HMM");
				double totalStateCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

				// prints the timing metrics to std out
				if (printTimingMetrics)
				{
					double endTime = System.currentTimeMillis();
					totalTime = (endTime - startTime)/1000F;
					System.out.println("Total Time: " + totalTime);
					System.out.println("Time for reading exons file       : " + totalExonReadTime + ", " + Exons.size());
					System.out.println("Time for calculating States       : " + totalStateCalcTime);
				}

				// Print confusion matrix
				ArrayList<Exon> groundTruth = new ArrayList<Exon>();
				String line = null; 
				try
				{
					BufferedReader br = new BufferedReader(new FileReader(exonFileName));
					while( (line = br.readLine()) != null)
					{
						Exon exon = new Exon(line, true);
						exon.state.stateName = Integer.parseInt(line.split("\t")[6]);
						if (exon.state.stateName > 6)
							exon.state.stateName = -1;
						groundTruth.add(exon);
					}
				} catch (IOException e)
				{
					System.err.println("Error: Unable to process exon file");
					e.printStackTrace();
					System.exit(0);
				}

				// Print output
				if (outputFileName.equals(""))
				{
					for(int i = 0; i<Exons.size(); i++)
					{
						System.out.println(Exons.get(i) + "\t" + ExpectedValues.get(i).SNPs + "\t" + ExpectedValues.get(i).FPKM + "\t" + groundTruth.get(i).state.stateName);
					}
				}
				else
				{
					Writer output = new BufferedWriter(new FileWriter(outputFileName));
					for(int i = 0; i<Exons.size(); i++)
						output.write(Exons.get(i) + "\t" + ExpectedValues.get(i).SNPs + "\t" + ExpectedValues.get(i).FPKM + "\t" + groundTruth.get(i).state.stateName + "\n");
					/*for(Exon e : Exons)
						output.write(e + "\n");
					 */
					output.close();
				}

				/*
				// calculate confusion matrix
				int[][] confusionMatrix = new int[HiddenMarkovModel.States.size()][HiddenMarkovModel.States.size()];
				for(int exonIndex = 0; exonIndex < Exons.size(); exonIndex++)
				{
					confusionMatrix[groundTruth.get(exonIndex).state.stateName][Exons.get(exonIndex).state.stateName]++;
				}
				System.out.println("RealValue\\CalcValue");
				System.out.println("\t0\t1\t2\t3\t4");
				System.out.println("-----------------------------------------------");
				for(int i=0; i<HiddenMarkovModel.States.size(); i++)
				{
					System.out.print(i + " |\t");
					for(int j=0; j<HiddenMarkovModel.States.size(); j++)
						System.out.print(confusionMatrix[i][j] + "\t");
					System.out.println();
				}

				// TP-FP matrix
				System.out.println("ROC values");
				System.out.println("State\tTrue_Pos\tFalse_Pos\tTotal_Ground\tTotal_Called");
				System.out.println("--------------------------------------------------------");
				for(int i=0; i<HiddenMarkovModel.States.size(); i++)
				{
					System.out.print("TP_"+i+" | \t");
					System.out.print(confusionMatrix[i][i]+ "\t");
					int totalReal = 0, false_pos = 0, total_bp_called = 0;
					for(int j=0; j<HiddenMarkovModel.States.size(); j++)
					{
						totalReal += confusionMatrix[i][j];
						total_bp_called += confusionMatrix[j][i]; 
						if(i==j)
							continue;
						else
							false_pos += confusionMatrix[j][i];
					}
					System.out.println(false_pos + "\t" + totalReal + '\t' + total_bp_called);
				}

				// # cnvs correctly called
				// called if >= 50% overlap between the ground truth and the call.
				HashMap<Integer, ArrayList<CNV>> groundTruthCNVs = new HashMap<Integer, ArrayList<CNV>>();
				// index == state
				// state: posLeft:posRight of CNV
				//for(State s: HiddenMarkovModel.States.values())
				State prevState = HiddenMarkovModel.States.get("NORMAL");
				Exon prevExon = null;
				for(Exon g : groundTruth)
				{
					if (prevState.stateName != g.state.stateName) // state transition
					{
						if(prevState.equals(HiddenMarkovModel.States.get("NORMAL")))
						{
							ArrayList<CNV> oldArray = groundTruthCNVs.get(g.state.stateName);
							if (oldArray==null)
								oldArray= new ArrayList<CNV>();
							CNV newCNV = new CNV(g.posLeft);
							oldArray.add(newCNV);
							groundTruthCNVs.put(g.state.stateName, oldArray);
						}
						else if(g.state.equals(HiddenMarkovModel.States.get("NORMAL")))
						{
							ArrayList<CNV> oldArray = groundTruthCNVs.get(prevState.stateName);
							CNV oldCNV = oldArray.get(oldArray.size()-1); // last element of old array
							if (oldCNV == null)
								System.out.println("?");
							if (oldCNV.posEnd !=0)
								System.out.println("Error in getting the CNVs " + oldCNV);
							oldCNV.posEnd = prevExon.posRight;
							oldArray.set(oldArray.size()-1, oldCNV);
							groundTruthCNVs.put(prevState.stateName, oldArray);
						}
						else
						{
							System.out.println("Weird states: " + prevState.stateName + "\t" + g.state.stateName + "\t" + g.posLeft +"--");
						}
					}
					prevState = g.state;
					prevExon = g;
				}

				HashMap<Integer, ArrayList<CNV>> calledCNVs = new HashMap<Integer, ArrayList<CNV>>();
				prevState = HiddenMarkovModel.States.get("NORMAL");
				prevExon = null;

				for(Exon g : Exons)
				{
					if (prevState.stateName != g.state.stateName) // state transition
					{
						if(prevState.equals(HiddenMarkovModel.States.get("NORMAL")))
						{
							ArrayList<CNV> oldArray = calledCNVs.get(g.state.stateName);
							if (oldArray==null)
								oldArray= new ArrayList<CNV>();
							CNV newCNV = new CNV(g.posLeft);
							newCNV.state = g.state.stateName;
							newCNV.exons ++;
							oldArray.add(newCNV);
							calledCNVs.put(g.state.stateName, oldArray);
						}
						else if(g.state.equals(HiddenMarkovModel.States.get("NORMAL")))
						{
							ArrayList<CNV> oldArray = calledCNVs.get(prevState.stateName);
							CNV oldCNV = oldArray.get(oldArray.size()-1); // last element of old array
							if (oldCNV == null)
								System.out.println("?");
							if (oldCNV.posEnd !=0)
								System.out.println("Error in getting the CNVs");
							oldCNV.posEnd = prevExon.posRight;
							oldArray.set(oldArray.size()-1, oldCNV);
							calledCNVs.put(prevState.stateName, oldArray);
						}
						else
						{
							System.out.println("Weird states: " + prevState.stateName + "\t" + g.state.stateName + "\t" + g.posLeft + "--");
						}
					}
					prevState = g.state;
					prevExon = g;
				}

				/*
				System.out.println("Ground Truth CNVs");
				for(int state : groundTruthCNVs.keySet())
					for(CNV cnv : groundTruthCNVs.get(state))
						System.out.println(state + "\t" + cnv.posStart + "\t" + cnv.posEnd);
				System.out.println("Called CNVs");
				for(int state : calledCNVs.keySet())
					for(CNV cnv : calledCNVs.get(state))
						System.out.println(state + "\t" + cnv.posStart + "\t" + cnv.posEnd);
				
				// compare how many have 50% overlap
				int[] correctlyCalled = new int[HiddenMarkovModel.States.size()];
				int[] totalCalled = new int[HiddenMarkovModel.States.size()];
				int[] totalGround = new int[HiddenMarkovModel.States.size()];
				int[] falsePositive = new int[HiddenMarkovModel.States.size()];
				for(int cnvState : calledCNVs.keySet())
				{
					//if (cnvState == 1)
						//System.out.println("@@");
					ArrayList<CNV> called = calledCNVs.get(cnvState);
					//for(int groundState : groundTruthCNVs.keySet())
					ArrayList<CNV> ground = groundTruthCNVs.get(cnvState); // only interested in those which are the same calls
					totalGround[cnvState] = ground.size();
					for(CNV c : called)
					{
						boolean noCall=true;
						//if (c.posStart==26384129)
							//System.out.println("^^");
						for(CNV g : ground)
						{
							if (c.overlaps(g) > .5)
							{
								correctlyCalled[cnvState]++;
								//System.out.println(c + "\t" + g);
								noCall=false;
								ground.remove(g);
								break;
							}
						}
						if(noCall)
							falsePositive[cnvState]++;
						totalCalled[cnvState]++;
					}
					if(totalCalled[cnvState] != called.size())
						System.out.println("Sizes dont match " + totalCalled[cnvState] + "\t" + called.size());
					if(correctlyCalled[cnvState] > totalGround[cnvState] )
						System.out.println("More correct than ground truth");
					if(falsePositive[cnvState] > totalCalled[cnvState] )
						System.out.println("More FP than called");
					System.out.println("State_" + cnvState + "_numberCalled\t" + totalCalled[cnvState]);
					System.out.println("State_" + cnvState + "_numberGround\t" + totalGround[cnvState]);
					System.out.println("State_" + cnvState + "_falsePositive\t" + falsePositive[cnvState]);
					System.out.println("State_" + cnvState + "_numberCorrect\t" + correctlyCalled[cnvState]);
					//System.out.println("State_" + cnvState + "_numberCorrect\t" + correctlyCalled[cnvState]);
				}

				// computing in terms of base pairs sequenced that are called correctly
				System.out.println("In sequenced base pairs");
				// calculate confusion matrix
				confusionMatrix = new int[6][6];
				for(int exonIndex = 0; exonIndex < Exons.size(); exonIndex++)
				{
					confusionMatrix[groundTruth.get(exonIndex).state.stateName][Exons.get(exonIndex).state.stateName]+=Exons.get(exonIndex).length();
				}
				System.out.println("RealValue\\CalcValue");
				System.out.println("\t0\t1\t2\t3\t4");
				System.out.println("-----------------------------------------------");
				for(int i=0; i<HiddenMarkovModel.States.size(); i++)
				{
					System.out.print("BP_" + i + " |\t");
					for(int j=0; j<HiddenMarkovModel.States.size(); j++)
						System.out.print(confusionMatrix[i][j] + "\t");
					System.out.println();
				}

				// TP-FP matrix
				System.out.println("ROC values");
				System.out.println("State\tTrue_Pos\tFalse_Pos\tTotal_Ground\tTotal_Called");
				System.out.println("----------------------------------------------------------");
				for(int i=0; i<HiddenMarkovModel.States.size(); i++)
				{
					System.out.print("BP_TP_"+i+" | \t");
					System.out.print(confusionMatrix[i][i]+ "\t");
					int totalReal = 0, false_pos = 0, total_bp_called = 0;
					for(int j=0; j<HiddenMarkovModel.States.size(); j++)
					{
						totalReal += confusionMatrix[i][j];
						total_bp_called += confusionMatrix[j][i]; 
						if(i==j)
							continue;
						else
							false_pos += confusionMatrix[j][i];
					}
					System.out.println(false_pos + "\t" + totalReal + '\t' + total_bp_called);
				}
				
				System.out.println("Threshold: " + threshold);
				
				/*
				System.out.println("Num_correct; Num_exons; prevExon state; prevExon.posRight; State");
				int num_exons = 0;
				int num_tp = 0;
				prevExon=groundTruth.get(0);
				for(int exonIndex = 1; exonIndex < Exons.size(); exonIndex++)
				{
					Exon e=Exons.get(exonIndex);
					Exon g=groundTruth.get(exonIndex);
					num_exons++;
					if (e.state.stateName == g.state.stateName)
						num_tp++;
					if (prevExon.state.stateName != g.state.stateName)
					{
						//num_cnvs++;
						System.out.println("*\t" + num_tp + "\t" + num_exons + "\t" + prevExon.state.stateName + "\t" + prevExon.posRight + "\t" + g.state.stateName);
						num_exons=0;
						num_tp=0;
					}
					prevExon=g;
				}
				*/
			}
		} catch (Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		} 
	}

	/**
	 * Parses arguments to program
	 * @param arguments list of arguments
	 */
	private static void parseArguments(String arguments[])
	{
		/*
		 * format:
		 * -e <exonFileName>
		 * -ex <expressions file>
		 * -sam <mapped reads, sam file>
		 * -ref <reference directory>
		 * -pileup <pileup file>
		 * -c <chromosomes (start-end)>
		 * -t (prints timing metrics)
		 * -o <output file name>
		 * -b <comma separated baseline file names>
		 */
		try {
			for(int index = 0; index < arguments.length; index ++)
			{
				String arg = arguments[index];
				if (arg.equals("-e"))
					exonFileName = arguments[index + 1];
				else if (arg.equals("-t"))
					printTimingMetrics = true;
				else if (arg.equals("-o"))
					outputFileName = arguments[index + 1];
				else if (arg.equals("-c"))
				{
					String[] fields = arguments[index+1].split("-");
					for (int c = Integer.parseInt(fields[0]); c <= Integer.parseInt(fields[1]); c++)
					{
						chromosomes.add(c);
					}
				}
				else if (arg.equals("-b"))
				{
					String[] fields = arguments[index+1].split(",");
					for (String f : fields)
					{
						baselineExonFileNames.add(f.trim());
					}					
				}
				else if (arg.equals("-p"))
				{
					parameterFileName = arguments[index + 1].trim();
				}
				else if ( arg.equals("-threshold"))
				{
					threshold = Double.parseDouble(arguments[index+1].trim());
				}
			}

			if (exonFileName.equals("") 
					|| chromosomes.isEmpty() || parameterFileName.equals("")) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java CONDR -e <exonFileName> " + 
						"-c <chromosomes (start-end)> " + "[-t] " + "[-o <output file name>]");
				System.exit(0);
			}
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

}

/**
 * Class to compute where the CNVs are. Mostly useful for verification of simulations
 * @author arthiramachandran
 *
 */
class CNV
{
	public CNV(int posStart)
	{
		this.posStart= posStart; 
	}
	public int length()
	{
		return(this.posEnd - this.posStart);
	}
	public double overlaps(CNV g)
	{
		// more stringent conditions
		int intersect = Math.min(this.posEnd, g.posEnd) - Math.max(this.posStart, g.posStart);
		//int union = Math.max(this.posEnd, g.posEnd) - Math.min(this.posStart, g.posStart);
		//return (double)intersect/union; 
		return (double)intersect/Math.min(g.length(), this.length());
	}
	int posStart=0;
	int posEnd=0;
	int state=-1;
	int exons=0;
	
	public String toString()
	{
		return(this.posStart + "-" + this.posEnd + " [" + this.state + "]");
	}
}
