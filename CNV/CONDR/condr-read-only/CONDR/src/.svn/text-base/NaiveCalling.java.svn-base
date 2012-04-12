import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

import cern.jet.random.Gamma;


public class NaiveCalling
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

	HiddenMarkovModel hmm;
	/*
	 * compute average/std dev of the reference files
	 */
	public static void main(String args[])
	{
		parseArguments( args );
		try
		{
			HiddenMarkovModel.initialize(parameterFileName);
			int nSamples = baselineExonFileNames.size();
			System.out.println("Number of baseline: " + nSamples);
			if (nSamples < 10)
				for(String s : baselineExonFileNames)
					System.out.println(s);
			for (int chromosome : chromosomes )
			{
				//System.out.println("Chromosome " + chromosome);
				System.out.println("Calculating expected values from given files....");
				ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, chromosome);
				StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, chromosome, ExpectedValues);
				Exon.sortExons(ExpectedValues);
				Exon.sortExons(StdDeviations);
				
				System.out.println("Reading exon with measurements file....");
				Exons = Exon.readAndStoreExonFile(exonFileName, chromosome);
				Exon.sortExons(Exons);
				
				double genomeWideAve = getGenomeWideAverage(Exons);
				double genomeWideStdDev = getGenomeWideStdDev(Exons, genomeWideAve);
				
				for(int i=0; i<Exons.size(); i++)
				{
					Exon e = ExpectedValues.get(i);
					Exon s = StdDeviations.get(i);
					e.FPKM = genomeWideAve;
					s.FPKM = genomeWideStdDev;
				}
				System.out.println("Calculating States....");
				getStates(Exons, ExpectedValues, StdDeviations);
				System.out.println("Finished calling");

			}
		} catch (Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		} 
	}


	private static double getGenomeWideStdDev(ArrayList<Exon> exons, double avg)
	{
		int total = 0;
		for (Exon e: exons)
		{
			total += (e.FPKM - avg)*(e.FPKM - avg);
		}
		return Math.sqrt(total/exons.size());
	}


	private static double getGenomeWideAverage(ArrayList<Exon> exons)
	{
		int total = 0;
		for (Exon e: exons)
		{
			total += e.FPKM;
		}
		return total/exons.size();
	}


	private static void getStates(ArrayList<Exon> exons, ArrayList<Exon> expectedValues,
			ArrayList<Exon> stdDeviations)
	{
		for(int i = 0; i < exons.size(); i++)
		{
			Exon e = exons.get(i);
			Exon mean = expectedValues.get(i);
			Exon stddev = stdDeviations.get(i);

			double normalPrior = 0.8;
			double p[] = new double[5];
			double p_normal  = Normal(e.FPKM, mean.FPKM, stddev.FPKM);
			p[0]=p_normal*normalPrior;
			double p_homodel = Normal(e.FPKM, 0, stddev.FPKM);
			p[1] = p_homodel*(1-normalPrior)/4;
			double p_hetdel  = Normal(e.FPKM, mean.FPKM * 0.5, stddev.FPKM);
			p[2] = p_hetdel*(1-normalPrior)/4;
			double p_amp1    = Normal(e.FPKM, mean.FPKM * 1.5, stddev.FPKM);
			p[3] = p_amp1*(1-normalPrior)/4;
			double p_amp2 = Normal(e.FPKM, mean.FPKM * 2, stddev.FPKM);
			p[4] = p_amp2*(1-normalPrior)/4;
			if (e.posLeft == 288236)
				System.out.println("f");

			//String state[] = {"Homo", "Het", "Normal", "Amp1", "Amp2"};
			double maxprob = Double.NEGATIVE_INFINITY;
			int maxstate = -1;
			for(int s = 0; s<5; s++)
			{
				//System.out.println(p[s] + "\t" + e.FPKM + " " + mean.FPKM + " " + stddev.FPKM);
				if (p[s] > maxprob)
				{
					maxstate = s;
					maxprob = p[s];
				}
			}
			if (p[0] >= maxprob)
			{
				maxstate = 0;
				maxprob = p[0];
			}
			//System.out.println(maxstate);
			e.state = HiddenMarkovModel.getStateFromIndex(maxstate);
		}

		if (outputFileName.equals(""))
		{
			for(int i = 0; i<Exons.size(); i++)
			{
				Exon e = exons.get(i);
				Exon mean = expectedValues.get(i);
				Exon stddev = stdDeviations.get(i);
				System.out.println(e + "\t" + mean.FPKM + "\t" + stddev.FPKM);
			}
		}
		else
		{
			Writer output;
			try
			{
				output = new BufferedWriter(new FileWriter(outputFileName));
				for(int i = 0; i<Exons.size(); i++)
				{
					Exon e = exons.get(i);
					Exon mean = expectedValues.get(i);
					Exon stddev = stdDeviations.get(i);
					output.write(e + "\t" + mean.FPKM + "\t" + stddev.FPKM + "\n");
				}
				output.close();
			} catch (IOException e1)
			{
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}


	}

	private static double Normal(double x, double mean, double stddev)
	{
		
		double var = Math.pow(stddev, 2);
		double prob = Math.exp(-Math.pow(x-mean, 2)/(2*var));
		prob = prob / Math.sqrt(2 * Math.PI * var);
		if (stddev == 0.0)
			if (x==mean)
				prob = 1;
			else
				prob = 0;
		return prob;
		
		//double prob = (Math.pow(mean, x) * Math.exp(- mean)) / (Factorial(x));
		/*int m = (int)mean;
		double prob = (x * Math.log(mean)) - mean - State.logFactorial(m);
		return prob;
		*/
	}


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