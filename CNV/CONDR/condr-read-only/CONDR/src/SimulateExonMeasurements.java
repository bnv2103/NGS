import java.util.ArrayList;

import cern.jet.random.engine.RandomEngine;
import cern.jet.random.*;

public class SimulateExonMeasurements
{
	/*
	 * -c "16-16" -b "WZ186.region16.out.exon.partial,WZ740.region16.out.exon.partial,WZ3313.region16.out.exon.partial,WZ561389.region16.out.exon.partial"
	 * 
	 */
	// W3 created with variable lengths with means based on the parameter file
	static ArrayList<String> baselineExonFileNames = new ArrayList<String>();
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static ArrayList<Exon> exons = new ArrayList<Exon>();

	public static void main(String args[])
	{
		parseArguments( args );
		HiddenMarkovModel.initialize( args[0] );

		ArrayList<Exon> ExpectedValues = new ArrayList<Exon>();
		//ArrayList<Exon> NormalizationFactor = new ArrayList<Exon>();
		ArrayList<Exon> StdDeviations = new ArrayList<Exon>(); 

		ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, 1);
		StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, 1, ExpectedValues);

		// assigning state values
		// TODO: switch the ordering and such while testing since some states are smaller than others as a result
		int index = 0;
		// continue until we've exhausted all the exons with some data
		while( true)
		{
			double randomNumber = Math.random();
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index == -1)
				break;
			if (index > ExpectedValues.size() - 10) break;
			if (randomNumber < 1.0/HiddenMarkovModel.States.size())
				index = SimulateRegion(ExpectedValues, StdDeviations, "HOMOZYGOUS_DELETE", index);
			else if  (randomNumber < 2.0/HiddenMarkovModel.States.size())
				index = SimulateRegion(ExpectedValues, StdDeviations, "HETEROZYGOUS_DELETE", index);
			//else if  (randomNumber < 3.0/HiddenMarkovModel.States.size())
				//index = SimulateRegion(ExpectedValues, StdDeviations, "COPY_NEUTRAL_LOH", index);
			else if (randomNumber < 3.0/HiddenMarkovModel.States.size())
				index = SimulateRegion(ExpectedValues, StdDeviations, "INSERTION", index);
			else
				index = SimulateRegion(ExpectedValues, StdDeviations, "DOUBLE_INSERT", index);
			if (index == -1)
				break;
			//else 
				//index = SimulateRegion(ExpectedValues, StdDeviations, "CHROMATIN_CHANGE", index);
		}

		for(int i=0; i<exons.size(); i++)
		{
			Exon e = exons.get(i);
			double origFPKM = ExpectedValues.get(i).FPKM;
			e.SNPs = getLikelyValue(ExpectedValues.get(i).SNPs*e.state.snpRatio, 
					StdDeviations.get(i).SNPs);
			e.FPKM = getLikelyValue(ExpectedValues.get(i).FPKM*e.state.rpkmRatio, 
					StdDeviations.get(i).FPKM);
			System.out.println(e + "\t||" + origFPKM);
		}

	}

	private static int SimulateRegion(ArrayList<Exon> ExpectedValues,
			ArrayList<Exon> StdDeviations, String stateName, int index)
	{
		// TODO: check about Poisson distribution for this
		if (index >= ExpectedValues.size())
			return -1;
		int regionStartPosition = ExpectedValues.get(index).posLeft;
		//System.out.println(HiddenMarkovModel.States.get(stateName));
		Poisson dist = new Poisson(HiddenMarkovModel.States.get(stateName).E_LengthOfState, RandomEngine.makeDefault());
		int regionLength = dist.nextInt();
		Exon e = ExpectedValues.get(index);
		int genomicLength = 0;
		//for (int i=0; i<10; i++)
		//for(int seqLength = 0; seqLength <= regionLength; seqLength += (e.posRight - e.posLeft))
		for(genomicLength = 0; genomicLength <= regionLength; genomicLength = (e.posRight - regionStartPosition))
		{
			//if (e.posLeft == 33866544)
				//System.out.println("##");
			if (e.SNPs==0.0 && e.FPKM==0.0)
			{
				// did this to avoid long stretches of no data where nothing would be called anyways
				//System.out.println("!!" + e.geneName);
				e.state = HiddenMarkovModel.States.get("NORMAL");
			}
			else
				e.state = HiddenMarkovModel.States.get(stateName);
			//System.out.println("#\t"+e);
			exons.add(e);
			index ++;
			if (index >= ExpectedValues.size())
			{
				break;
			}
			e = new Exon();
			e = ExpectedValues.get(index);
		}
		return index;
	}

	private static double getLikelyValue(double mean, double stddev)
	{
		// TODO: check calculation of lambda parameter
		// TODO: don't need stddev anymore?
		// pick a value from poisson distribution with those parameters
		Poisson dist = new cern.jet.random.Poisson(mean, RandomEngine.makeDefault());
		//Normal dist = new Normal(mean, stddev, RandomEngine.makeDefault());
		double value = Double.NEGATIVE_INFINITY;
		while (value < 0)
			value = dist.nextDouble();
		return value;
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
				if (arg.equals("-c"))
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
			}

			if (chromosomes.isEmpty()) 
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
