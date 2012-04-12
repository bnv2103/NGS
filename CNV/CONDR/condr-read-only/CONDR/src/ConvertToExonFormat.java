import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;

/**
 * File Format of output (tab delimited)
 * chromosome exonStart exonEnd geneName FPKM/coverage SNP levels/BAF
 * 
 * .pileup -> .exon
 * @author arthiramachandran
 *
 */
public class ConvertToExonFormat
{
	// TODO: Convert to script maybe? Or leave as Java?
	static ArrayList<Exon> Exons = new ArrayList<Exon>();

	static String exonCaptureArrayFileName = "";
	static String mappedReadsFileName = "";
	static String outputFileName = "";
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;

	/*
	 * arguments:
	 * -e "./exonList/exonParsedUniqChr16" -pileup "WZ1034.region16.out" -c "16-16" -o "WZ1034.region16.exon" -t 
	 */

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();


		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal/genomic order
		 */
		try
		{
			BufferedReader br = null;
			if (mappedReadsFileName == "")
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(mappedReadsFileName));

			double currentTime = 0, totalExonReadTime = 0, totalReadsReadTime = 0;

			for (int chromosome : chromosomes )
			{
				System.out.println("Chromosome " + chromosome);
				System.out.println("Reading exons file....");
				currentTime = System.currentTimeMillis();
				String exonFileNameByChr = "";
				//if (chromosomes.size() == 1)
				//exonFileNameByChr = exonCaptureArrayFileName;
				//else
				/*
				 * need to be able to handle exons by separate files
				 * need to be able to handle exons and one full file
				 * any specification of names should take place outside of this program
				 */
				// TODO: hardcoding in here
				exonFileNameByChr = exonCaptureArrayFileName + ".chr" + chromosome;
				Exons = Exon.readExon(exonFileNameByChr, chromosome);
				totalExonReadTime += (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);

				// TODO add chromosome checks
				//if (usingPileup)
				{
					System.out.println("Reading pileup file, calculating coverage and SNPs....");
					currentTime = System.currentTimeMillis();
					Pileup.readData(Exons, br);
					totalReadsReadTime += (System.currentTimeMillis() - currentTime)/1000F;
				}

				// Print output
				Writer output = new BufferedWriter(new FileWriter(outputFileName+".chr"+chromosome));
				System.out.println(outputFileName+".chr"+chromosome);
				for(Exon e : Exons)
				{
					//System.out.println(e);
					output.write(e + "\t" + e.numPileupPositions + "\n");
				}
				output.close();
			}
			// prints the timing metrics to std out
			if (printTimingMetrics)
			{
				double endTime = System.currentTimeMillis();
				totalTime = (endTime - startTime)/1000F;
				System.out.println("Total Time: " + totalTime);
				System.out.println("Time for reading exons file       : " + totalExonReadTime + ", " + Exons.size());
				System.out.println("Time for reading mapped reads file: " + totalReadsReadTime);
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
					exonCaptureArrayFileName = arguments[index + 1];
				else if (arg.equals("-sam"))
					mappedReadsFileName = arguments[index + 1];
				else if (arg.equals("-pileup"))
				{
					usingPileup = true;
					mappedReadsFileName = arguments[index + 1];
				}
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
						//chromosomes.add(fields[0].charAt(0));

					}
				}
			}

			if (exonCaptureArrayFileName.equals("") 
					|| chromosomes.isEmpty()
					|| outputFileName.equals("")) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java ConvertToExonFormat -e <exonFileName> " +  
						"-pileup <pileup fileName> " + 
						"-c <chromosomes (start-end)> " + "-o <output file name>");
				System.out.println(exonCaptureArrayFileName);
				System.out.println(chromosomes);
				System.out.println(outputFileName);
				System.out.println(mappedReadsFileName);
				System.exit(0);
				//throw(new Exception("Improper Usage"));
			}
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

}
