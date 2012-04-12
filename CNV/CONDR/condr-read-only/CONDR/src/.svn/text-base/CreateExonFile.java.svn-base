import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.io.Writer;
import java.util.ArrayList;

public class CreateExonFile
{
	// TODO: Convert to script maybe? Or leave as Java?
	/*
	 * File Format of (tab delimited)
	 * chromosome exonStart exonEnd geneName FPKM/coverage SNP levels/BAF
	 * 
	 * .expr, .sam / .pileup -> .exon
	 * 
	 */
	static ArrayList<Exon> Exons = new ArrayList<Exon>();

	static String exonFileName = "";
	static String expressionFileName = "";
	static String mappedReadsFileName = "";
	static String referenceDirectory = "";
	static String outputFileName = "";
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;

	/*
	 * arguments:
	 * -e "./exonList/exonParsedUniqChr16" -ref "./chromFa" -pileup "WZ1034.region16.out" -c "16-16" -o "WZ1034.region16.exon" -t 
	 */

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();


		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal order
		 */
		try
		{
			for (int chromosome : chromosomes )
			{
				double currentTime = 0, totalExonReadTime = 0, totalExprReadTime = 0, totalFPKMCalcTime = 0, totalReadsReadTime = 0, totalSNPsCalcTime = 0, totalRefCalcTime = 0;
				int numberOfExpr = 0, numberOfReads = 0;

				System.out.println("Chromosome " + chromosome);
				System.out.println("Reading exons file....");
				currentTime = System.currentTimeMillis();
				Exons = Exon.readExon(exonFileName, chromosome);
				totalExonReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);

				// TODO add chromosome checks
				if (usingPileup)
				{
					BufferedReader br = new BufferedReader(new FileReader(mappedReadsFileName));
					System.out.println("Reading expression file....");
					ArrayList<Expression> Expressions = new ArrayList<Expression>();
					currentTime = System.currentTimeMillis();
					Expressions = Expression.readExon(expressionFileName, chromosome);
					totalExprReadTime = (System.currentTimeMillis() - currentTime)/1000F;
					numberOfExpr = Expressions.size();

					System.out.println("Calculating FPKMs....");
					currentTime = System.currentTimeMillis();
					Exon.getFPKM(Expressions, Exons);
					totalFPKMCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
					Expressions.removeAll(Expressions); // explicitly deleting to free up memory

					System.out.println("Reading pileup file, calculating coverage and SNPs....");
					currentTime = System.currentTimeMillis();
					Pileup.readData(Exons, br);
					totalReadsReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				}
				else
				{
					System.out.println("Reading expression file....");
					ArrayList<Expression> Expressions = new ArrayList<Expression>();
					currentTime = System.currentTimeMillis();
					Expressions = Expression.readExon(expressionFileName, chromosome);
					totalExprReadTime = (System.currentTimeMillis() - currentTime)/1000F;
					numberOfExpr = Expressions.size();

					System.out.println("Calculating FPKMs....");
					currentTime = System.currentTimeMillis();
					Exon.getFPKM(Expressions, Exons);
					totalFPKMCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
					Expressions.removeAll(Expressions); // explicitly deleting to free up memory

					System.out.println("Reading mapped reads SAM file....");
					ArrayList<MappedReads> mappedReads = new ArrayList<MappedReads>();
					currentTime = System.currentTimeMillis();
					mappedReads = MappedReads.readMappedReads(mappedReadsFileName, chromosome);
					totalReadsReadTime = (System.currentTimeMillis() - currentTime)/1000F;
					MappedReads.sort(mappedReads);
					numberOfReads = mappedReads.size();

					System.out.println("Reading reference genome file....");
					String referenceFileName = referenceDirectory + "/chr" + chromosome + ".fa";
					currentTime = System.currentTimeMillis();
					RandomAccessFile inputReference = new RandomAccessFile(referenceFileName, "r");
					for(Exon e : Exons)
					{
						e.getReferenceSequence(inputReference);
					}
					totalRefCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

					System.out.println("Calculating SNPs....");
					currentTime = System.currentTimeMillis();
					Exon.getSNPs(Exons, mappedReads); 
					totalSNPsCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
					mappedReads.removeAll(mappedReads);
				}

				// Print output
				Writer output = new BufferedWriter(new FileWriter(outputFileName));
				for(Exon e : Exons)
				{
					//System.out.println(e);
					output.write(e + "\n");
				}
				output.close();

				// prints the timing metrics to std out
				if (printTimingMetrics)
				{
					double endTime = System.currentTimeMillis();
					totalTime = (endTime - startTime)/1000F;
					System.out.println("Total Time: " + totalTime);
					System.out.println("Time for reading exons file       : " + totalExonReadTime + ", " + Exons.size());
					System.out.println("Time for reading expression file  : " + totalExprReadTime + ", " + numberOfExpr);
					System.out.println("Time for reading mapped reads file: " + totalReadsReadTime + ", " + numberOfReads);
					System.out.println("Time for getting reference seq    : " + totalRefCalcTime);
					System.out.println("Time for calculating FPKM         : " + totalFPKMCalcTime);
					System.out.println("Time for calculating Num of SNPs  : " + totalSNPsCalcTime);
				}
			}
		} catch (Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		} 
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
				else if (arg.equals("-ex"))
					expressionFileName = arguments[index + 1];
				else if (arg.equals("-sam"))
					mappedReadsFileName = arguments[index + 1];
				else if (arg.equals("-pileup"))
				{
					usingPileup = true;
					mappedReadsFileName = arguments[index + 1];
				}
				else if (arg.equals("-ref"))
					referenceDirectory = arguments[index + 1];
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
			}

			if (exonFileName.equals("") 
					|| chromosomes.isEmpty()
					|| outputFileName.equals("")
					|| ((!usingPileup && (expressionFileName.equals("") || mappedReadsFileName.equals(""))) 
							&& (usingPileup && mappedReadsFileName.equals("")))) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java CreateExonFile -e <exonFileName> " + "-ex <expressions file> " + 
						"-sam <mapped reads, sam file> " + "-ref <reference directory> " + 
						"-c <chromosomes (start-end)> " + "-o <output file name>");
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
