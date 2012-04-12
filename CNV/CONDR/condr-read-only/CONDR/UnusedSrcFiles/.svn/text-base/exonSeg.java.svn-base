import java.io.*;
import java.util.*;


public class exonSeg
{
	static ArrayList<Exon> Exons = new ArrayList<Exon>();

	static String exonFileName = "";
	static String expressionFileName = "";
	static String mappedReadsFileName = "";
	static String referenceDirectory = "";
	static String outputFileName = "";
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();


		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal order
		 */

		System.out.println("Preprocessing input files to find chromosome boundaries...");
		System.out.println("\tExons file..");
		ArrayList<Integer> exonFileChrBoundaries = getChromosomalBoundaries(exonFileName, 0);
		System.out.println("\tExpression file..");
		ArrayList<Integer> exprFileChrBoundaries = getChromosomalBoundaries(expressionFileName, 2);
		System.out.println("\tMapped reads file..");
		ArrayList<Integer> readsFileChrBoundaries = getChromosomalBoundaries(mappedReadsFileName, 2);


		try
		{
			BufferedReader inputExons = new BufferedReader(new FileReader(exonFileName));
			BufferedReader inputExpr = new BufferedReader(new FileReader(expressionFileName));
			BufferedReader inputSAMData = new BufferedReader(new FileReader(mappedReadsFileName));

			for (int chromosome : chromosomes )
			{
				int arrayPosition = chromosome - chromosomes.get(0) + 1;

				System.out.println("Chromosome " + chromosome);
				System.out.println("Reading exons file....");
				int numberOfLines = exonFileChrBoundaries.get(arrayPosition) - exonFileChrBoundaries.get(arrayPosition - 1);
				double currentTime = System.currentTimeMillis();
				Exons = Exon.readExon(inputExons, chromosome, numberOfLines);
				double totalExonReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);

				System.out.println("Reading expression file....");
				numberOfLines = exprFileChrBoundaries.get(arrayPosition) - exprFileChrBoundaries.get(arrayPosition - 1);
				ArrayList<Expression> Expressions = new ArrayList<Expression>();
				currentTime = System.currentTimeMillis();
				Expressions = Expression.readExon(inputExpr, chromosome, numberOfLines);
				double totalExprReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				int numberOfExpr = Expressions.size();

				System.out.println("Calculating FPKMs....");
				currentTime = System.currentTimeMillis();
				Exon.getFPKM(Expressions, Exons);
				double totalFPKMCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
				Expressions.removeAll(Expressions); // explicitly deleting to free up memory

				System.out.println("Reading mapped reads SAM file....");
				numberOfLines = readsFileChrBoundaries.get(arrayPosition) - readsFileChrBoundaries.get(arrayPosition - 1);
				ArrayList<MappedReads> mappedReads = new ArrayList<MappedReads>();
				currentTime = System.currentTimeMillis();
				mappedReads = MappedReads.readMappedReads(inputSAMData, chromosome, numberOfLines);
				double totalReadsReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				MappedReads.sort(mappedReads);
				int numberOfReads = mappedReads.size();

				System.out.println("Reading reference genome file....");
				String referenceFileName = referenceDirectory + "/chr" + chromosome + ".fa";
				currentTime = System.currentTimeMillis();
				RandomAccessFile inputReference = new RandomAccessFile(referenceFileName, "r");
				for(Exon e : Exons)
				{
					e.getReferenceSequence(inputReference);
				}
				double totalRefCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

				System.out.println("Calculating SNPs....");
				currentTime = System.currentTimeMillis();
				Exon.getSNPs(Exons, mappedReads); 
				double totalSNPsCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
				mappedReads.removeAll(mappedReads);

				System.out.println("Calculating States....");
				currentTime = System.currentTimeMillis();
				HiddenMarkovModel.getStates(Exons);
				double totalStateCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

				// Print output
				if (outputFileName.equals(""))
				{
					// print to stdout
					for(Exon e : Exons)
						System.out.println(e);
				}
				else
				{
					Writer output = new BufferedWriter(new FileWriter(outputFileName));
					for(Exon e : Exons)
						output.write(e + "\n");
					output.close();
				}

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
					System.out.println("Time for calculating States       : " + totalStateCalcTime);
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
		 * -c <chromosomes (start-end)>
		 * -t (prints timing metrics)
		 * -o <output file name>
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
			
			if (exonFileName.equals("") || expressionFileName.equals("") || mappedReadsFileName.equals("") || referenceDirectory.equals("") || chromosomes.isEmpty())
			{
				System.err.println("Improper Usage:");
				System.err.println("java exonSeg -e <exonFileName> " + "-ex <expressions file> " + 
						"-sam <mapped reads, sam file> " + "-ref <reference directory> " + 
						"-c <chromosomes (start-end)> " + "[-t] " + "[-o <output file name>]");
				throw(new Exception("Improper Usage"));
			}
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

	private static ArrayList<Integer> getChromosomalBoundaries(String fileName, int chrColumnNumber)
	{
		ArrayList<Integer> boundaries = new ArrayList<Integer>();
		String line = "";

		try
		{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			int prevChr = 0;
			int lineNumber = 1;

			while( (line = br.readLine()) != null)
			{
				int chr;
				try
				{
					chr = Integer.parseInt((line.split("\t")[chrColumnNumber].substring(3)));
				}catch(NumberFormatException e)
				{
					chr = 0; // for now, ignore X, Y, M chromosomes
				}

				if (chr != prevChr)
					boundaries.add(lineNumber-1);
				lineNumber ++;
				prevChr = chr;
			}
			boundaries.add(lineNumber);

		} catch (FileNotFoundException e)
		{
			System.err.println("File Not Found: " + e.getMessage());
			e.printStackTrace();
		} catch (IOException e)
		{
			System.err.println("IO Exception: " + e.getMessage());
			e.printStackTrace();
		}


		return boundaries;
	}

}
