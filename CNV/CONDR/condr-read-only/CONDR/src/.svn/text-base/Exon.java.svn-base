import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 * Stores and manipulates relevant exon information
 * @author arthiramachandran
 *
 */
public class Exon
{
	int chr;
	int posLeft;
	int posRight;
	String geneName;
	double FPKM;
	String referenceSequence;
	double SNPs;
	int numberOfOverlappingReads;
	ArrayList<MappedReads> reads;
	HashMap<Integer, Integer> SNPPositions;
	State state;
	public double numPileupPositions;
	public double likelihoodRatio;
	public double pValue;

	public static final int READLENGTH = 50;

	Exon()
	{
		chr = 0;
		state = null;
		this.FPKM = 0;
		this.SNPs = 0;
		this.reads = new ArrayList<MappedReads>();
		this.numberOfOverlappingReads = 0;
		this.numPileupPositions = 0;
		this.SNPPositions = new HashMap<Integer, Integer>();
		this.likelihoodRatio = 0;
	}

	/**
	 * parse a line of data (modified SAM format) and input appropriate fields
	 * @param dataLine
	 * @param fromExonFile
	 */
	Exon(String dataLine, boolean fromExonFile)
	{
		String[] fields = dataLine.split("\t");
		try
		{
			//this.chr  = Integer.parseInt(fields[0]);
			this.chr  = Integer.parseInt(fields[0].substring(3));
			//this.chr = fields[0].charAt(fields[0].length()-1);

		}catch(NumberFormatException e)
		{
			this.chr = 0; // for now, ignore X, Y, M chromosomes
		}catch(Exception e)
		{
			System.err.println("Error in reading exon file: "+ e.getMessage());
			System.err.println("Errored exon data: " + dataLine);
			e.printStackTrace();	
		}
		try
		{
			if (fromExonFile)
			{
				/*
			this.posLeft = Integer.parseInt(fields[1]);
			this.posRight = Integer.parseInt(fields[1]);
			this.geneName = fields[0];
			this.FPKM = Double.parseDouble(fields[4]);
			this.SNPs = Double.parseDouble(fields[3]);
				 */
				this.posLeft = Integer.parseInt(fields[1]);
				this.posRight = Integer.parseInt(fields[2]);
				this.geneName = fields[3];
				//double length = this.posRight - this.posLeft;
				this.FPKM = Double.parseDouble(fields[5]); //length;
				this.SNPs = Double.parseDouble(fields[4]); //length;
			}
			else
			{
				/*
			this.posLeft = Integer.parseInt(fields[5]);
			this.posRight = Integer.parseInt(fields[6]);
			this.geneName = fields[7];
				 */
				this.posLeft = Integer.parseInt(fields[1]);
				this.posRight = Integer.parseInt(fields[2]);
				//this.geneName = fields[0];
				this.FPKM = 0;
				this.SNPs = 0;
				this.reads = new ArrayList<MappedReads>();
				this.numberOfOverlappingReads = 0;
				this.SNPPositions = new HashMap<Integer, Integer>();
			}
		}catch(Exception e)
		{
			System.err.println("Error in reading exon file "+ e.getMessage());
			System.err.println("Errored exon data: " + dataLine);
			e.printStackTrace();
		}
		this.state = new State(); 
	}

	// Pretty printing
	public String toString()
	{
		return("chr" + this.chr + "\t" + this.posLeft + "\t" + this.posRight + "\t" + this.geneName + "\t" + 
				this.SNPs + "\t" + this.FPKM + "\t" + this.state.stateName + "\t" + this.likelihoodRatio + "\t" + this.pValue);
	}

	/**
	 * Parses the file (br)
	 * Reads file until the next chromosome
	 * Loads into an array of Exon objects
	 * @param exonFileName
	 * @param chromosome
	 * @return Array of exons
	 */
	public static ArrayList<Exon> readExon(String exonFileName, int chromosome)
	{
		ArrayList<Exon> exons = new ArrayList<Exon>();
		String line = null; 
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(exonFileName));
			while( (line = br.readLine()) != null)
			{
				Exon exon = new Exon(line, false);
				exons.add(exon);
			}
			br.close();
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process exon file");
			e.printStackTrace();
			System.exit(0);
		}

		/* do this only when reading the exon files. not the pileup */
		
		//normalize coverage by the genome-wide average
		/*
		double genomeWideAverageCoverage = 0;
		for(Exon e : exons)
			genomeWideAverageCoverage += e.FPKM;
		System.out.println("Genome wide Avg (from normalization): " + genomeWideAverageCoverage);
		genomeWideAverageCoverage = genomeWideAverageCoverage/exons.size();
		for(Exon e : exons)
			e.FPKM = e.FPKM/genomeWideAverageCoverage;
			*/
		return exons;
		
	}

	/**
	 * Given an array of exons and expression data, calculate the FPKM levels per exon <p>
	 * Expressions & Exons should be sorted in order of beginning position <p>
	 * Think its deprecated
	 * @param Expressions
	 * @param Exons
	 */
	static void getFPKM(ArrayList<Expression> Expressions, ArrayList<Exon> Exons)
	{
		double count = 0;
		int exprIndex = 0;

		// Helps to account for transcripts which are present in multiple exons
		int numberOfMovedExprPositions = 0;
		int numberOfMovedExonPositions = 0;

		for (int exonIndex = 0; exonIndex < Exons.size(); exonIndex++)
		{
			// until we are indexing the same chromosome, increment indices
			// get to the same chromosome in exons & expressions
			/*
			while(exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr != Expressions.get(exprIndex).chr)
			{
				if(Exons.get(exonIndex).chr > Expressions.get(exprIndex).chr)
				{
					exprIndex++;
				}
				if(Exons.get(exonIndex).chr < Expressions.get(exprIndex).chr)
				{
					exonIndex++;	
					numberOfMovedExonPositions++;
				}
			}
			 */

			// go to the first transcript that in the exon
			while( exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr == Expressions.get(exprIndex).chr 
					&& Exons.get(exonIndex).posLeft > Expressions.get(exprIndex).transcriptStart)
			{
				exprIndex++;
				numberOfMovedExprPositions++;
			}


			// for all the transcripts within the exon
			// add up the expression levels as long is its within the boundaries of the exon
			while ( exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr  == Expressions.get(exprIndex).chr 
					&& Exons.get(exonIndex).posLeft <= Expressions.get(exprIndex).transcriptStart
					&& Exons.get(exonIndex).posRight >= Expressions.get(exprIndex).transcriptEnd )
			{
				numberOfMovedExprPositions++;
				count += Expressions.get(exprIndex).rpkm;
				exprIndex++;
			}

			Exons.get(exonIndex).FPKM = count;

			exonIndex = exonIndex - numberOfMovedExonPositions;
			exprIndex = exprIndex - numberOfMovedExprPositions;
			numberOfMovedExprPositions=0;
			numberOfMovedExonPositions=0;
			count = 0;

		}
	}

	/**
	 * Reads reference file and computes the reference sequence for the exon <p>
	 * Deprecated
	 * @param inputReference
	 */
	public void getReferenceSequence(RandomAccessFile inputReference)
	{
		String refRead = "";
		try
		{
			inputReference.seek(this.posLeft - 1);
			for(int i=0; i<(this.posRight - this.posLeft + 1); i++)
			{
				byte b = inputReference.readByte();
				char c = (char)b;
				if (c == '\n')
					i--;
				else if (c == '>') // ignore the header lines
				{
					i--;
					inputReference.readLine();
				}
				else
					refRead += ((char)b);
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process reference genome file: " + e.getMessage());
			e.printStackTrace();
			System.exit(0);
		}

		this.referenceSequence = refRead;
	}

	/**
	 * Gets the number of differences (ie SNPs) between the two sequences <p>
	 * Deprecated
	 * @param s1
	 * @param s2
	 * @return
	 */
	public static int getNumberOfDifferences(String s1, String s2)
	{
		if (s1.length() > READLENGTH)
			System.out.println("Read length > expected read length of " + READLENGTH); 
		// check if the lengths are equal
		if (s1.length() != s2.length())
			return -1;

		int count = 0;
		for(int pos = 0; pos < s1.length(); pos++)
		{
			if (s1.charAt(pos) != s2.charAt(pos))
				count++;
		}
		return count;
	}

	/**
	 * Gets the number of SNPs per exon <p>
	 * Deprecated
	 * @param exons
	 * @param mappedReads
	 */
	public static void getSNPs(ArrayList<Exon> exons, ArrayList<MappedReads> mappedReads)
	{
		for (int exonIndex = 0; exonIndex < exons.size(); exonIndex++)
		{
			Exon exon = exons.get(exonIndex);
			if (mappedReads.isEmpty())
			{
				System.err.println("No mapped reads present");
				return;
			}

			int readsIndex = 0;
			MappedReads read = null;
			while(!mappedReads.isEmpty())
			{
				read = mappedReads.get(readsIndex);
				if (exon.posLeft > (read.startPosition + read.read.length()))
				{
					//System.out.println("Read removed: " + read.startPosition);
					mappedReads.remove(readsIndex);
				}
				else
					break;
			}

			while ( exon.overlaps(read) ) // mappedRead is within exon
			{
				calculateSNPPositions(exon, read);
				exon.numberOfOverlappingReads ++;
				readsIndex ++;
				if (readsIndex >= mappedReads.size())
					break;
				read = mappedReads.get(readsIndex);
				//System.out.println( "Read position: " + read.startPosition );
			}
			exon.SNPs = (double)exon.SNPPositions.keySet().size()/exon.length();
		}
	}

	int length()
	{
		return (this.posRight-this.posLeft+1);
	}

	/**
	 * Calculated the positions of the SNPs <p>
	 * Deprecated
	 * @param exon
	 * @param read
	 */
	private static void calculateSNPPositions(Exon exon, MappedReads read)
	{
		String referenceSeq = "";
		String readSeq = "";
		int exonStart = Math.max(exon.posLeft, read.startPosition) - exon.posLeft;
		int exonEnd = Math.min(exon.posRight, read.startPosition+read.read.length()) - exon.posLeft;
		int readStart = Math.max(exon.posLeft, read.startPosition) - read.startPosition;
		int readEnd = Math.min(exon.posRight, read.startPosition+read.read.length()) - read.startPosition;

		referenceSeq = exon.referenceSequence.toLowerCase().substring(exonStart, exonEnd);
		readSeq = read.read.toLowerCase().substring(readStart, readEnd);

		for(int pos = 0; pos < referenceSeq.length(); pos++)
		{
			if (referenceSeq.charAt(pos) != readSeq.charAt(pos))
				exon.SNPPositions.put(pos+Math.max(exon.posLeft, read.startPosition), 1);
		}

	}

	/**
	 * Does exon (this) overlap with mapped reads? <p>
	 * Deprecated
	 * @param mappedRead
	 * @return
	 */
	private boolean overlaps(MappedReads mappedRead)
	{
		if (this.chr == mappedRead.chr)
		{
			int overlapDistance = Math.min(this.posRight, mappedRead.startPosition+mappedRead.read.length())
			- Math.max(this.posLeft, mappedRead.startPosition);
			if (overlapDistance > 0)
				return true;
			else
				return false;

		}
		else
			return false;
	}

	/**
	 * Sort exons in order of their genomic position
	 * @param data list of exons
	 */
	static void sortExons(ArrayList<Exon> data)
	{
		final Comparator<Exon> order =
			new Comparator<Exon>() 
			{
			public int compare(Exon e1, Exon e2) {
				// chr diff & position diff
				int positionDiff = e1.chr - e2.chr;
				if (positionDiff == 0)
					positionDiff = e1.posLeft - e2.posLeft;
				if (positionDiff == 0)
					positionDiff = e1.posRight - e2.posRight;
				return positionDiff;
			}
			};
			Collections.sort(data, order);
	}

	/**
	 * Does this exon contain the specified genomic position?
	 * @param position
	 * @return true/false
	 */
	public boolean containsPosition(int position)
	{
		if (this.posLeft <= position && position <= this.posRight)
			return true;
		else
			return false;
	}

	/**
	 * Huh?
	 * @return
	 */
	public String inputData()
	{
		return(this.chr + "\t" + this.posLeft + "\t" + this.posRight + "\t" + this.geneName + "\t" + 
				this.SNPs + "\t" + this.FPKM);

	}

	/**
	 * Reads exon file; stores necessary info in array of exons
	 * @param exonFileName name of file with list of exons
	 * @param chromosome
	 * @return array of exons
	 */
	public static ArrayList<Exon> readAndStoreExonFile(String exonFileName, int chromosome)
	{
		ArrayList<Exon> exons = new ArrayList<Exon>();
		String line = null; 
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(exonFileName));
			while( (line = br.readLine()) != null)
			{
				Exon exon = new Exon(line, true);
				exons.add(exon);
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process exon file");
			e.printStackTrace();
			System.exit(0);
		}

		return exons;
	}

	/**
	 * Computes the expected value for each exon by average the values in the baseline files
	 * @param baselineExonFileNames array of names of the baseline files (the samples that to be compared against)
	 * @param chromosome chromosome
	 * @return list of exons with values of coverage/heterozygosity being the averages for the baselines
	 */
	public static ArrayList<Exon> calculateExpectedValues(ArrayList<String> baselineExonFileNames, int chromosome)
	{
		ArrayList<Exon> baselineExonValues = new ArrayList<Exon>();

		// for each file, 'add' the values to the baseline
		for(String fileName : baselineExonFileNames)
		{
			String line = null; 
			int exonIndex = 0;
			try
			{
				//System.out.println(fileName);
				BufferedReader br = new BufferedReader(new FileReader(fileName));
				while( (line = br.readLine()) != null)
				{
					// first file read
					Exon e = new Exon(line, true);
					if (fileName.equals(baselineExonFileNames.get(0)))
					{
						baselineExonValues.add(e);
					}
					else
					{
						Exon baselineE = baselineExonValues.get(exonIndex);
						baselineE.FPKM += e.FPKM;
						baselineE.SNPs += e.SNPs;
						baselineExonValues.set(exonIndex, baselineE);
					}
					exonIndex ++;
				}
			} catch (IOException e)
			{
				System.err.println("Error: Unable to process exon file");
				e.printStackTrace();
				System.exit(0);
			}
		}

		// average the values
		for(Exon e: baselineExonValues)
		{
			e.SNPs = e.SNPs/baselineExonFileNames.size();
			e.FPKM = e.FPKM/baselineExonFileNames.size();
		}

		return baselineExonValues;
	}

	/**
	 * Computes the std deviation values for each exon by average the values in the baseline files
	 * @param baselineExonFileNames array of names of the baseline files (the samples that to be compared against)
	 * @param chromosome chromosome
	 * @param expectedValues list of expected values per exon
	 * @return list of exons with values of coverage/heterozygosity being the std deviation for the baselines
	 */
	public static ArrayList<Exon> calculateStdDevValues(ArrayList<String> baselineExonFileNames, int chromosome, ArrayList<Exon> expectedValues)
	{
		ArrayList<Exon> baselineExonStdDevValues = new ArrayList<Exon>();

		// for each file, 'add' the values to the baseline
		for(String fileName : baselineExonFileNames)
		{
			String line = null; 
			int exonIndex = 0;
			try
			{
				BufferedReader br = new BufferedReader(new FileReader(fileName));
				while( (line = br.readLine()) != null)
				{
					// first file read
					Exon e = new Exon(line, true);					
					Exon expected = expectedValues.get(exonIndex);
					if (fileName.equals(baselineExonFileNames.get(0)))
					{
						e.SNPs = (e.SNPs - expected.SNPs)*(e.SNPs - expected.SNPs);
						e.FPKM = (e.FPKM - expected.FPKM)*(e.FPKM - expected.FPKM);
						baselineExonStdDevValues.add(e);						
					}
					else
					{
						Exon baselineE = baselineExonStdDevValues.get(exonIndex);
						baselineE.FPKM += (e.FPKM - expected.FPKM)*(e.FPKM - expected.FPKM);
						baselineE.SNPs += (e.SNPs - expected.SNPs)*(e.SNPs - expected.SNPs);
						baselineExonStdDevValues.set(exonIndex, baselineE);
					}
					exonIndex ++;
				}
			} catch (IOException e)
			{
				System.err.println("Error: Unable to process exon file");
				e.printStackTrace();
				System.exit(0);
			}
		}

		// average the values
		for(Exon e: baselineExonStdDevValues)
		{
			//TODO: remove System.out.println(e);
			e.SNPs = Math.sqrt(e.SNPs/(baselineExonFileNames.size()-1));
			e.FPKM = Math.sqrt(e.FPKM/(baselineExonFileNames.size()-1));
		}

		return baselineExonStdDevValues;
	}

	/**
	 * Deprecated
	 * @param baselineExonFileNames
	 * @param i
	 * @return
	 */
	public static ArrayList<Exon> calculateNormalizationFactor(	ArrayList<String> baselineExonFileNames, int i)
	{
		ArrayList<Exon> minimalValues = new ArrayList<Exon>();	

		// assumes that the values are multiples of some common factor
		// makes sense since they are all computed by dividing a whole number by the exon length
		for(String fileName : baselineExonFileNames)
		{
			String line = null; 
			int exonIndex = 0;
			try
			{
				BufferedReader br = new BufferedReader(new FileReader(fileName));
				while( (line = br.readLine()) != null)
				{
					// first file read
					Exon e = new Exon(line, true);
					if (fileName.equals(baselineExonFileNames.get(0)))
					{
						if (e.SNPs == 0.0)
							e.SNPs = 1;
						if (e.FPKM == 0.0)
							e.FPKM = 1;
						minimalValues.add(e);
					}
					else
					{
						Exon minimalValue = minimalValues.get(exonIndex);
						if (e.FPKM < minimalValue.FPKM && e.FPKM != 0.0)
							minimalValue.FPKM = e.FPKM;
						if (e.SNPs < minimalValue.SNPs && e.SNPs != 0.0)
							minimalValue.SNPs = e.SNPs;
						minimalValues.set(exonIndex, minimalValue);
					}
					exonIndex ++;
				}
			} catch (IOException e)
			{
				System.err.println("Error: Unable to process exon file");
				e.printStackTrace();
				System.exit(0);
			}
		}
		for(Exon e : minimalValues)
		{
			if (e.SNPs <= 0.0)
				e.SNPs = 1;
			if (e.FPKM <= 0.0)
				e.FPKM = 1;
		}
		return minimalValues;
	}

	/**
	 * Find the k parameter for the Gamma distribution describing the genomewide coverage distribution
	 * <a href="http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation">MLE for Gamma</a>
	 * @param exons list of observation values per exon
	 * @return value of k parameter
	 */
	public static double getGenomeWideGammaParametersK(ArrayList<Exon> exons)
	{
		int N = exons.size();
		double sumExonCoverage = 0, sumLnExonCoverage = 0;
		for (Exon e : exons)
		{
			sumExonCoverage += e.FPKM;
			if (e.FPKM != 0)
				sumLnExonCoverage += Math.log(e.FPKM);
		}
		System.out.println("Genomewide avg coverage:" + sumExonCoverage/N);
		double s = Math.log(sumExonCoverage/N) - 1/N * sumLnExonCoverage;
		double k = ( 3 - s + Math.sqrt( Math.pow(s-3, 2) + 24*s) ) / ( 12 * s );
		return k;
	}

	/**
	 * Find the Theta parameter for the Gamma distribution describing the genomewide coverage distribution
	 * <a href="http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation">MLE for Gamma</a>
	 * @param exons list of observation values per exon
	 * @param k gamma parameter k
	 * @return value of Theta parameter
	 */
	public static double getGenomeWideGammaParametersTheta(ArrayList<Exon> exons, double k)
	{
		double sumExonCoverage = 0;
		for (Exon e : exons)
			sumExonCoverage += e.FPKM;
		int N = exons.size();
		double theta = 1 / (k * N) * sumExonCoverage;
		return theta;
	}

	/**
	 * Find the k parameter for the Gamma distribution describing the genomewide heterozygosity distribution
	 * <a href="http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation">MLE for Gamma</a>
	 * @param exons list of observation values per exon
	 * @return value of k parameter
	 */
	public static double getGenomeWideGammaParametersSNPsK(ArrayList<Exon> exons)
	{
		int N = exons.size();
		double sumExonSNPs = 0, sumLnExonSNPs = 0;
		for (Exon e : exons)
		{
			sumExonSNPs += e.SNPs;
			if (e.SNPs != 0)
				sumLnExonSNPs += Math.log(e.SNPs);
		}
		System.out.println("Genomewide avg snps:" + sumExonSNPs/N);
		double s = Math.log(sumExonSNPs/N) - 1/N * sumLnExonSNPs;
		double k = ( 3 - s + Math.sqrt( Math.pow(s-3, 2) + 24*s) ) / ( 12 * s );
		return k;
	}

	/**
	 * Find the Theta parameter for the Gamma distribution describing the genomewide heterozygosity distribution
	 * <a href="http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation">MLE for Gamma</a>
	 * @param exons list of observation values per exon
	 * @param k gamma parameter k
	 * @return value of Theta parameter
	 */
	public static double getGenomeWideGammaParametersSNPsTheta(ArrayList<Exon> exons, double k)
	{
		double sumExonSNPs = 0;
		for (Exon e : exons)
			sumExonSNPs += e.SNPs;
		int N = exons.size();
		double theta = 1 / (k * N) * sumExonSNPs;
		return theta;
	}

}
