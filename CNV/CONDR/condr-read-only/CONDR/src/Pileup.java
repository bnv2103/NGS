
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class to handle pileup files
 * @author arthiramachandran
 *
 */
public class Pileup
{
	int position;
	int heterozygous;
	int homozygous;
	int coverage;
	int chromosome;

	public Pileup(String line)
	{
		String[] fields = line.split("\t");
		position = Integer.parseInt(fields[1]);
		coverage = Integer.parseInt(fields[7]);
		try{
			chromosome = Integer.parseInt(fields[0]);
			//chromosome = fields[0].charAt(fields[0].length()-1);
		} catch(NumberFormatException nfe)
		{
			try{
				chromosome = Integer.parseInt(fields[0].substring(3));
			}catch(NumberFormatException nfe2)
			{
				chromosome = 0;
			}
		}
		if (fields[2].equalsIgnoreCase(fields[3]))
			heterozygous = 0;
		else if (fields[3].equals("R") || fields[3].equals("Y") 
				|| fields[3].equals("K") || fields[3].equals("M") 
				|| fields[3].equals("S") || fields[3].equals("W"))
			heterozygous = 1;
		else if (fields[3].equals("A") || fields[3].equals("C")
				|| fields[3].equals("G") || fields[3].equals("T"))
			homozygous = 1;
	}

	/**
	 * Reads and parses the pileup file
	 * 
	 * @param br buffered reader to the pileup file
	 * @return hashmap of the pileups 
	 */
	public static HashMap<Integer, Pileup> readData(BufferedReader br)
	{
		HashMap<Integer, Pileup> data = new HashMap<Integer, Pileup>();
		String line;
		try
		{
			while( (line = br.readLine()) != null)
			{
				Pileup p = new Pileup(line);
				data.put(p.position, p);
			}
		} catch (IOException e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}
		return data;
	}

	/**
	 * Reads the pileup and simultaneously updates the exon lists with appropriate values
	 * 
	 * @param exons list of exons to be populated with values from the pileup
	 * @param br buffered reader to pileup file
	 */
	public static void readData(ArrayList<Exon> exons, BufferedReader br)
	{
		int exonIndex = 0;
		// TODO: move this up so pileups are stored between iterations
		// TODO: might be missing some pileups here 
		ArrayList<Pileup> Pileups = new ArrayList<Pileup>();

		String line = null;
		try
		{
			int count = 0;
			Exon e = exons.get(exonIndex);
			while( (line = br.readLine()) != null)
			{
				count ++;
				Pileup p = new Pileup(line);
				Pileups.add(p);
				if (count%1000000==0)	
					System.out.println(p + "\t" + Pileups.size());
				if (e.chr > p.chromosome) // reached end of chromosome
				{
					//System.out.println(e + "==" + p);
					Pileups.clear();
					continue;
				}
				//else if (e.chr < p.chromosome)
					//break;
				while(!Pileups.isEmpty() && e.posLeft > Pileups.get(0).position)
					// remove until that point
					Pileups.remove(0);

				if (e.containsPosition(p.position))
				{
					if (p.heterozygous != 0 || p.homozygous != 0)
						e.SNPPositions.put(p.position, p.heterozygous);
					e.FPKM += p.coverage;
					e.numPileupPositions ++;
				}

				while (e.posRight < p.position)
				{
					// done with the current exon
					exonIndex ++;
					if (exonIndex >= exons.size())
						break;

					// get next exon. initialize with values already loaded from pileup
					if (e.posLeft == 867769)
						System.out.println("***" + e + "\t" + e.numPileupPositions + "\t" + e.FPKM/e.numPileupPositions);
					e = exons.get(exonIndex);
					for(Pileup pileup : Pileups)
					{
						if ( e.containsPosition(pileup.position))
						{
							e.FPKM += pileup.coverage;
							e.numPileupPositions ++;
							if (pileup.heterozygous != 0 || pileup.homozygous != 0)
								e.SNPPositions.put(pileup.position, pileup.heterozygous);
						}
					}
				}
			}

			for(Exon ex : exons)
			{
				ex.SNPs = (double)ex.SNPPositions.size(); ///ex.length();
				if (ex.numPileupPositions == 0)
					ex.numPileupPositions = 1;
				ex.FPKM = (double)ex.FPKM/ex.numPileupPositions;
				//ex.FPKM = (double)ex.FPKM/ex.length();
			}
		} catch (IOException e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}
	}

	public String toString()
	{
		return ("chr" + chromosome + "\t" + position + ": cov=" + coverage + "\t het=" + heterozygous);
	}
}
