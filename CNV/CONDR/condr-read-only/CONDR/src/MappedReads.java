
/*
 * Stores necessary information about RNA short read fragments that have already been 
 * mapped to the reference genome
 * 
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class MappedReads
{
	Integer chr;
	Integer startPosition;
	String read;
	Integer coverage;

	/* 
	 * Parses a line (in SAM format) into the components
	 */
	MappedReads(String line) 
	{
		String[] fields = line.split("\t");
		chr = Integer.parseInt(fields[2].substring(3));
		startPosition = Integer.parseInt(fields[3]);
		coverage = Integer.parseInt(fields[4]);
		read = fields[9];
		
	}

	// Pretty Printing
	public String toString()
	{
		return( this.chr + ": " + this.startPosition + ": " + this.coverage + ": " + this.read );
	}

	/*
	 * Goes through input file (br) which is in SAM format, obtained from mapping,
	 * parses and stores the data in an array
	 */
	public static ArrayList<MappedReads> readMappedReads(String mappedReadsFileName, int chromosome)
	{
		ArrayList<MappedReads> reads = new ArrayList<MappedReads>();
		String line = null;
		//int currentChromosome = chromosome;
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(mappedReadsFileName));
			while( (line = br.readLine()) != null)
			{
				if (line.startsWith("@")) // header
				{
					//lineCount ++;
					continue;
				}
				MappedReads read = new MappedReads(line);
				//currentChromosome = read.chr;
				reads.add(read);
			}

			/*
			while( currentChromosome == chromosome )
			{
				line = br.readLine();
				if (line == null)
					return reads; 
				MappedReads read = new MappedReads(line);
				currentChromosome = read.chr;
				reads.add(read);				
			}
			*/
			// TODO: at the switch of chromosomes, move the file pointer back			
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process short read file");
			e.printStackTrace();
			System.exit(0);
		}

		return reads;
	}

	public static void sort(ArrayList<MappedReads> data)
	{
		// TODO Auto-generated method stub
		
			final Comparator<MappedReads> order =
				new Comparator<MappedReads>() 
				{
				public int compare(MappedReads m1, MappedReads m2) {
					// chr diff & position diff
					int positionDiff = m1.chr - m2.chr;
					if (positionDiff == 0)
						positionDiff = m1.startPosition - m2.startPosition;
					return positionDiff;
				}
				};
				Collections.sort(data, order);
	}

}
