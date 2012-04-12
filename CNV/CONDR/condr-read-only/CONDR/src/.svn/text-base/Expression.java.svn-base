
/*
 * Stores necessary information about RNA expression data for the short read fragments 
 * that have already been mapped to the reference genome
 * 
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Expression
{
	Integer chr;
	Integer transcriptStart;
	Integer transcriptEnd;
	Double rpkm;
	
	/* 
	 * Parses a line (of .expr file) into the components
	 */
	Expression(String dataLine)
	{
		String[] fields = dataLine.split("\t");
		try
		{
			this.chr  = Integer.parseInt(fields[2].substring(3));
		}catch(NumberFormatException e)
		{
			this.chr = 0; // for now, ignore X, Y, M chromosomes
		}
		this.transcriptStart = Integer.parseInt(fields[3]);
		this.transcriptEnd = Integer.parseInt(fields[4]);
		this.rpkm = Double.parseDouble(fields[5]);
	}

	/*
	 * Goes through input file (br) obtained from mapping,
	 * parses and stores the data in an array
	 */
	public static ArrayList<Expression> readExon(String expressionFileName, int chromosome)
	{
		ArrayList<Expression> exprs = new ArrayList<Expression>();
		String line = null;
		try
		{
			System.out.println("Expression File Name: " + expressionFileName);
			BufferedReader br = new BufferedReader(new FileReader(expressionFileName));

			while( (line = br.readLine()) != null)
			{
				Expression expr = new Expression(line);
				exprs.add(expr);
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process expression file");
			e.printStackTrace();
			System.exit(0);
		}

		return exprs;
	}

	// Pretty Printing
	public String toString()
	{
		return("chr" + this.chr + " (" + this.transcriptStart + ", " + this.transcriptEnd + ") " + this.rpkm);
	}

}
