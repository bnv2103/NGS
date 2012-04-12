import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;


public class TabulateResults
{
	public static void main(String args[])
	{
		File dir = new File ("./Data/Plagnol_Test_Data/SimulationResults");
		try
		{
			System.out.println(dir.getCanonicalPath());
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int num_correct[][] = new int[6][30]; 
		int num_exons[][] = new int[6][30]; 
		File[] files = dir.listFiles(new FilenameFilter() {
	           public boolean accept(File dir, String name) {
	                return (name.toLowerCase().endsWith(".results.summary") && name.contains("11012010_RNA2"));
	                }
	           }
	        );
		for(File f : files)
		{
			System.out.println(f.getName());
			try
			{
				BufferedReader br = new BufferedReader(new FileReader(f.getAbsolutePath()));
				String line = null;
				while((line = br.readLine()) != null)
				{
					if (line.startsWith("-\t"))
					{
						String[] fields = line.split("\t");
						int groundstate = Integer.parseInt(fields[3]);
						if (groundstate == 0)
							continue;
						int length = Integer.parseInt(fields[2]);
						if (groundstate==4 && length==18)
							System.out.println("$$");
						int calledState = Integer.parseInt(fields[4]);
						num_correct[groundstate][length]+= Integer.parseInt(fields[1]);
						//if (calledState+groundstate==9) // equate 4&5
							//num_correct[groundstate][length]+= length/2;
						
						num_exons[groundstate][length]+=length;
					}
				}
			} catch (FileNotFoundException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//System.out.println("HI");
		}
		
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<30; j++)
			{
				System.out.println(i+"\t" + j+"\t" + num_correct[i][j] + "\t" + num_exons[i][j]);
			}
		}
	}
}
