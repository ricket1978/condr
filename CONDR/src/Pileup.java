import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class Pileup
{
	int position;
	int heterozygous;
	int homozygous;
	int coverage;

	public Pileup(String line)
	{
		String[] fields = line.split("\t");
		position = Integer.parseInt(fields[1]);
		coverage = Integer.parseInt(fields[4]);
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

	/*
	 * Reads the pileup and simultaneously updates the exon lists with appropriate values
	 */
	public static void readData(ArrayList<Exon> exons, String pileupFileName)
	{
		int exonIndex = 0;
		ArrayList<Pileup> Pileups = new ArrayList<Pileup>();

		String line = null;
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(pileupFileName));

			Exon e = exons.get(exonIndex);
			while( (line = br.readLine()) != null)
			{
				Pileup p = new Pileup(line);
				Pileups.add(p);

				while(!Pileups.isEmpty() && e.posLeft > Pileups.get(0).position)
					// remove until that point
					Pileups.remove(0);

				if (e.containsPosition(p.position))
				{
					if (p.heterozygous != 0 || p.homozygous != 0)
						e.SNPPositions.put(p.position, p.heterozygous);
					e.FPKM += p.coverage;
				}

				while (e.posRight < p.position)
				{
					// done with the current exon
					exonIndex ++;
					if (exonIndex >= exons.size())
						break;

					// get next exon. initialize with values already loaded from pileup
					e = exons.get(exonIndex);
					for(Pileup pileup : Pileups)
					{
						if ( e.containsPosition(pileup.position))
						{
							e.FPKM += pileup.coverage;
							if (pileup.heterozygous != 0 || pileup.homozygous != 0)
								e.SNPPositions.put(pileup.position, pileup.heterozygous);
						}
					}
				}
			}
			
			for(Exon ex : exons)
			{
				ex.SNPs = (double)ex.SNPPositions.size()/ex.length();
				ex.FPKM = (double)ex.FPKM/ex.length();
			}
		} catch (IOException e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	public String toString()
	{
		return (position + ", " + coverage + ", " + heterozygous);
	}
}
