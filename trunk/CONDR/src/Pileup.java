import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class Pileup
{
	int position;
	int heterozygous;
	int coverage;

	public Pileup(String line)
	{
		String[] fields = line.split("\t");
		position = Integer.parseInt(fields[1]);
		coverage = Integer.parseInt(fields[4]);
		if (fields[2].equalsIgnoreCase(fields[3]))
			heterozygous = 1;
		else
			heterozygous = 0;
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return data;
	}

	public static void readData(ArrayList<Exon> exons, BufferedReader br)
	{
		int exonIndex = 0;
		ArrayList<Pileup> data = new ArrayList<Pileup>();

		String line;
		try
		{
			Exon e = exons.get(exonIndex);
			while( (line = br.readLine()) != null)
			{
				Pileup p = new Pileup(line);
				data.add(p);

				int pileupIndex = 0;
				while(e.posLeft > data.get(pileupIndex).position)
					// remove until that point
					data.remove(pileupIndex);

				if (e.posRight < p.position)
				{
					// done with the current exon
					e.SNPs = e.SNPPositions.size()/e.length();
					exonIndex ++;

					// get next exon. initialize with values already loaded from pileup
					e = exons.get(exonIndex);
					for(Pileup pileup : data)
					{
						if ( e.posLeft < pileup.position && e.posRight > pileup.position)
						{
							if (pileup.heterozygous != 0)
								e.SNPPositions.put(pileup.position, pileup.heterozygous);
							e.FPKM += pileup.coverage;
						}
					}
				}
				else
				{
					if (p.heterozygous != 0)
						e.SNPPositions.put(p.position, p.heterozygous);
					e.FPKM += p.coverage;
				}
			}
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
