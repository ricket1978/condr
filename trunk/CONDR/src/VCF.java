import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;


public class VCF
{

	public static void readData(ArrayList<Exon> exons, BufferedReader br)
	{
		int exonIndex = 0;
		for (Exon e : exons)
			e.SNPs = 0;

		String line = null;
		try
		{
			Exon e = exons.get(exonIndex);
			while( (line = br.readLine()) != null && exonIndex < exons.size())
			{
				if (line.startsWith("#"))
					continue;

				String[] fields = line.split("\t");

				int chr;
				try
				{
					chr = Integer.parseInt(fields[0].substring(3));
				}catch(NumberFormatException nfe)
				{
					System.out.println(e);
					System.out.println(line);
					chr = 0;
				}
				if (chr != e.chr)
					continue;
				int position = Integer.parseInt(fields[1]);

				while (position > e.posLeft)
				{
					//System.out.println(position + "\t" + e.posLeft);
					exonIndex ++;
					if (exonIndex >= exons.size())
						break;
					e = exons.get(exonIndex);
				}
				if (position > e.posLeft && position < e.posRight)
				{
					//within the exon
					// might need to change fields later
					String genotype = fields[2].split(":")[0];
					//System.out.println(position + "\t"+ genotype + "\t" + e);
					if (genotype.matches("0?1") || genotype.matches("1?0"))
						e.SNPs += 1;
				}

			}

		} catch (IOException e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}
	}

}
