import java.io.*;
import java.util.*;


public class CONDR
{
	static ArrayList<Exon> Exons = new ArrayList<Exon>();
	static ArrayList<Exon> ExpectedValues = new ArrayList<Exon>();
	static ArrayList<Exon> StdDeviations = new ArrayList<Exon>(); // TODO reconsider a format for these so they're not arrays of exons
	// maybe extend exon to take n types of inputs?

	static String exonFileName = "";
	static String outputFileName = "";
	static ArrayList<String> baselineExonFileNames = new ArrayList<String>();
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();

		/*
		 * -c "16-16" -e "WZ1034.region16.out.exon" -t -b "WZ186.region16.out.exon,WZ740.region16.out.exon,WZ3313.region16.out.exon,WZ561389.region16.out.exon"
		 */

		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal order
		 */
		try
		{
			for (int chromosome : chromosomes )
			{
				double currentTime = 0, totalExonReadTime = 0;
				
				System.out.println("Chromosome " + chromosome);
				System.out.println("Calculating expected values from given files....");
				ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, chromosome);
				StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, chromosome, ExpectedValues);
				Exon.sortExons(ExpectedValues);
				Exon.sortExons(StdDeviations);
				
				System.out.println("Reading exon with measurements file....");
				currentTime = System.currentTimeMillis();
				Exons = Exon.readAndStoreExonFile(exonFileName, chromosome);
				totalExonReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);

				// TODO add chromosome checks
				System.out.println("Calculating States....");
				currentTime = System.currentTimeMillis();
				HiddenMarkovModel.getStates(Exons, ExpectedValues, StdDeviations);
				double totalStateCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

				// Print output
				if (outputFileName.equals(""))
					for(int i = 0; i<Exons.size(); i++) //for(Exon e : Exons)
						System.out.println(Exons.get(i) + "\t" + ExpectedValues.get(i).SNPs + "\t" + ExpectedValues.get(i).FPKM);
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
				else if (arg.equals("-b"))
				{
					String[] fields = arguments[index+1].split(",");
					for (String f : fields)
					{
						baselineExonFileNames.add(f);
					}
					
				}
			}

			if (exonFileName.equals("") 
					|| chromosomes.isEmpty()
			) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java CONDR -e <exonFileName> " + 
						"-c <chromosomes (start-end)> " + "[-t] " + "[-o <output file name>]");
				System.exit(0);
			}
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

}
