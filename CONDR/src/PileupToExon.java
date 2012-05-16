import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;


public class PileupToExon
{
	// TODO: Convert to script maybe? Or leave as Java?
	/*
	 * File Format of output (tab delimited)
	 * chromosome exonStart exonEnd geneName FPKM/coverage SNP levels/BAF
	 * 
	 * .pileup -> .exon
	 * 
	 */
	static ArrayList<Exon> Exons = new ArrayList<Exon>();

	static String exonCaptureArrayFileName = "";
	static String pileupFileName = "";
	static String outputFileName = "";
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;

	/*
	 * arguments:
	 * -e "./exonList/exonParsedUniqChr16" -pileup "WZ1034.region16.out" -c "16-16" -o "WZ1034.region16.exon" -t 
	 */

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();

		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal/genomic order
		 */
		try
		{
                    	double currentTime = 0, totalExonReadTime = 0, totalReadsReadTime = 0;
			BufferedReader br = null;
			
			if (pileupFileName == "")
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(pileupFileName));

			
			for (int chromosome : chromosomes )
			{
				System.out.println("Chromosome " + chromosome);
				System.out.println("Reading exon capture file....");
				currentTime = System.currentTimeMillis();
				String exonFileNameByChr = "";
				if (chromosomes.size() == 1)
					exonFileNameByChr = exonCaptureArrayFileName;
				else
					/*
					 * need to be able to handle exons by separate files
					 * need to be able to handle exons and one full file
					 * any specification of names should take place outside of this program
					 */
					// TODO: hardcoding in here
					exonFileNameByChr = exonCaptureArrayFileName + ".chr" + chromosome;
                                System.out.println("Exon file name: "+ exonFileNameByChr);
                                Exons = new ArrayList<Exon>();
				Exons = Exon.readExon(exonFileNameByChr, chromosome);
				totalExonReadTime += (System.currentTimeMillis() - currentTime)/1000F;
				//Exon.sortExons(Exons);

				// TODO add chromosome checks
				if (usingPileup)
				{
					System.out.println("Reading pileup file, calculating coverage and SNPs....");
					currentTime = System.currentTimeMillis();
					Pileup.readData(Exons, br);
					totalReadsReadTime += (System.currentTimeMillis() - currentTime)/1000F;
				}

				// Print output
				Writer output = new BufferedWriter(new FileWriter(outputFileName+".chr"+chromosome));
				System.out.println("Output file name:" + outputFileName+".chr"+chromosome);
				for(Exon e : Exons)
				{
					//System.out.println(e);
					output.write(e + "\n");
				}
				output.close();
			}
			// prints the timing metrics to std out
			if (printTimingMetrics)
			{
				double endTime = System.currentTimeMillis();
				totalTime = (endTime - startTime)/1000F;
				System.out.println("Total Time: " + totalTime);
				System.out.println("Time for reading exons file       : " + totalExonReadTime + ", " + Exons.size());
				System.out.println("Time for reading mapped reads file: " + totalReadsReadTime);
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
		 * -pileup <pileup file>
		 * -c <chromosomes (start-end)>
		 * -t (prints timing metrics)
		 * -o <output file name>
		 */
		try {
			for(int index = 0; index < arguments.length; index ++)
			{
				String arg = arguments[index];
				if (arg.equals("-e"))
					exonCaptureArrayFileName = arguments[index + 1];
				else if (arg.equals("-pileup"))
				{
					usingPileup = true;
					pileupFileName = arguments[index + 1];
				}
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
						//chromosomes.add(fields[0].charAt(0));
					}
				}
			}

			if (exonCaptureArrayFileName.equals("") 
					|| chromosomes.isEmpty()
					|| outputFileName.equals("")) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java ConvertToExonFormat -e <exonFileName> " +  
						"-pileup <pileup fileName> " + 
						"-c <chromosomes (start-end)> " + "-o <output file name>");
				System.out.println(exonCaptureArrayFileName);
				System.out.println(chromosomes);
				System.out.println(outputFileName);
				System.out.println(pileupFileName);
				System.exit(0);
				//throw(new Exception("Improper Usage"));
			}
			// default is the entire genome
			//if (chromosomes.isEmpty())
				//for(int c = 1; c<= 22; c++)
					//chromosomes.add(c);
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

}
