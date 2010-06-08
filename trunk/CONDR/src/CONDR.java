import java.io.*;
import java.util.*;


public class CONDR
{
	static ArrayList<Exon> Exons = new ArrayList<Exon>();

	static String exonFileName = "";
	static String expressionFileName = "";
	static String mappedReadsFileName = "";
	static String referenceDirectory = "";
	static String outputFileName = "";
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static boolean printTimingMetrics = false;
	static boolean usingPileup = false;

	public static void main(String args[])
	{
		double totalTime;
		double startTime = System.currentTimeMillis();


		parseArguments( args );

		// preprocess input files to find the chromosomal boundaries in terms of line number
		/*
		 * All files should be sorted in chromosomal order
		 */
		try
		{
			for (int chromosome : chromosomes )
			{
				double currentTime = 0, totalExonReadTime = 0, totalExprReadTime = 0, totalFPKMCalcTime = 0, totalReadsReadTime = 0, totalSNPsCalcTime = 0, totalRefCalcTime = 0;
				int numberOfExpr = 0, numberOfReads = 0;

				System.out.println("Chromosome " + chromosome);
				System.out.println("Reading exons file....");
				currentTime = System.currentTimeMillis();
				Exons = Exon.readExon(exonFileName, chromosome);
				totalExonReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				Exon.sortExons(Exons);

				// TODO add chromosome checks
				if (usingPileup)
				{
					System.out.println("Reading pileup file, calculating coverage and SNPs....");
					currentTime = System.currentTimeMillis();
					Pileup.readData(Exons, mappedReadsFileName);
					totalReadsReadTime = (System.currentTimeMillis() - currentTime)/1000F;
				}
				else
				{
					System.out.println("Reading expression file....");
					ArrayList<Expression> Expressions = new ArrayList<Expression>();
					currentTime = System.currentTimeMillis();
					Expressions = Expression.readExon(expressionFileName, chromosome);
					totalExprReadTime = (System.currentTimeMillis() - currentTime)/1000F;
					numberOfExpr = Expressions.size();

					System.out.println("Calculating FPKMs....");
					currentTime = System.currentTimeMillis();
					Exon.getFPKM(Expressions, Exons);
					totalFPKMCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
					Expressions.removeAll(Expressions); // explicitly deleting to free up memory

					System.out.println("Reading mapped reads SAM file....");
					ArrayList<MappedReads> mappedReads = new ArrayList<MappedReads>();
					currentTime = System.currentTimeMillis();
					mappedReads = MappedReads.readMappedReads(mappedReadsFileName, chromosome);
					totalReadsReadTime = (System.currentTimeMillis() - currentTime)/1000F;
					MappedReads.sort(mappedReads);
					numberOfReads = mappedReads.size();

					System.out.println("Reading reference genome file....");
					String referenceFileName = referenceDirectory + "/chr" + chromosome + ".fa";
					currentTime = System.currentTimeMillis();
					RandomAccessFile inputReference = new RandomAccessFile(referenceFileName, "r");
					for(Exon e : Exons)
					{
						e.getReferenceSequence(inputReference);
					}
					totalRefCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

					System.out.println("Calculating SNPs....");
					currentTime = System.currentTimeMillis();
					Exon.getSNPs(Exons, mappedReads); 
					totalSNPsCalcTime = (System.currentTimeMillis() - currentTime)/1000F;
					mappedReads.removeAll(mappedReads);
				}

				System.out.println("Calculating States....");
				currentTime = System.currentTimeMillis();
				HiddenMarkovModel.getStates(Exons);
				double totalStateCalcTime = (System.currentTimeMillis() - currentTime)/1000F;

				// Print output
				if (outputFileName.equals(""))
				{
					// print to stdout
					for(Exon e : Exons)
						System.out.println(e);
				}
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
					System.out.println("Time for reading expression file  : " + totalExprReadTime + ", " + numberOfExpr);
					System.out.println("Time for reading mapped reads file: " + totalReadsReadTime + ", " + numberOfReads);
					System.out.println("Time for getting reference seq    : " + totalRefCalcTime);
					System.out.println("Time for calculating FPKM         : " + totalFPKMCalcTime);
					System.out.println("Time for calculating Num of SNPs  : " + totalSNPsCalcTime);
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
		 * -c <chromosomes (start-end)>
		 * -t (prints timing metrics)
		 * -o <output file name>
		 */
		try {
			for(int index = 0; index < arguments.length; index ++)
			{
				String arg = arguments[index];
				if (arg.equals("-e"))
					exonFileName = arguments[index + 1];
				else if (arg.equals("-ex"))
					expressionFileName = arguments[index + 1];
				else if (arg.equals("-sam"))
					mappedReadsFileName = arguments[index + 1];
				else if (arg.equals("-pileup"))
				{
					usingPileup = true;
					mappedReadsFileName = arguments[index + 1];
				}
				else if (arg.equals("-ref"))
					referenceDirectory = arguments[index + 1];
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
			}

			if (exonFileName.equals("") 
					|| referenceDirectory.equals("") || chromosomes.isEmpty()
					|| ((!usingPileup && (expressionFileName.equals("") || mappedReadsFileName.equals(""))) 
							&& (usingPileup && mappedReadsFileName.equals("")))) 
			{
				System.err.println("Improper Usage:");
				System.err.println("java exonSeg -e <exonFileName> " + "-ex <expressions file> " + 
						"-sam <mapped reads, sam file> " + "-ref <reference directory> " + 
						"-c <chromosomes (start-end)> " + "[-t] " + "[-o <output file name>]");
				System.exit(0);
				//throw(new Exception("Improper Usage"));
			}
		} catch(Exception e)
		{
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}

	}

}
