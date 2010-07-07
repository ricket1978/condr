import java.util.ArrayList;

import cern.jet.random.engine.RandomEngine;


public class SimulateExonMeasurements
{
	static ArrayList<String> baselineExonFileNames = new ArrayList<String>();
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();

	public static void main(String args[])
	{
		ArrayList<Exon> exons = new ArrayList<Exon>();
		HiddenMarkovModel.initialize();

		parseArguments( args );

		for(int i=0; i<60; i++)
		{
			Exon e = new Exon();

			if (i<10)
			{
				e.state = HiddenMarkovModel.States.get(0);
			}
			else if (i<20)
			{
				e.state = HiddenMarkovModel.States.get(1);
			}
			else if (i<30)
			{
				e.state = HiddenMarkovModel.States.get(2);
			}
			else if (i<40)
			{
				e.state = HiddenMarkovModel.States.get(3);
			}
			else if (i<50)
			{
				e.state = HiddenMarkovModel.States.get(4);
			}
			else if (i<60)
			{
				e.state = HiddenMarkovModel.States.get(5);
			}			
			exons.add(e);
		}

		ArrayList<Exon> ExpectedValues = new ArrayList<Exon>();
		ArrayList<Exon> StdDeviations = new ArrayList<Exon>(); 

		ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, 1);
		StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, 1, ExpectedValues);

		for(int i=0; i<exons.size(); i++)
		{
			Exon e = exons.get(i);
			e.SNPs = getLikelyValue(ExpectedValues.get(i).SNPs*e.state.snpRatio, StdDeviations.get(i).SNPs);
			e.FPKM = getLikelyValue(ExpectedValues.get(i).FPKM*e.state.rpkmRatio, StdDeviations.get(i).FPKM);
			//			System.out.println(e + "\t" + ExpectedValues.get(i).SNPs + "\t" + ExpectedValues.get(i).FPKM + "\t" + ExpectedValues.get(i).FPKM*e.state.rpkmRatio);
			System.out.println(e);
		}

	}

	private static double getLikelyValue(double mean, double stddev)
	{
		// TODO: check calculation of lambda parameter
		// pick a value from poisson distribution with those parameters
		cern.jet.random.Poisson dist = new cern.jet.random.Poisson(mean, RandomEngine.makeDefault());
		return dist.nextDouble();
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
				if (arg.equals("-c"))
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

			if (chromosomes.isEmpty()) 
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
