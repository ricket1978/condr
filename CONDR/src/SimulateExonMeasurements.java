import java.util.ArrayList;

import cern.jet.random.engine.RandomEngine;
import cern.jet.random.*;

public class SimulateExonMeasurements
{
	// W3 created with variable lengths with means based on the parameter file
	static ArrayList<String> baselineExonFileNames = new ArrayList<String>();
	static ArrayList<Integer> chromosomes = new ArrayList<Integer>();
	static ArrayList<Exon> exons = new ArrayList<Exon>();

	public static void main(String args[])
	{
		HiddenMarkovModel.initialize( "SimulationParameterFile" );

		parseArguments( args );
		ArrayList<Exon> ExpectedValues = new ArrayList<Exon>();
		ArrayList<Exon> StdDeviations = new ArrayList<Exon>(); 

		ExpectedValues = Exon.calculateExpectedValues(baselineExonFileNames, 1);
		StdDeviations = Exon.calculateStdDevValues(baselineExonFileNames, 1, ExpectedValues);

		// assigning state values
		// TODO: switch the ordering and such while testing since some states are smaller than others as a result
		int index = 0;
		while( true)
		{
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "HOMOZYGOUS_DELETE", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "HETEROZYGOUS_DELETE", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "COPY_NEUTRAL_LOH", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "INSERTION", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "NORMAL", index);
			if (index > ExpectedValues.size() - 10) break;
			index = SimulateRegion(ExpectedValues, StdDeviations, "CHROMATIN_CHANGE", index);
			if (index > ExpectedValues.size() - 10) break;
		}
		for(int i=0; i<exons.size(); i++)
		{
			Exon e = exons.get(i);
			e.SNPs = getLikelyValue(ExpectedValues.get(i).SNPs*e.state.snpRatio, StdDeviations.get(i).SNPs);
			e.FPKM = getLikelyValue(ExpectedValues.get(i).FPKM*e.state.rpkmRatio, StdDeviations.get(i).FPKM);
			//System.out.println(e + "\t" + ExpectedValues.get(i).SNPs + "\t" + ExpectedValues.get(i).FPKM*e.state.rpkmRatio
			//	+ "\t" + StdDeviations.get(i).SNPs + "\t" + StdDeviations.get(i).FPKM);
			System.out.println(e);
		}

	}

	private static int SimulateRegion(ArrayList<Exon> ExpectedValues,
			ArrayList<Exon> StdDeviations, String stateName, int index)
	{
		// TODO: check about Poisson distribution for this
		int regionStartPosition = ExpectedValues.get(index).posLeft;
		Poisson dist = new Poisson(HiddenMarkovModel.States.get(stateName).E_LengthOfState, RandomEngine.makeDefault());
		int regionLength = dist.nextInt();
		Exon e = ExpectedValues.get(index);
		int genomicLength = 0;
		for(genomicLength = 0; genomicLength <= regionLength; genomicLength = (e.posRight - regionStartPosition))
		{
			e.state = HiddenMarkovModel.States.get(stateName);
			exons.add(e);
			index ++;
			e = new Exon();
			e = ExpectedValues.get(index);
		}
		return index;
	}

	private static double getLikelyValue(double mean, double stddev)
	{
		// TODO: check calculation of lambda parameter
		// pick a value from poisson distribution with those parameters
		//cern.jet.random.Poisson dist = new cern.jet.random.Poisson(mean, RandomEngine.makeDefault());
		Normal dist = new Normal(mean, stddev, RandomEngine.makeDefault());
		double value = Double.NEGATIVE_INFINITY;
		while (value < 0)
			value = dist.nextDouble();
		return value;
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
