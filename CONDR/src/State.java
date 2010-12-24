
/*
 * This class defines the states included in the Markov Model
 * They are based on underlying biology and as such we define six states
 * 
 * NORMAL
 * HOMOZYGOUS_DELETE
 * HETEROZYGOUS_DELETE
 * COPY_NEUTRAL_LOH
 * INSERTION
 * CHROMATIN_CHANGE
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;

public class State
{
	int stateName;
	double rpkmRatio;
	double snpRatio;
	int E_LengthOfState;
	//int E_exons;

	static double rateOfOccurenceOfCNV = 0;

	static HashMap<State, HashMap<State, Double>> transitionProbabilities = new HashMap<State, HashMap<State, Double>>();

	// Associating state names with constants
	public static final int NORMAL = 0, HOMOZYGOUS_DELETE = 1, HETEROZYGOUS_DELETE = 2, COPY_NEUTRAL_LOH = 3, 
	INSERTION = 4, CHROMATIN_CHANGE = 5, UNKNOWN_TYPE = 6;
	// Decides whether or not we compute the transition probability at run time or can pre-compute
	public static final int computeAtRunTime = -1;
	public static HashMap<Integer, Double> LOG_FACTORIAL = new HashMap<Integer, Double>();

	/* 
	 * Gets the transition probability from state s1 -> state s2
	 * For CNV state -> normal state, transition probability is dependent on length of intron
	 */
	public static double getTransitionProbability(State s1, State s2, int lengthOfIntron, int lengthOfExon)
	{
		// TODO see if this an appropriate fix
		double intronLength = 0;
		if ( lengthOfIntron < 0 )
			intronLength = 0.0000001;
		else
			intronLength = (double)lengthOfIntron;
		//System.out.println("Rate of occurence: " + rateOfOccurenceOfCNV );
		Double transitionProb = 0.0;
		if (s1 == null || s2 == null)
		{
			System.err.println("Invalid states in transition Probabilities");
			System.exit(0);
		}

		if (s1.equals(s2) && s1.stateName!=NORMAL) // self loop in Markov Model
		{
			//TODO: sigh yet another hack
			transitionProb = (Math.max(1 - ( (double)intronLength / s1.E_LengthOfState ), 0));
			// returning a positive number in the event that the length of the intron > expected length of the CNV
		}	
		else if (s1.equals(s2) && s1.stateName == NORMAL) // normal -> normal
		{
			// see below part
			// Poisson process with lambda = rate of occurrence of CNV
			double lambda = rateOfOccurenceOfCNV * lengthOfExon;
			transitionProb = lambda*Math.exp(-lambda);
			transitionProb = 1 - 5*transitionProb;
			if (transitionProb < 0) // TODO is this correct?
				transitionProb = 0.0;
		}
		else if (! s1.equals(s2))
		{
			if (s2.stateName == NORMAL)
			{
				transitionProb = (double)intronLength / s1.E_LengthOfState;
				// sanity checking in the event that length of the intron > expected length of the CNV
				if (transitionProb < 0)
					transitionProb = 0.0;
				else if (transitionProb > 1.0)
					transitionProb = 1.0;
			}
			else if (s1.stateName == NORMAL)
			{
				// TODO check formula/equation
				// Poisson process with lambda = rate of occurrence of CNV
				double lambda = rateOfOccurenceOfCNV * lengthOfExon;
				transitionProb = lambda*Math.exp(-lambda);

				/*
				if (transitionProb < 0)
					transitionProb = 0.0;
				else if (transitionProb > 1.0)
					transitionProb = 1.0;
				 */
			}
			else transitionProb = 0.0;
		}
		else
			transitionProb = 0.0;

		return transitionProb;
	}

	/*
	 * Initializes all the state information: expected RPKM, SNP values and the allowed deviations from expectation
	 * This information is based on prior observations
	 * Contains Hard coded assumptions about the data
	 */
	public static HashMap<String, State> initializeStates(String parameterFile)
	{
		HashMap<String, State> states = new HashMap<String, State>();

		// Read parameter file and input states
		String line = "";
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(parameterFile));
			while( (line = br.readLine()) != null)
			{
				if (line.startsWith("//") || line.equals("")) // ignore comments
					continue;
				else if (line.startsWith("rate"))
				{
					rateOfOccurenceOfCNV = Double.parseDouble(line.split("=")[1]);
				}
				else
				{
					String[] lineElements = line.split("\\."); 
					String stateName = lineElements[0];
					State s = null;
					if (states.containsKey(stateName))
						s = states.get(stateName);
					else
					{
						s = new State();
						if (stateName.equals("NORMAL"))
							s.stateName = State.NORMAL;
						else if (stateName.equals("HOMOZYGOUS_DELETE"))
							s.stateName = State.HOMOZYGOUS_DELETE;
						else if (stateName.equals("HETEROZYGOUS_DELETE"))
							s.stateName = State.HETEROZYGOUS_DELETE;
						else if (stateName.equals("COPY_NEUTRAL_LOH"))
							s.stateName = State.COPY_NEUTRAL_LOH;
						else if (stateName.equals("INSERTION"))
							s.stateName = State.INSERTION;
						else if (stateName.equals("CHROMATIN_CHANGE"))
							s.stateName = State.CHROMATIN_CHANGE;
					}
					//System.out.println(stateName);
					String fieldName = lineElements[1].split("=")[0].trim();
					String fieldValue = line.split("=")[1].trim();
					// TODO make this neater
					if (fieldName.equals("rpkmRatio"))
						s.rpkmRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("snpRatio"))
						s.snpRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("E_LengthOfState"))
						s.E_LengthOfState = Integer.parseInt(fieldValue);
					states.put(stateName, s);
				}			
			}
			// TODO: maybe not correct? Change so it different for each state
			rateOfOccurenceOfCNV = (double)states.get("HOMOZYGOUS_DELETE").E_LengthOfState/3000000000.0	;//states.get("NORMAL").E_LengthOfState;
			//System.out.println("Rate of occurrence: " + State.rateOfOccurenceOfCNV + "\t" + 
			//		states.get("HOMOZYGOUS_DELETE").E_LengthOfState + "\t" + states.get("NORMAL").E_LengthOfState);
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process parameter file");
			e.printStackTrace();
			System.exit(0);
		}

		return states;
	}

	// For pretty-printing and for serving coffee with the code
	public String getStateName()
	{
		switch(this.stateName)
		{
			case State.NORMAL: 				return "Normal"; 
			case State.HOMOZYGOUS_DELETE: 	return "Homozygous Deletion"; 
			case State.HETEROZYGOUS_DELETE: return "Heterozygous Deletion"; 
			case State.COPY_NEUTRAL_LOH: 	return "Copy Neutral Loss of Heterozygosity"; 
			case State.INSERTION: 			return "Insertion"; 
			case State.CHROMATIN_CHANGE: 	return "Chromatin Change"; 
			//case State.UNKNOWN_TYPE: 		return "Unknown Type"; 
		}
		return "";
	}

	// Printing functions
	public String toString()
	{
		String stateName = this.getStateName();
		return stateName + "\t" + this.rpkmRatio + "\t" + this.snpRatio + "\t";
	}

	// Defining our equality
	@Override public boolean equals(Object aThat) {
		if ( this == aThat ) return true;

		if ( !(aThat instanceof State) ) return false;
		State that = (State)aThat;

		return	(this.stateName == that.stateName);
	}

	// TODO: remove.. most likely useless function
	public static double getEmissionProbability(double observedFPKM, double observedSNPs, double expectedFPKM,
			double expectedSNPs, double stdDevFPKM, double stdDevSNPs, double normalizationFPKM, double normalizationSNPs, State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming normal. Hence for other states, we multiple the mean
		// by the ratios for that state. We assume the same std dev for all the states distributions

		// TODO: ChECK THE NROMALIZATION!!!!!

		// Prob(value x | expected, stdDev, state) -> assuming this is a NORMAL DISTRUIBUTION

		/*
		//cern.jet.random.Normal distFPKM = new cern.jet.random.Normal (s.rpkmRatio*expectedFPKM, stdDevFPKM, RandomEngine.makeDefault());
		Poisson distFPKM = new Poisson ((int)(s.rpkmRatio*expectedFPKM), RandomEngine.makeDefault());
		//System.out.println(observedFPKM/normalizationFPKM);
		double probFPKM = distFPKM.pdf((int)observedFPKM);
		//cern.jet.random.Normal distSNPs = new cern.jet.random.Normal (s.rpkmRatio*expectedSNPs, stdDevSNPs, RandomEngine.makeDefault());
		Poisson distSNPs = new Poisson ((int)(s.rpkmRatio*expectedSNPs), RandomEngine.makeDefault());
		//System.out.println(observedSNPs/normalizationSNPs + "\t" + observedSNPs + "\t" + normalizationSNPs + "\t" + (int)(observedSNPs/normalizationSNPs) + "\t" + (int)(s.rpkmRatio*expectedSNPs/normalizationSNPs));
		//System.out.println(observedSNPs + "\t" + normalizationSNPs + "\t" + (int)(s.rpkmRatio*expectedSNPs));
		double probSNPs = distSNPs.pdf((int)observedSNPs);
		 */

		observedFPKM = observedFPKM/normalizationFPKM;
		expectedFPKM = expectedFPKM/normalizationFPKM;
		observedSNPs = observedSNPs/normalizationSNPs;
		expectedSNPs = expectedSNPs/normalizationSNPs;

		//System.out.println(observedFPKM + "\t" + expectedFPKM + "\t" + observedSNPs + "\t" + expectedSNPs);
		// STILL NEED TO FIX NORMALIZATION
		// e^(-lambda) * lambda^k / k!
		//System.out.println(observedSNPs + 		"\t" + normalizationSNPs + "\t" + (int)(s.rpkmRatio*expectedSNPs));
		int lambda = (int)(s.rpkmRatio*expectedSNPs);
		double probSNPs = Math.exp((double)lambda) * Math.pow(lambda, (int)observedSNPs) / cern.jet.math.Arithmetic.factorial((int)observedSNPs);
		lambda = (int)(s.rpkmRatio*expectedFPKM);
		double probFPKM = Math.exp((double)lambda) * Math.pow(lambda, (int)observedFPKM) / cern.jet.math.Arithmetic.factorial((int)observedFPKM);

		//System.out.println(probFPKM + "\t" + observedFPKM + "\t" + probSNPs + "\t" + observedSNPs);

		// handling case where all the baselines are equal
		if (stdDevFPKM == 0.0)
			if (observedFPKM == s.rpkmRatio*expectedFPKM)
				probFPKM = 1;
			else 
				probFPKM = 0;
		if (stdDevSNPs == 0.0)
			if (observedSNPs == s.snpRatio*expectedSNPs)
				probSNPs = 1;
			else
				probSNPs = 0;

		// TODO: how to combine the two?
		return (probFPKM * probSNPs);
	}

	public static Double getEmissionProbability(Exon exon, Exon expected, Exon stdDev, State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming poisson. Hence for other states, we multiple the mean
		// by the ratios for that state. 

		double observedSNPs = exon.SNPs;
		double expectedSNPs = expected.SNPs; 
		double observedFPKM = exon.FPKM;
		double expectedFPKM = expected.FPKM;

		// TODO: Are std dev's necessary?
		double stdDevFPKM = stdDev.FPKM;
		double stdDevSNPs = stdDev.SNPs;

		Poisson distFPKM = new Poisson (s.rpkmRatio*expectedFPKM, RandomEngine.makeDefault());
		double probFPKM = distFPKM.pdf((int)observedFPKM);
		Poisson distSNPs = new Poisson (s.rpkmRatio*expectedSNPs, RandomEngine.makeDefault());
		double probSNPs = distSNPs.pdf((int)observedSNPs);

		// Since the package handles 0! = 0 instead of 0! = 1, we have to do that manually		
		if (observedSNPs == 0.0)
		{
			double lambda = s.rpkmRatio*expectedSNPs;
			probSNPs = Math.exp(-lambda) * Math.pow(lambda, (int)observedSNPs) / cern.jet.math.Arithmetic.factorial((int)observedSNPs);
		}
		if (observedFPKM == 0.0)
		{
			double lambda = s.rpkmRatio*expectedFPKM;
			probFPKM = Math.exp(-lambda) * Math.pow(lambda, (int)observedFPKM) / cern.jet.math.Arithmetic.factorial((int)observedFPKM);
		}	

		/*
		double lambdaSNPs = s.rpkmRatio*expectedSNPs;
		double probSNPs = -lambdaSNPs * ((int)observedSNPs * Math.log(lambdaSNPs)) - logFactorial((int)observedSNPs); //cern.jet.math.Arithmetic.factorial((int)observedSNPs);
		double lambdaFPKM = s.rpkmRatio*expectedFPKM;
		double probFPKM = -lambdaFPKM * ((int)observedFPKM * Math.log(lambdaFPKM)) - logFactorial((int)observedFPKM); //cern.jet.math.Arithmetic.factorial((int)observedFPKM);
		 */

		// handling case where all the baselines are equal
		if (s.rpkmRatio*expectedFPKM == 0.0)
			if (observedFPKM == s.rpkmRatio*expectedFPKM)
				probFPKM = 1;
			else 
				probFPKM = 0;
		if (s.snpRatio*expectedSNPs == 0.0)
			if (observedSNPs == s.snpRatio*expectedSNPs)
				probSNPs = 1;
			else
				probSNPs = 0;


		if (Double.isInfinite(probSNPs) || Double.isNaN(probSNPs)  || probSNPs == 0.0)
		{
			//System.out.println("!!!");
			probSNPs = 10e-30;
		}
		if (Double.isInfinite(probFPKM) || Double.isNaN(probFPKM) || probFPKM == 0.0)
		{
			//System.out.println("!!!");
			probFPKM = 10e-30;
		}


		// TODO: how to combine the two?
		return (probFPKM * probSNPs);
	}

	private static double logFactorial(int number)
	{
		// storing some numbers for quicker computation
		Double value = LOG_FACTORIAL.get(number);
		if (value == null)
		{
			value = 0.0;
			for (int i=1; i<=number; i++)
			{
				value += Math.log(i);
			}
			LOG_FACTORIAL.put(number, value);
		}
		return value;
	}

	public static double getLogEmissionProbability(Exon exon, Exon expected, Exon stdDev,
			State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming poisson. Hence for other states, we multiple the mean
		// by the ratios for that state. 

		double observedSNPs = exon.SNPs;
		double expectedSNPs = expected.SNPs; 
		double observedFPKM = exon.FPKM;
		double expectedFPKM = expected.FPKM;

		double lambdaSNPs = s.snpRatio*expectedSNPs;
		double logProbSNPs = -lambdaSNPs + ((int)observedSNPs * Math.log(lambdaSNPs)) - logFactorial((int)observedSNPs);
		double lambdaFPKM = s.rpkmRatio*expectedFPKM;
		double logProbFPKM = -lambdaFPKM + ((int)observedFPKM * Math.log(lambdaFPKM)) - logFactorial((int)observedFPKM);


		// handling case where all the baselines are equal, since log(0) = Inf 
		if (s.rpkmRatio*expectedFPKM == 0.0)
			if (observedFPKM == s.rpkmRatio*expectedFPKM)
				logProbFPKM = 0; // log 1 = 0; expected = observed = 0 so it fits with the expectation
			else
				logProbFPKM = Double.NEGATIVE_INFINITY;
		//else 
		//logProbFPKM = -4000000; // log 0 = ? TODO
		if (s.snpRatio*expectedSNPs == 0.0)
			if (observedSNPs == s.snpRatio*expectedSNPs)
				logProbSNPs = 0;
			else
				logProbSNPs = Double.NEGATIVE_INFINITY;
		//else
		//logProbSNPs = -4000000; // log 0 = ?

		if (exon.posLeft == 47950286)
			System.out.println("$$\t" + s.stateName + "\t" + (logProbFPKM + logProbSNPs));
		/*
		if (Double.isInfinite(probSNPs) || Double.isNaN(probSNPs) /* || probSNPs == 0.0)
		{
			System.out.println("!!!");
			probSNPs = 10e-30;
		}
		if (Double.isInfinite(probFPKM) || Double.isNaN(probFPKM) /* || probFPKM == 0.0*)
		{
			System.out.println("!!!");
			probFPKM = 10e-30;
		}
		 */
		// TODO: how to combine the two?
		return (logProbFPKM + logProbSNPs);
	}

}
