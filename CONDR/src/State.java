
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
import java.util.ArrayList;
import java.util.HashMap;
import cern.jet.random.engine.RandomEngine;

public class State
{
	int stateName;
	double rpkmRatio;
	double snpRatio;
	int E_LengthOfState;

	static double rateOfOccurenceOfCNV = 0;

	static HashMap<State, HashMap<State, Double>> transitionProbabilities = new HashMap<State, HashMap<State, Double>>();

	// Associating state names with constants
	public static final int NORMAL = 0, HOMOZYGOUS_DELETE = 1, HETEROZYGOUS_DELETE = 2, COPY_NEUTRAL_LOH = 3, 
	INSERTION = 4, CHROMATIN_CHANGE = 5, UNKNOWN_TYPE = 6;
	// Decides whether or not we compute the transition probability at run time or can pre-compute
	public static final int computeAtRunTime = -1;

	/* 
	 * Gets the transition probability from state s1 -> state s2
	 * For CNV state -> normal state, transition probability is dependent on length of intron
	 */
	public static double getTransitionProbability(State s1, State s2, int lengthOfIntron, int lengthOfExon)
	{
		Double transitionProb = 0.0;
		if (s1 == null || s2 == null)
		{
			System.err.println("Invalid states in transition Probabilities");
			System.exit(0);
		}

		if (s1.equals(s2) && s1.stateName!=NORMAL) // self loop in Markov Model
		{
			transitionProb = (Math.max(1 - ( (double)lengthOfIntron / s1.E_LengthOfState ), 0.0));
			// returning a positive number in the event that the length of the intron > expected length of the CNV
		}	
		else if (s1.equals(s2) && s1.stateName == NORMAL) // normal -> normal
		{
			// see below part
			// Poisson process with lambda = rate of occurrence of CNV
			double lambda = rateOfOccurenceOfCNV * lengthOfExon;
			transitionProb = lambda*Math.exp(-lambda);
			transitionProb = 1 - 5*transitionProb;
			if (transitionProb < 0)
				transitionProb = 0.0;
		}
		else if (! s1.equals(s2))
		{
			if (s2.stateName == NORMAL)
			{
				transitionProb = (double)lengthOfIntron / s1.E_LengthOfState;
				// sanity checking in the event that length of the intron > expected length of the CNV
				if (transitionProb < 0)
					transitionProb = 0.0;
				else if (transitionProb > 1.0)
					transitionProb = 1.0;
			}
			else if (s1.stateName == NORMAL)
			{
				// TODO check formula/equation
				double rate = .0001; // TODO: MAKE this into a parameter/user driven input
				// Poisson process with lambda = rate of occurrence of CNV
				double lambda = rate * lengthOfExon;
				transitionProb = lambda*Math.exp(-lambda);

				if (transitionProb < 0)
					transitionProb = 0.0;
				else if (transitionProb > 1.0)
					transitionProb = 1.0;

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
	public static HashMap<String, State> initializeStates()
	{
		//ArrayList<State> states = new ArrayList<State>();
		HashMap<String, State> states = new HashMap<String, State>();

		// Read parameter file and input states
		String parameterFile = "ParameterFile"; // TODO: pass in as command line argument
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
					String fieldName = lineElements[1].split("=")[0].trim();
					String fieldValue = line.split("=")[1].trim();
					// TODO make this neater
					if (fieldName.equals("rpkmRatio"))
						s.rpkmRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("snpRatio"))
						s.snpRatio = Double.parseDouble(fieldValue);
					else if (fieldName.equals("E_LengthOfState"))
						s.E_LengthOfState = Integer.parseInt(fieldValue);
					/*else if (fieldName.equals("E_NumberOfOccurrencesOfState"))
						s.E_NumberOfOccurrencesOfState = Integer.parseInt(fieldValue);
					 */states.put(stateName, s);
				}
			}
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
			case State.UNKNOWN_TYPE: 		return "Unknown Type"; 
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

	public static double getEmissionProbability(double observedFPKM, double observedSNPs, double expectedFPKM,
			double expectedSNPs, double stdDevFPKM, double stdDevSNPs, State s)
	{
		// Emission probability is the probability of observing that value given the distribution parameters
		// the expected values are computed assuming normal. Hence for other states, we multiple the mean
		// by the ratios for that state. We assume the same std dev for all the states distributions


		// Prob(value x | expected, stdDev, state) -> assuming this is a NORMAL DISTRUIBUTION

		cern.jet.random.Normal distFPKM = new cern.jet.random.Normal (s.rpkmRatio*expectedFPKM, stdDevFPKM, RandomEngine.makeDefault());
		double probFPKM = distFPKM.pdf(observedFPKM);
		cern.jet.random.Normal distSNPs = new cern.jet.random.Normal (s.rpkmRatio*expectedSNPs, stdDevSNPs, RandomEngine.makeDefault());
		double probSNPs = distSNPs.pdf(observedSNPs);

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

}
