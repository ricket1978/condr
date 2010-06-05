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

import java.util.ArrayList;
import java.util.HashMap;

public class State
{
	int stateName;
	double rpkmRatio;
	double snpRatio;
	int E_LengthOfState;
	int E_NumberOfOccurrencesOfState;

	static double deltaRPKM;
	static double deltaSNPs;
	static double stdDevRPKM;
	static double E_RPKM;
	static double E_SNPs;

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
	public static double getTransitionProbability(State s1, State s2, int lengthOfIntron)
	{
		Double transitionProb = State.transitionProbabilities.get(s1).get(s2);

		// Check for invalid state transition. Prob(invalid state transition) = 0
		if (transitionProb == null)
			transitionProb = 0.0;
		
		if ( transitionProb == State.computeAtRunTime )
		{
			if (s1.equals(s2)) // self loop in Markov Model
			{
				transitionProb = (Math.max(1 - ( (double)lengthOfIntron / s1.E_LengthOfState ), 0.0));
				// returning a positive number in the event that the length of the intron > expected length of the CNV
			}	
			else
			{
				transitionProb = (double)lengthOfIntron / s1.E_LengthOfState;
				// sanity checking in the event that length of the intron > expected length of the CNV
				if (transitionProb < 0)
					transitionProb = 0.0;
				else if (transitionProb > 1.0)
					transitionProb = 1.0;
			}
		}
		return transitionProb;
	}

	/* 
	 * Gets the probability of emitting the observations for RPKM and SNP from state s
	 */
	public static double getEmissionProbability(double observedRPKM, double observedSNPs, State s)
	{
		double observedRPKMRatio = observedRPKM/State.E_RPKM;
		double observedSNPRatio = observedSNPs/State.E_SNPs;

		// Bin the observed values between (Expectation - Delta) and (Expectation + Delta)
		// We have allowed a Delta deviation from the expectation to allow for 
		// less stringent cut-off values
		// Expectation & Delta values are computed from prior observations such that the 
		// probability of being within that range is 95%
		// Thus, if the observation is within Delta of the Expectation, there is a 95% 
		// probability of it having been from a distribution with the state's parameters
		// If not, it is 5%
		if (observedRPKMRatio >= s.rpkmRatio - State.deltaRPKM && 
				observedRPKMRatio <= s.rpkmRatio + State.deltaRPKM &&
				observedSNPRatio >= s.snpRatio - State.deltaSNPs &&
				observedSNPRatio <= s.snpRatio + State.deltaSNPs )
			return 0.95;
		else
			return 0.05;
	}

	/*
	 * Initializes all the state information: expected RPKM, SNP values and the allowed deviations from expectation
	 * This information is based on prior observations
	 * Contains Hard coded assumptions about the data
	 */
	public static ArrayList<State> initializeStates(double expected_RPKM, double expected_SNPs)
	{
		ArrayList<State> states = new ArrayList<State>();
		State.E_RPKM = expected_RPKM;
		State.E_SNPs = expected_SNPs;
		
		State.deltaRPKM = 0.25;
		State.deltaSNPs = 0.05;

		State s = new State();
		s.stateName = State.NORMAL;
		s.rpkmRatio = 1;
		s.snpRatio = 1;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 1000;
		states.add(s);

		s = new State();
		s.stateName = State.HOMOZYGOUS_DELETE;
		s.rpkmRatio = 0;
		s.snpRatio = 0;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 200;
		states.add(s);

		s = new State();
		s.stateName = State.HETEROZYGOUS_DELETE;
		s.rpkmRatio = .5;
		s.snpRatio = 0;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 200;
		states.add(s);

		s = new State();		
		s.stateName = State.COPY_NEUTRAL_LOH;
		s.rpkmRatio = 1;
		s.snpRatio = 0;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 200;
		states.add(s);

		s = new State();		
		s.stateName = State.INSERTION;
		s.rpkmRatio = 1.5;
		s.snpRatio = 1.5;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 200;
		states.add(s);

		s = new State();		
		s.stateName = State.CHROMATIN_CHANGE;
		s.rpkmRatio = 1.5;
		s.snpRatio = 1;
		s.E_LengthOfState = 2000;
		s.E_NumberOfOccurrencesOfState = 200;
		states.add(s);

		intializeTransitionProb(states);
		return states;
	}

	/* 
	 * Computes the Transition Probability matrix
	 * 
	 * If we can compute ahead of time, we put into the matrix.
	 * Otherwise, Prob = -1 and we compute when we have the necessary information
	 */
	private static void intializeTransitionProb(ArrayList<State> states)
	{
		int totalNumberOfOccurrencesOfStates = 0;
		for(int i=0; i<states.size(); i++)
			totalNumberOfOccurrencesOfStates += states.get(i).E_NumberOfOccurrencesOfState;

		HashMap<State, Double> transProbNormalToCNV = new HashMap<State, Double>();
		for (State s : states)
		{
			double transProb = 0;
			transProb = (double)s.E_NumberOfOccurrencesOfState / totalNumberOfOccurrencesOfStates;
			transProbNormalToCNV.put(s, transProb);
		}

		HashMap<State, Double> transProbCNVToNormal = new HashMap<State, Double>();
		for (State s : states)
		{
			// Normal -> Normal taken care of in previous transition probability list
			if (s.stateName == NORMAL)
				continue;
			double transProb = State.computeAtRunTime;
			transProbCNVToNormal.put(s, transProb);
		}

		// insert into transition probability matrix
		State normalState = states.get(NORMAL);
		for(State s1 : states)
		{
			if (s1.stateName == NORMAL)
				State.transitionProbabilities.put(s1, transProbNormalToCNV);
			else
			{
				HashMap<State, Double> transProbStateToOthers = new HashMap<State, Double>();			
				transProbStateToOthers.put(normalState, transProbCNVToNormal.get(s1));
				transProbStateToOthers.put(s1, (double)State.computeAtRunTime);			
				State.transitionProbabilities.put(s1, transProbStateToOthers);
			}
		}
	}

	// For pretty-printing and for serving coffee with the code
	public String getStateName()
	{
		switch(this.stateName)
		{
			case State.NORMAL: return "Normal"; 
			case State.HOMOZYGOUS_DELETE: return "Homozygous Deletion"; 
			case State.HETEROZYGOUS_DELETE: return "Heterozygous Deletion"; 
			case State.COPY_NEUTRAL_LOH: return "Copy Neutral Loss of Heterozygosity"; 
			case State.INSERTION: return "Insertion"; 
			case State.CHROMATIN_CHANGE: return "Chromatin Change"; 
			case State.UNKNOWN_TYPE: return "Unknown Type"; 
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
}
