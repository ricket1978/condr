/*
 * Calculates and stores states according to the exon data
 */

import java.util.ArrayList;

public class HiddenMarkovModel
{
	static ArrayList<State> States = new ArrayList<State>();

	// constant: helps maintain minimum value for probability that we want to store
	public static final double UNDERFLOW = 10E-10;
	public static final double UNDERFLOW_factor = 10E10;

	// Initialize all the states
	public static void initialize(ArrayList<Exon> exons)
	{
		ArrayList<Double> avgs = calculateAverages(exons);
		double expected_RPKM = avgs.get(0);
		double expected_SNPs = avgs.get(1);
		States = State.initializeStates(expected_RPKM, expected_SNPs);
	}

	/*
	// given the observations (RPKM and SNPs) and the current state
	// get the next state
	public static State getNextState(int RPKM, double SNPs, State currentState, int lengthOfIntron)
	{
		State maxState = null;
		double maxProbability = 0;
		for (State s: States)
		{
			double prob = State.getTransitionProbability(currentState, s, lengthOfIntron) / State.getEmissionProbability(RPKM, SNPs, s);
			if ( prob > maxProbability )
			{
				maxState = s;
				maxProbability = prob;
			}
		}
		if (Math.random() < maxProbability)
			return maxState;
		else
			return currentState;
	}
	 */

	private static ArrayList<Double> calculateAverages(ArrayList<Exon> exons)
	{
		double rpkm = 0;
		double snps = 0;
		for(Exon e : exons)
		{
			rpkm += e.FPKM;
			snps += e.SNPs;
		}
		ArrayList<Double> returnValues = new ArrayList<Double>();
		returnValues.add(rpkm/exons.size());
		returnValues.add(snps/exons.size());
		return(returnValues);
	}

	public static ArrayList<ArrayList<Double>> computeForwardProbabilities(ArrayList<Exon> exons)
	{
		ArrayList<ArrayList<Double>> ForwardProbability = new ArrayList<ArrayList<Double>>();
		// Initialize so that we're always starting out in normal state
		ArrayList<Double> initialProbabilities = new ArrayList<Double>();
		initialProbabilities.add(1.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		ForwardProbability.add(initialProbabilities);

		for(int exonIndex = 1; exonIndex < exons.size(); exonIndex ++)
		{
			Exon e = exons.get(exonIndex);
			ArrayList<Double> prob = new ArrayList<Double>();
			// Fill in exon's state array
			for(int state=0; state<HiddenMarkovModel.States.size(); state++) 
			{
				double sum = 0;
				State currentState = HiddenMarkovModel.States.get(state);
				for(int prevExonState = 0; prevExonState < HiddenMarkovModel.States.size(); prevExonState++)
				{
					State prevState = HiddenMarkovModel.States.get(prevExonState);
					int lengthOfIntron = e.posLeft - exons.get(exonIndex-1).posRight;
					sum += ForwardProbability.get(exonIndex-1).get(prevExonState) * State.getTransitionProbability(prevState, currentState, lengthOfIntron);
				}
				prob.add(sum * State.getEmissionProbability(e.FPKM, e.SNPs, HiddenMarkovModel.States.get(state)));					
			}

			// Deals with underflow issues
			// as we iterate through the array, the probabilities become increasingly small 
			boolean flag = false;
			for (double p : prob)
			{
				if (p < UNDERFLOW && p > 0.0)
					flag = true;
			}
			if (flag)
			{
				for (int i=0; i<prob.size(); i++)
					prob.set(i, prob.get(i) * UNDERFLOW_factor);
			}

			ForwardProbability.add(prob);
		}	
		return (ForwardProbability);
	}

	public static ArrayList<ArrayList<Double>> computeBackwardProbabilities(ArrayList<Exon> exons)
	{
		ArrayList<ArrayList<Double>> BackwardProbability = new ArrayList<ArrayList<Double>>();
		for (int i=0; i<exons.size(); i++)
			BackwardProbability.add(new ArrayList<Double>());
		// Initialize so that we're always starting out in normal state
		ArrayList<Double> initialProbabilities = new ArrayList<Double>();
		initialProbabilities.add(1.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		initialProbabilities.add(0.0);
		BackwardProbability.set(exons.size()-1, initialProbabilities);

		for(int exonIndex = exons.size() - 2; exonIndex >= 0; exonIndex --)
		{
			Exon e = exons.get(exonIndex);
			ArrayList<Double> prob = new ArrayList<Double>();
			for(int state=0; state<HiddenMarkovModel.States.size(); state++) // for each of the states
			{
				double sum = 0;
				State currentState = HiddenMarkovModel.States.get(state);
				for(int nextExonState = 0; nextExonState < HiddenMarkovModel.States.size(); nextExonState++)
				{
					State nextState = HiddenMarkovModel.States.get(nextExonState);
					int lengthOfIntron = exons.get(exonIndex+1).posLeft - e.posRight;
					sum += BackwardProbability.get(exonIndex+1).get(nextExonState) 
					* State.getTransitionProbability(currentState, nextState, lengthOfIntron)
					* State.getEmissionProbability(exons.get(exonIndex+1).FPKM, exons.get(exonIndex+1).SNPs, nextState);
				}
				prob.add(sum); 				
			}
			
			// Deals with underflow issues
			// as we iterate through the array, the probabilities become increasingly small 
			boolean flag = false;
			for (double p : prob)
			{
				if (p < UNDERFLOW && p > 0.0)
					flag = true;
			}
			if (flag)
			{
				for (int i=0; i<prob.size(); i++)
					prob.set(i, prob.get(i) * UNDERFLOW_factor);
			}

			BackwardProbability.set(exonIndex, prob);
		}

		return (BackwardProbability);
	}

	/*
	 * Given the list of exons (rpkm value & number of SNPs)
	 * compute the most likely state for each using forward-backward algorithm
	 */			
	public static void getStates(ArrayList<Exon> exons)
	{
		HiddenMarkovModel.initialize(exons);

		ArrayList<ArrayList<Double>> ForwardProbability = HiddenMarkovModel.computeForwardProbabilities(exons);
		ArrayList<ArrayList<Double>> BackwardProbability = HiddenMarkovModel.computeBackwardProbabilities(exons);

		// for each exon
		// most likely state = argmax(forward prob * back prob)
		for(int exonIndex = 0; exonIndex < exons.size(); exonIndex ++)
		{
			double maxProb = 0;
			State maxState = null;
			for(int stateIndex = 0; stateIndex < HiddenMarkovModel.States.size(); stateIndex ++)
			{
				double prob = ForwardProbability.get(exonIndex).get(stateIndex) * BackwardProbability.get(exonIndex).get(stateIndex);
				if ( prob > maxProb )
				{
					maxProb = prob;
					maxState = HiddenMarkovModel.States.get(stateIndex);
				}
			}
			Exon e = exons.get(exonIndex);
			e.state = maxState;
			exons.set(exonIndex, e);
		}
	}

}
