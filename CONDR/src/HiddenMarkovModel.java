package CONDR.src;

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
		States = State.initializeStates();
	}

	public static void getStates(ArrayList<Exon> exons, ArrayList<Exon> expectedValues)
	{
		HiddenMarkovModel.initialize(exons);

		ArrayList<ArrayList<Double>> ForwardProbability = HiddenMarkovModel.computeForwardProbabilities(exons, expectedValues);
		ArrayList<ArrayList<Double>> BackwardProbability = HiddenMarkovModel.computeBackwardProbabilities(exons, expectedValues);

		// for each exon
		// most likely state = argmax(forward prob * back prob)
		for(int exonIndex = 0; exonIndex < exons.size(); exonIndex ++)
		{
			double maxProb = 0;
			State maxState = null;
			for(int stateIndex = 0; stateIndex < HiddenMarkovModel.States.size(); stateIndex ++)
			{
				double prob = ForwardProbability.get(exonIndex).get(stateIndex) * BackwardProbability.get(exonIndex).get(stateIndex);
				//System.out.println(ForwardProbability.get(exonIndex).get(stateIndex) + "\t" + 
					//	BackwardProbability.get(exonIndex).get(stateIndex) + "\t" + prob + "\t" + stateIndex);
				if ( prob > maxProb )
				{
					maxProb = prob;
					maxState = HiddenMarkovModel.States.get(stateIndex);
					//System.out.println(stateIndex);
				}
			}
			Exon e = exons.get(exonIndex);
			e.state = maxState;
			if (maxState == null)
				System.out.println("!!!");
			exons.set(exonIndex, e);
		}
	}

	private static ArrayList<ArrayList<Double>> computeBackwardProbabilities(ArrayList<Exon> exons, ArrayList<Exon> expectedValues)
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
			Exon expected = expectedValues.get(exonIndex);
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
					* State.getTransitionProbability(currentState, nextState, lengthOfIntron, e.length())
					* State.getEmissionProbability(exons.get(exonIndex+1).FPKM, exons.get(exonIndex+1).SNPs, expected.FPKM, expected.SNPs, nextState);
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

	private static ArrayList<ArrayList<Double>> computeForwardProbabilities(ArrayList<Exon> exons, ArrayList<Exon> expectedValues)
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
			Exon expected = expectedValues.get(exonIndex);

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
					sum += ForwardProbability.get(exonIndex-1).get(prevExonState) 
					* State.getTransitionProbability(prevState, currentState, lengthOfIntron, e.length());
				}
				prob.add(sum * State.getEmissionProbability(e.FPKM, e.SNPs, expected.FPKM, expected.SNPs, HiddenMarkovModel.States.get(state)));					
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

}
