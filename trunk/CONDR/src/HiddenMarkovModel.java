
/*
 * Calculates and stores states according to the exon data
 */

import java.util.ArrayList;
import java.util.HashMap;

public class HiddenMarkovModel
{
	static HashMap<String, State> States = new HashMap<String, State>();

	// constant: helps maintain minimum value for probability that we want to store
	public static final double UNDERFLOW = 10E-10;
	public static final double UNDERFLOW_factor = 10E10;

	public static void getStates(ArrayList<Exon> exons, ArrayList<Exon> expectedValues, ArrayList<Exon> stdDeviations, String parameterFileName)
	{
		States = State.initializeStates(parameterFileName);

		ArrayList<HashMap<State, Double>> ForwardProbability = HiddenMarkovModel.computeForwardProbabilities(exons, expectedValues, stdDeviations);
		ArrayList<HashMap<State, Double>> BackwardProbability = HiddenMarkovModel.computeBackwardProbabilities(exons, expectedValues, stdDeviations);

		// for each exon
		// most likely state = argmax(forward prob * back prob)
		for(int exonIndex = 0; exonIndex < exons.size(); exonIndex ++)
		{
			double maxProb = Double.NEGATIVE_INFINITY;
			State maxState = null;
			for(int stateIndex = 0; stateIndex < HiddenMarkovModel.States.size(); stateIndex ++)
			{
				double prob = ForwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex)) 
						* BackwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex)); //stateIndex);
				if ( prob > maxProb )
				{
					maxProb = prob;
					maxState = HiddenMarkovModel.getStateFromIndex(stateIndex); 
				}
			}
			Exon e = exons.get(exonIndex);
			e.state = maxState;
			exons.set(exonIndex, e);
		}
	}

	private static State getStateFromIndex(int stateIndex)
	{
		for(State s : States.values())
		{
			if (s.stateName == stateIndex)
				return s;
		}
		return null;
	}

	private static ArrayList<HashMap<State, Double>> computeBackwardProbabilities(ArrayList<Exon> exons, 
			ArrayList<Exon> expectedValues, ArrayList<Exon> stdDeviations)
			{
		ArrayList<HashMap<State, Double>> BackwardProbability = new ArrayList<HashMap<State, Double>>();

		for (int i=0; i<exons.size(); i++)
			BackwardProbability.add(new HashMap<State, Double>());
		// Initialize so that we're always starting out in normal state
		HashMap<State, Double> initialProbabilities = new HashMap<State, Double>();
		for(State s : HiddenMarkovModel.States.values())
			if (s.stateName == State.NORMAL)
				initialProbabilities.put(s, 1.0);
			else
				initialProbabilities.put(s, 0.0);
		BackwardProbability.set(exons.size()-1, initialProbabilities);

		for(int exonIndex = exons.size() - 2; exonIndex >= 0; exonIndex --)
		{
			Exon e = exons.get(exonIndex);
			Exon expected = expectedValues.get(exonIndex+1);
			Exon stdDev = stdDeviations.get(exonIndex+1);
			HashMap<State, Double> prob = new HashMap<State, Double>();
			
			for( State currentState : HiddenMarkovModel.States.values() )
			{
				double sum = 0;
				for( State nextState : HiddenMarkovModel.States.values()) 
				{
					int lengthOfIntron = exons.get(exonIndex+1).posLeft - e.posRight;
					sum += BackwardProbability.get(exonIndex+1).get(nextState)  
					* State.getTransitionProbability(currentState, nextState, lengthOfIntron, e.length())
					* State.getEmissionProbability(exons.get(exonIndex+1).FPKM, exons.get(exonIndex+1).SNPs, 
							expected.FPKM, expected.SNPs, stdDev.FPKM, stdDev.SNPs, nextState);
				}
				prob.put(currentState, sum); 				
			}

			// Deals with underflow issues
			// as we iterate through the array, the probabilities become increasingly small 
			boolean flag = false;
			for (double p : prob.values())
			{
				if (p < UNDERFLOW && p > 0.0)
					flag = true;
			}
			if (flag)
				for (State s : prob.keySet())
					prob.put(s, prob.get(s) * UNDERFLOW_factor);

			BackwardProbability.set(exonIndex, prob);
		}

		return (BackwardProbability);
			}

	private static ArrayList<HashMap<State, Double>> computeForwardProbabilities(ArrayList<Exon> exons, 
			ArrayList<Exon> expectedValues, ArrayList<Exon> stdDeviations)
			{
		ArrayList<HashMap<State, Double>> ForwardProbability = new ArrayList<HashMap<State, Double>>();
		HashMap<State, Double> initialProbabilities = new HashMap<State, Double>();
		for(State s : HiddenMarkovModel.States.values())
			if (s.stateName == State.NORMAL)
				initialProbabilities.put(s, 1.0);
			else
				initialProbabilities.put(s, 0.0);
		ForwardProbability.add(initialProbabilities);

		for(int exonIndex = 1; exonIndex < exons.size(); exonIndex ++)
		{
			Exon e = exons.get(exonIndex);
			Exon expected = expectedValues.get(exonIndex);
			Exon stdDev = stdDeviations.get(exonIndex);

			HashMap<State, Double> prob = new HashMap<State, Double>();
			// Fill in exon's state array

			for(State currentState : HiddenMarkovModel.States.values()) 
			{
				double sum = 0;
				for(State prevState : HiddenMarkovModel.States.values())
				{
					int lengthOfIntron = e.posLeft - exons.get(exonIndex-1).posRight;
					sum += ForwardProbability.get(exonIndex-1).get(prevState) 
					* State.getTransitionProbability(prevState, currentState, lengthOfIntron, e.length());
				}
				prob.put(currentState, (sum * State.getEmissionProbability(e.FPKM, e.SNPs, expected.FPKM, 
						expected.SNPs, stdDev.FPKM, stdDev.SNPs, currentState)));	
			}

			// Deals with underflow issues
			// as we iterate through the array, the probabilities become increasingly small 
			boolean flag = false;
			for (double p : prob.values())
			{
				if (p < UNDERFLOW && p > 0.0)
					flag = true;
			}
			if (flag)
				for (State s : prob.keySet())
					prob.put(s, prob.get(s) * UNDERFLOW_factor);

			ForwardProbability.add(prob);
		}	
		return (ForwardProbability);
			}

}
