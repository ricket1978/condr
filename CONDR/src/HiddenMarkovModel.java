
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

	public static void initialize(String parameterFileName)
	{
		States = State.initializeStates(parameterFileName);		
	}

	public static void getStates(ArrayList<Exon> exons, ArrayList<Exon> expectedValues, ArrayList<Exon> stdDeviations, String parameterFileName, double threshold)
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
			//double threshold = 1000.0;

			State normalState = HiddenMarkovModel.States.get("NORMAL");
			double normalProb = ForwardProbability.get(exonIndex).get(normalState) 
			+ BackwardProbability.get(exonIndex).get(normalState);
			maxProb = normalProb;
			maxState = normalState;
			//System.out.println(exons.get(exons.size()-1));
			for(int stateIndex = 0; stateIndex < HiddenMarkovModel.States.size(); stateIndex ++)
			{
				double prob = ForwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex)) 
				+ BackwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex));
				if (exons.get(exonIndex).posLeft==15258731 ||
						exons.get(exonIndex).posLeft==33846423 ||
						exons.get(exonIndex).posLeft==37901635 ||
						exons.get(exonIndex).posLeft==39171084 ||
						exons.get(exonIndex).posLeft==42280008 ||
						exons.get(exonIndex).posLeft==42284072)
					System.out.println(exons.get(exonIndex).posLeft+"\t"+stateIndex + "\t" + ForwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex)) 
							+"\t"+ BackwardProbability.get(exonIndex).get(HiddenMarkovModel.getStateFromIndex(stateIndex)) + "\t" + prob);

				if ( (prob > maxProb) && (prob > (normalProb + threshold)) )
				{
					maxProb = prob;
					maxState = HiddenMarkovModel.getStateFromIndex(stateIndex); 
				}
			}
			Exon e = exons.get(exonIndex);
			e.state = maxState;
			if (exons.get(exonIndex).posLeft==15258731 ||
					exons.get(exonIndex).posLeft==33846423 ||
					exons.get(exonIndex).posLeft==37901635 ||
					exons.get(exonIndex).posLeft==39171084 ||
					exons.get(exonIndex).posLeft==42280008 ||
					exons.get(exonIndex).posLeft==42284072)
				System.out.println(maxState);
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
				initialProbabilities.put(s, 0.0);
			else
				initialProbabilities.put(s, Double.NEGATIVE_INFINITY);
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
				if (exons.get(exonIndex).posLeft == 37901635)
					System.out.println(currentState);
				for( State nextState : HiddenMarkovModel.States.values()) 
				{
					int lengthOfIntron = exons.get(exonIndex+1).posLeft - e.posRight;

					double b_prime = BackwardProbability.get(exonIndex+1).get(nextState);
					double logTransProb = Math.log(State.getTransitionProbability(currentState, nextState, lengthOfIntron, e.length()));
					double e_prime = State.getLogEmissionProbability(exons.get(exonIndex+1), 
							expected, stdDev, nextState);
					double p_prime = b_prime + logTransProb + e_prime;
					//if (exonIndex == exons.size()-2)
					//System.out.println("!!!!!"+lengthOfIntron + "\t" + currentState.stateName + "\t" + nextState.stateName + "\t" + sum + "\t" + b_prime + "\t" + e_prime + "\t" + logTransProb + "\t" + p_prime);
					//if (exons.get(exonIndex).posLeft == 37901635 && (currentState.stateName==0 || currentState.stateName == 4))
					//System.out.println(currentState.stateName + "\t" + nextState.stateName + "\t" + b_prime + "\t" + e_prime 
					//	+ "\t" + logTransProb + "\t" + p_prime);
					if (sum == 0) // initial condition.. otherwise sum is translated as starting out with a probability of 1
						sum = p_prime;
					else
						sum = Probability.logSum(sum, p_prime);
					//sum = Probability.logSum(sum, p_prime);
					if (sum > 0)
						System.out.println("--- No log(prob) should be > 0");						
				}
				//if (exons.get(exonIndex).posLeft == 37901635 && (currentState.stateName==0 || currentState.snpRatio == 4))
				//System.out.println(currentState.stateName + "\t" + sum);
				if (Double.isNaN(sum)) // deals with the case that probability = 0 and hence log probability is undefined
					prob.put(currentState, Double.NEGATIVE_INFINITY);
				else
					prob.put(currentState, sum); 		
				if (sum > 0)
					System.out.println("Bwd: No log(prob) should be > 0");					
				if (sum==0)
					System.out.println("Bwd: Sum == 0");
			}
			//TODO case of 0,0 is not 3 when those are the expected

			//if (exonIndex == 106)
			//System.out.println("Weird exon");
			// if all prob are -Inf, then its useless
			boolean allInf = true;
			for(Double p : prob.values())
				if (!Double.isInfinite(p))
					allInf = false;
			if (allInf)
				System.out.println(e + "\tAll backward prob are -Inf");

			BackwardProbability.set(exonIndex, prob);

		}

		return (BackwardProbability);
			}

	private static ArrayList<HashMap<State, Double>> computeForwardProbabilities(ArrayList<Exon> exons, ArrayList<Exon> expectedValues, ArrayList<Exon> stdDeviations)
	{
		ArrayList<HashMap<State, Double>> ForwardProbability = new ArrayList<HashMap<State, Double>>();
		HashMap<State, Double> initialProbabilities = new HashMap<State, Double>();
		for(State s : HiddenMarkovModel.States.values())
			if (s.stateName == State.NORMAL)
				initialProbabilities.put(s, 0.0);
			else
				initialProbabilities.put(s, Double.NEGATIVE_INFINITY);
		ForwardProbability.add(initialProbabilities);

		for(int exonIndex = 1; exonIndex < exons.size(); exonIndex ++)
		{
			Exon e = exons.get(exonIndex);
			Exon expected = expectedValues.get(exonIndex);
			Exon stdDev = stdDeviations.get(exonIndex);

			//if (exonIndex == 106)
			//System.out.println("$$");
			HashMap<State, Double> prob = new HashMap<State, Double>();
			// Fill in exon's state array

			for(State currentState : HiddenMarkovModel.States.values()) 
			{
				int lengthOfIntron = e.posLeft - exons.get(exonIndex-1).posRight;
				// TODO: hack for normal numbers. need to figure out this case:
				/*
				 * chr16	48820563	48820609	TRF4-2
					chr16	48820563	48826720	PAPD5
				 */
				lengthOfIntron = Math.abs(lengthOfIntron);
				Double sum = null;

				for(State prevState : HiddenMarkovModel.States.values())
				{
					/*
					 * log(p + q) = p_prime + log( 1 + exp(q_prime - p_prime) )
					 * p_prime = log p
					 * q_prime = log q
					 * 
					 * p = f*transition prob
					 * p_prime = f_prime + log transition prob
					 */
					double f_prime = ForwardProbability.get(exonIndex-1).get(prevState);
					double logTransProb = Math.log(State.getTransitionProbability(prevState, currentState, lengthOfIntron, e.length()));
					double p_prime = f_prime + logTransProb;
					//Double prevSum = sum;
					//if (sum == 0)
					if (sum==null)
						sum = p_prime;
					else
						sum = Probability.logSum(sum, p_prime);
					//if (e.posLeft==15117160 && currentState.stateName == 2)
						//System.out.println(e.posLeft + "\t" + currentState.stateName + "\t" + prevState.stateName + "\t" +  
							//	f_prime + "\t" + logTransProb + "\t" + p_prime + "\t" + prevSum + "\t" + sum);
					//firstTime = false;
				}
				if (sum > 0)
					System.out.println("Fwd: No log(prob) should be > 0");
				if (sum==0)
					System.out.println("Fwd: Sum == 0");

				//if (e.posLeft == 10033308)
				//System.out.println(currentState + "\t" + sum + "\t" + State.getLogEmissionProbability(e, expected, stdDev, currentState)
				//	+ "\tSum: " + (sum + State.getLogEmissionProbability(e, expected, stdDev, currentState)));
				if (Double.isNaN(sum + State.getLogEmissionProbability(e, expected, stdDev, currentState))) // deals with the case that probability = 0 and hence log probability is undefined
					prob.put(currentState, Double.NEGATIVE_INFINITY);
				else
					prob.put(currentState, (sum + State.getLogEmissionProbability(e, expected, stdDev, currentState)));
			}

			//if (exonIndex == 56)
			//System.out.println("Weird exon");

			// if all prob are -Inf, then its useless
			boolean allInf = true;
			for(Double p : prob.values())
				if (!Double.isInfinite(p))
				{allInf = false;break;}
			if (allInf)
				System.out.println(e + "\tAll forward prob are -Inf");
			ForwardProbability.add(prob);
			//System.out.println("-----");
		}	
		return (ForwardProbability);
	}

}
