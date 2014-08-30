package model;

import java.io.IOException;
import java.util.ArrayList;

public class TestRig {
	private FBAreader fbaReader;
	
	public void test(){
		fbaReader = new FBAreader("resources/in.xls","resources/out.csv");
	
		fbaReader.createSmatrix();
		boolean [] inBiomass = new boolean[fbaReader.noBiomass];
	
		ArrayList<Integer> onesNotIn = new ArrayList<Integer>();
		for(int k = 0;k < inBiomass.length;k++) {
			inBiomass[k] = fbaReader.biomassIn[k];
			if(!inBiomass[k]) {
				onesNotIn.add(k);
			}
		}
	
		double [] answer = new double[fbaReader.noReactions];
	
		fbaReader.createBiomassReaction(inBiomass);
	
		FBA fba = new FBA(fbaReader.noReactions,fbaReader.noCompounds);
	
		fba.loadS(fbaReader.S);			
		fba.loadReactions(fbaReader.reactions);			
		fba.loadNames(fbaReader.compoundNames);
	
		answer = fba.optimise();
	
		//Now remove anything futile -- secondary minimisation
	
		//Lock the growth rate
		fbaReader.reactions[fbaReader.noReactions-1].lowerBound =  answer[fbaReader.noReactions-1];
		fbaReader.reactions[fbaReader.noReactions-1].upperBound =  answer[fbaReader.noReactions-1];
		fbaReader.reactions[fbaReader.noReactions-1].optimisationCoefficient = 0.0;
	
		//Minimise everything else
		for(int k = 0;k < fbaReader.noReactions - 2;k++) {
			if(answer[k] > 0) {
				if(fbaReader.reactions[k].upperBound.equals(Double.valueOf(1000))) {
					fbaReader.reactions[k].optimisationCoefficient = -1.0;		
					if(fbaReader.reactions[k].lowerBound.doubleValue() < 0) {
						fbaReader.reactions[k].lowerBound = 0.0;
					}
				}	
			} else {
				if(fbaReader.reactions[k].lowerBound.equals(Double.valueOf(-1000))) {
					fbaReader.reactions[k].optimisationCoefficient = 1.0;	
					if(fbaReader.reactions[k].upperBound.doubleValue() > 0) {
						fbaReader.reactions[k].upperBound = 0.0;
					}
				}
			}
		}		
	
		FBA F2 = new FBA(fbaReader.noReactions,fbaReader.noCompounds);
	
		F2.loadS(fbaReader.S);           
		F2.loadReactions(fbaReader.reactions);   
	
		F2.loadNames(fbaReader.compoundNames);
	
		answer = F2.optimise();
	
		//To here
	
		System.out.println("Growth is " + answer[fbaReader.noReactions-1]);
	
		for(int k = 0;k < inBiomass.length;k++) {
			fbaReader.biomassIn[k] = inBiomass[k];		   
		}
	
		try {		    
			fbaReader.writeSmatrix(answer);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	
		/*int z = ((Number) lookSpin.getValue()).intValue() - 1;
	
	
		System.out.println("Examining reactant no " + z + " which is " + R.compoundNames[z]);
		for(int k = 0;k < answer.length;k++) {
			if(R.S[z][k] != 0) {
				if(Math.abs(answer[k]) > 0.01) {
					System.out.println("ACTIVE: " + R.reactionNames[k] + "," + R.reactions[k]+ "," + answer[k]);
				}
				else {
					System.out.println("INACTIVE: " + R.reactionNames[k] + "," + R.reactions[k]+ "," + answer[k]);
				}
			}	 
		} 
	
		System.out.println("And done");*/
	}
}
