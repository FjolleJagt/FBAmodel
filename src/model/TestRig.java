package model;

import java.io.IOException;
import java.util.ArrayList;

public class TestRig {
	private FBAreader R;
	
	public void test(){
		R = new FBAreader("resources/in.xls","resources/out.xls");
	
		R.createSmatrix();
		boolean [] inBiomass = new boolean[R.noBiomass];
	
		ArrayList<Integer> onesNotIn = new ArrayList<Integer>();
		for(int k = 0;k < inBiomass.length;k++) {
			inBiomass[k] = R.biomassIn[k];
			if(!inBiomass[k]) {
				onesNotIn.add(k);
			}
		}
	
		double [] answer = new double[R.noReactions];
	
		R.createBiomassReaction(inBiomass);
	
		FBA F = new FBA(R.noReactions,R.noCompounds);
	
		F.loadS(R.S);			
		F.loadVectors(R.a,R.b,R.c,R.lb,R.ub);			
		F.loadNames(R.compoundNames,R.reactionNames);
	
		answer = F.optimise();
	
		//Now remove anything futile -- secondary minimisation
	
		//Lock the growth rate
		R.a[R.noReactions-1] =  answer[R.noReactions-1];
		R.b[R.noReactions-1] =  answer[R.noReactions-1];
		R.c[R.noReactions-1] = 0;
		R.ub[R.noReactions-1] = true;
		R.lb[R.noReactions-1] = true;
	
		//Minimise everthing else
		for(int k = 0;k < R.noReactions - 2;k++) {
			if(answer[k] > 0) {
				if(R.b[k] == 1000) {
					R.c[k] = -1;		
					if(R.a[k] < 0) {
						R.a[k] = 0;
					}
					R.lb[k] = true;
				}	
			}
			else {
				if(R.a[k]==-1000) {
					R.c[k] = 1;	
					if(R.b[k] > 0) {
						R.b[k] = 0;
					}
					R.ub[k] = true;
				}
			}
		}		
	
		FBA F2 = new FBA(R.noReactions,R.noCompounds);
	
		F2.loadS(R.S);           
		F2.loadVectors(R.a,R.b,R.c,R.lb,R.ub);   
	
		F2.loadNames(R.compoundNames,R.reactionNames);
	
		answer = F2.optimise();
	
		//To here
	
		System.out.println("Growth is " + answer[R.noReactions-1]);
	
		for(int k = 0;k < inBiomass.length;k++) {
			R.biomassIn[k] = inBiomass[k];		   
		}
	
		try {		    
			R.writeSmatrix(answer);
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
