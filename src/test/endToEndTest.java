package test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import jxl.NumberCell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import model.FBA;
import model.FBAreader;

import org.junit.Ignore;
import org.junit.Test;

public class endToEndTest {
	private FBAreader fbaReader;
	
	@Ignore
	@Test
	public void runModel_getsSameResultsAsOriginalCode(){
		fbaReader = new FBAreader("resources/in/original.xls","resources/out.xls");
	
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
		fba.loadNames(fbaReader.compounds);
	
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
	
		F2.loadNames(fbaReader.compounds);
	
		answer = F2.optimise();
	
		System.out.println("Growth is " + answer[fbaReader.noReactions-1]);
	
		for(int k = 0;k < inBiomass.length;k++) {
			fbaReader.biomassIn[k] = inBiomass[k];		   
		}
	
		double[] actualFluxes = getCalculatedFluxes("resources/ref/original.xls");
		
		for(int i=0;i<actualFluxes.length; i++){
			assertTrue("Flux for " + fbaReader.reactions[i].name + " doesn't match: Expected <"
					+ actualFluxes[i] + ">, Actual <" + answer[i] + ">",actualFluxes[i] == answer[i]);
		}
	}
	
	@Test
	public void runModel_getsSameResultsAsDivergence1(){
		fbaReader = new FBAreader("resources/in/divergence1.xls","resources/out.xls");
	
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
		fba.loadNames(fbaReader.compounds);
	
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
	
		F2.loadNames(fbaReader.compounds);
	
		answer = F2.optimise();
	
		System.out.println("Growth is " + answer[fbaReader.noReactions-1]);
	
		for(int k = 0;k < inBiomass.length;k++) {
			fbaReader.biomassIn[k] = inBiomass[k];		   
		}
	
		double[] actualFluxes = getCalculatedFluxes("resources/ref/divergence1.xls");
		
		for(int i=0;i<actualFluxes.length; i++){
			assertTrue("Flux for " + fbaReader.reactions[i].name + " doesn't match: Expected <"
					+ actualFluxes[i] + ">, Actual <" + answer[i] + ">",actualFluxes[i] == answer[i]);
		}
	}

	private double[] getCalculatedFluxes(String filename) {
		Workbook input;
		try {
			input = Workbook.getWorkbook(new File(filename));
			Sheet reactionSheet = input.getSheet(1);
			int noReactions = reactionSheet.getRows();
			double[] fluxes = new double[noReactions];
			
			for(int j = 0;j < noReactions;j++) {
				NumberCell cell = (NumberCell) reactionSheet.getCell(5,j);				
	        	fluxes[j] = Double.valueOf(cell.getValue());
	        }
			return fluxes;
		} catch (BiffException | IOException e) {
			e.printStackTrace();
			return null;
		}
	}
}
