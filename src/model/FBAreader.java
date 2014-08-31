package model;

import java.io.*;

import jxl.*;
import jxl.read.biff.BiffException;
import jxl.write.*;

public class FBAreader {
    String inputFileName;
    String outputFileName;

    Workbook input;

    Sheet compoundSheet;
    Sheet biomassSheet;

    public int noReactions;
    public int noCompounds;
    public int noBiomass;

    public Compound [] compounds;

    public Reaction [] reactions;
    public Compound [] biomassNames;
    double [] biomassComp;
    public boolean [] biomassIn;
    int [] biomassNo;

    boolean [] usedCompounds;

    public double [][] S;

    public FBAreader(String in, String out) {
        inputFileName = in;
        outputFileName = out;

        try {
            read();
        }
        catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void read() throws IOException {
        try {
            input = Workbook.getWorkbook(new File(inputFileName));
        }
        catch (IOException | BiffException m) {
            m.printStackTrace();
        }

        compoundSheet = input.getSheet(0);
        biomassSheet = input.getSheet(2);

        noCompounds = compoundSheet.getRows();
        noBiomass = biomassSheet.getRows();

        compounds = new Compound [noCompounds+1];
        usedCompounds = new boolean [noCompounds+1];
        

        biomassNames = new Compound [noBiomass];
        biomassComp = new double [noBiomass];
        biomassIn = new boolean [noBiomass];

        for(int i = 0;i < noCompounds;i++) {
            usedCompounds[i] = false;
            compounds[i] = new Compound(compoundSheet.getCell(0,i).getContents());
        }

        usedCompounds[noCompounds]=false;

        loadReactions();

        for(int k = 0;k < noBiomass;k++) {
            biomassNames[k] = new Compound(biomassSheet.getCell(0,k).getContents());
            Cell nc1 = biomassSheet.getCell(1,k);
            NumberCell nc = (NumberCell) nc1;
            biomassComp[k] = Double.valueOf(nc.getValue());
            biomassIn[k] = false;
            if(Double.valueOf(biomassSheet.getCell(2,k).getContents())>0) {
                biomassIn[k] = true;
            }
        }

        biomassNo = new int [noBiomass];

        for(int i = 0;i < noCompounds;i++) {
            for(int k = 0;k < noBiomass;k++) {
                if(biomassNames[k].equals(compounds[i])) {
                    biomassNo[k] = i;
                }
            }
        }

        input.close();
        System.out.println(noCompounds +" compounds, " + noReactions + " reactions loaded");

    }

	private void loadReactions() {
		Sheet reactionSheet = input.getSheet(1);
		
		noReactions = reactionSheet.getRows();
		reactions = new Reaction[noReactions+2];
		
		for(int j = 0;j < noReactions;j++) {
        	reactions[j] = new Reaction();
            reactions[j].name = reactionSheet.getCell(0,j).getContents();
            reactions[j].equation = reactionSheet.getCell(1,j).getContents();
            reactions[j].lowerBound = getDoubleFromCell(reactionSheet.getCell(2,j));
            reactions[j].upperBound = getDoubleFromCell(reactionSheet.getCell(3,j));
            reactions[j].optimisationCoefficient = getDoubleFromCell(reactionSheet.getCell(4,j));
        }
	}

	private Double getDoubleFromCell(Cell cell) {
		if(cell.getType() == CellType.LABEL) {
		    String contents = cell.getContents();
		    if(contents.compareTo("infty") == 0 || contents.compareTo("-infty") == 0){
		        return null;
		    } else {
		    	throw new RuntimeException("Couldn't parse cell value as double: " + contents);
		    }
		} else {
		    NumberCell nc = (NumberCell) cell;
		    return new Double(nc.getValue());
		}
	}


    public void createSmatrix() {
        S = new double [noCompounds+1][noReactions+2];
        //+2, +1 for the biomass reaction, growth and biomass to be included later.

        for(int j = 0;j < noReactions;j++) {
            double[] stoichiometry = getStoichiometryFromEquation(reactions[j].equation);
            for(int i = 0; i < noCompounds; i++){
            	S[i][j] += stoichiometry[i];
            }
        }
        for(int i = 0;i < noCompounds;i++) {
            if(!usedCompounds[i]) {
                System.out.println("Compound no " + i + ", which is " + compounds[i].name + " has not been used");
            }
        }
    }

	public double[] getStoichiometryFromEquation(String equation) {
		String [] parts;

		//pad equation, else split("foo<==>") will return ["foo"] not ["foo",""]
		equation = equation + " ";
		
		if(equation.contains("-->")) {
		    parts = equation.split("-->");
		} else if(equation.contains("<==>")) {
		    parts = equation.split("<==>");
		} else {
			throw new RuntimeException("No arrow or unknown arrow type in equation: " + equation);
		}
		double[] leftStoichiometry = getStoichiometryFromSideOfEquation(parts[0]);
		double[] rightStoichiometry = getStoichiometryFromSideOfEquation(parts[1]);
		double[] stoichiometry = new double[noReactions];
		for(int i=0;i<noReactions;i++){
			stoichiometry[i] = rightStoichiometry[i] - leftStoichiometry[i];
		}
		return stoichiometry;
	}

	public double[] getStoichiometryFromSideOfEquation(String equationSide) {
		double[] stoichiometry = new double[noReactions];
		String [] compoundNames = splitIntoCompounds(equationSide);
		double [] compoundStoichiometry = new double [compoundNames.length];

		//check for stoichiometry different from 1
		for(int k = 0; k < compoundNames.length;k++) {
		    if(compoundNames[k].contains("*")) {
		    	String [] chop = compoundNames[k].split("\\*");
		        compoundNames[k] = chop[1].trim();
		        compoundStoichiometry[k] = Double.valueOf(chop[0].trim());
		    } else {
		        compoundStoichiometry[k] = 1.0;
		        compoundNames[k] = compoundNames[k].trim();
		    }

		    boolean nothingFound = true;

		    //now find them...
		    for(int i = 0;i < noCompounds;i++) {
		        if(compoundNames[k].equals(compounds[i].name)) {
		            stoichiometry[i] = compoundStoichiometry[k];
		            nothingFound = false;
		            usedCompounds[i] = true;
		        } else if(compoundNames[k].equals("")) {
	                nothingFound = false;
		        }
		    }
		    
		    if(nothingFound){
		    	throw new RuntimeException("Couldn't find compound '" + compoundNames[k] + "' in list of compound names");
		    }
		}
		return stoichiometry;
	}

	private String[] splitIntoCompounds(String equationSide) {
		String[] parts;
		if(equationSide.contains("+")) {
		    parts = equationSide.split("\\+");
		} else {
		    parts = new String [1];
		    parts[0] = equationSide;
		}
		return parts;
	}
    public void createBiomassReaction(boolean [] in) {

        noReactions++; //one for the biomass reaction
        noReactions++; //one for growth
        
        reactions[noReactions-2] = new Reaction();
        reactions[noReactions-1] = new Reaction();

        reactions[noReactions-2].name = "Biomass reaction";
        reactions[noReactions-1].name = "Growth";

        String pos="",neg="";
        boolean firstN = true;
        boolean firstP = true;

        for(int k = 0;k < noBiomass;k++) {
            if(in[k]) {
                System.out.println(compounds[biomassNo[k]].name + " is in with value " + biomassComp[k]);

                S[biomassNo[k]][noReactions-2] = biomassComp[k];

                if(Double.valueOf(biomassComp[k])< 0) {
                    if(firstN) {
                        neg+=String.valueOf(-1*biomassComp[k]) + "*" + compounds[biomassNo[k]].name;
                        firstN = false;
                    }
                    else {
                        neg+="+" + String.valueOf(-1*biomassComp[k]) + "*" + compounds[biomassNo[k]].name;
                    }
                }
                else if(Double.valueOf(biomassComp[k]) > 0) {
                    if(firstP) {
                        pos+=String.valueOf(biomassComp[k]) + "*" + compounds[biomassNo[k]].name;
                        firstP = false;
                    }
                    else {
                        pos+="+" + String.valueOf(biomassComp[k]) + "*" + compounds[biomassNo[k]].name;
                    }
                }
            }
            else {
                S[biomassNo[k]][noReactions-2] = 0.0;
            }
        }

        reactions[noReactions-2].equation = neg + "-->" + pos;
        reactions[noReactions-1].equation = "Biomass-->";


        noCompounds++; //one for the biomass

        compounds[noCompounds-1] = new Compound("Biomass");

        S[noCompounds-1][noReactions-2] = 1;
        S[noCompounds-1][noReactions-1] = -1;

        reactions[noReactions-2].lowerBound = 0.0;
        reactions[noReactions-2].upperBound = null;
        reactions[noReactions-2].optimisationCoefficient = 0.0;

        reactions[noReactions-1].lowerBound = 0.0;
        reactions[noReactions-1].upperBound = null;
        reactions[noReactions-1].optimisationCoefficient = 1.0;
    }

    public void writeSmatrix(double [] answer) throws IOException {
        try {
            input = Workbook.getWorkbook(new File(inputFileName));

        }
        catch (jxl.read.biff.BiffException m) {
            m.printStackTrace();
        }


        try {

            WritableWorkbook copy = Workbook.createWorkbook(new File(outputFileName), input);

            WritableSheet sheet2 = copy.getSheet(1);
            WritableSheet sheet3 = copy.getSheet(2);
            jxl.write.Label label;
            jxl.write.Number number;

            label = new jxl.write.Label(0,noReactions-2, reactions[noReactions-2].name);
            sheet2.addCell(label);
            label = new jxl.write.Label(1,noReactions-2, reactions[noReactions-2].equation);
            sheet2.addCell(label);
            if(reactions[noReactions-2].lowerBound != null){
            	number = new jxl.write.Number(2,noReactions-2, reactions[noReactions-2].lowerBound);
            	sheet2.addCell(number);
            } else {
            	sheet2.addCell(new jxl.write.Number(2,noReactions-2,0.0));
            }
            if(reactions[noReactions-2].upperBound != null){
            	number = new jxl.write.Number(3,noReactions-2, reactions[noReactions-2].upperBound);
            	sheet2.addCell(number);
            } else {
            	sheet2.addCell(new jxl.write.Number(3,noReactions-2,999));
            }
            number = new jxl.write.Number(4,noReactions-2, reactions[noReactions-2].optimisationCoefficient);
            sheet2.addCell(number);
            label = new jxl.write.Label(0,noReactions-1, reactions[noReactions-1].name);
            sheet2.addCell(label);
            label = new jxl.write.Label(1,noReactions-1, reactions[noReactions-1].equation);
            sheet2.addCell(label);
            if(reactions[noReactions-1].lowerBound != null){
            	number = new jxl.write.Number(2,noReactions-1, reactions[noReactions-1].lowerBound);
            	sheet2.addCell(number);
            }
            if(reactions[noReactions-1].upperBound != null){
            	number = new jxl.write.Number(3,noReactions-1, reactions[noReactions-1].upperBound);
            	sheet2.addCell(number);
            }
            number = new jxl.write.Number(4,noReactions-1, reactions[noReactions-1].optimisationCoefficient);
            sheet2.addCell(number);

            for(int j = 0;j < noReactions;j++) {
                number = new jxl.write.Number(5,j, answer[j]);
                sheet2.addCell(number);
            }

            WritableCellFormat goodCF = new WritableCellFormat();
            goodCF.setWrap(true);
            goodCF.setBackground(Colour.GREEN);
            WritableCellFormat badCF = new WritableCellFormat();
            badCF.setWrap(true);
            badCF.setBackground(Colour.RED);


            for(int k = 0;k < noBiomass;k++) {
                if(biomassIn[k]) {

                    number = new jxl.write.Number(3,k,1,goodCF);

                }
                else {
                    number = new jxl.write.Number(3,k,0,badCF);

                }
                sheet3.addCell(number);
            }


            copy.write();
            copy.close();
            input.close();
        }
        catch (jxl.write.WriteException m) {
            m.printStackTrace();
        }
    }

    public static void main(String[] args) {

        FBAreader fbaReader = new FBAreader(args[0],args[1]);

        fbaReader.createSmatrix();
        FBA fba = new FBA(fbaReader.noReactions,fbaReader.noCompounds);

        fba.loadS(fbaReader.S);

        fba.loadReactions(fbaReader.reactions);
        fba.loadNames(fbaReader.compounds);


        double [] answer = new double[fbaReader.noReactions];
        answer = fba.optimise();
        try {
            fbaReader.writeSmatrix(answer);
        }
        catch (IOException e) {
            e.printStackTrace();
        }


        PrintWriter pwout;
        FileWriter out;

        try {
            out = new FileWriter("sol.csv");
            pwout = new PrintWriter(out);

            for(int i = 0;i < fbaReader.noReactions;i++) {
                pwout.write(fbaReader.reactions[i].name + "," + fbaReader.reactions[i]+ "," + answer[i] + "\n");
            }

            pwout.close();
            try {
                out.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }


    }


}

