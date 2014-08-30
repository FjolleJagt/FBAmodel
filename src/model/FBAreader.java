package model;

import java.io.*;

import jxl.*;
import jxl.write.*;

public class FBAreader {
    String inputFileName;
    String outputFileName;

    Workbook input;

    Sheet compoundSheet;
    Sheet biomassSheet;

    int noReactions;
    int noCompounds;
    int noBiomass;

    String [] compoundNames;

    Reaction [] reactions;
    String [] biomassNames;
    double [] biomassComp;
    boolean [] biomassIn;
    int [] biomassNo;

    boolean [] usedCompounds;

    double [][] S;

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
        catch (jxl.read.biff.BiffException m) {
            m.printStackTrace();
        }

        compoundSheet = input.getSheet(0);
        biomassSheet = input.getSheet(2);

        noCompounds = compoundSheet.getRows();
        noBiomass = biomassSheet.getRows();

        compoundNames = new String [noCompounds+1];
        usedCompounds = new boolean [noCompounds+1];
        

        biomassNames = new String [noBiomass];
        biomassComp = new double [noBiomass];
        biomassIn = new boolean [noBiomass];

        for(int i = 0;i < noCompounds;i++) {
            usedCompounds[i] = false;
            compoundNames[i] = compoundSheet.getCell(0,i).getContents();
        }

        usedCompounds[noCompounds]=false;

        loadReactions();

        for(int k = 0;k < noBiomass;k++) {
            biomassNames[k] = biomassSheet.getCell(0,k).getContents();
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
                if(biomassNames[k].equals(compoundNames[i])) {
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

        //this should work, but no error checking, woe betide!!!


        for(int j = 0;j < noReactions;j++) {
            //pad to make sure bits splits it
            reactions[j].equation = reactions[j].equation + " ";

            //split the reactions in the two bits either side of the operator
            //note reversibility as we go
            String [] bits = new String [2];


            if(reactions[j].equation.lastIndexOf("<==>") < 0) {
                bits = reactions[j].equation.split("-->");
            }
            else {
                bits = reactions[j].equation.split("<==>");

                //This is causing the problem, its overwriting what we need...
		/*
		if(a[j] == 0) {		    
		    lb[j] = false;
		}
		*/
            }


            //then break up the parts in the neg side
            String [] negParts;

            if(bits[0].lastIndexOf("+") > 0) {
                negParts = bits[0].split("\\+");
            }
            else {
                negParts = new String [1];
                negParts[0] = bits[0];
            }
            double [] negStoich = new double [negParts.length];

            //check for stoichiometry different from 1
            for(int k = 0;k < negParts.length;k++) {

                String [] chop;
                if(bits[0].lastIndexOf("*") > 0) {
                    chop = negParts[k].split("\\*");
                }
                else {
                    chop = new String [1];
                    chop[0] = negParts[k];
                }
                if(chop.length > 1) {
                    negParts[k] = chop[1].trim();
                    negStoich[k] = Double.valueOf(chop[0].trim());
                }
                else {
                    negParts[k] = chop[0].trim();
                    negStoich[k] = 1.0;
                }

                boolean nothingFound = true;

                //now find them...
                for(int i = 0;i < noCompounds;i++) {
                    if(negParts[k].equals(compoundNames[i])) {
                        S[i][j] = -1*negStoich[k];
                        nothingFound = false;
                        usedCompounds[i] = true;
                    }
                    else {
                        if(negParts[k].equals("")) {
                            nothingFound = false;
                        }
                    }
                }
            }

            //then do the whole thing again for the positive side
            String [] posParts;
            if(bits[1].lastIndexOf("+") > 0) {
                posParts = bits[1].split("\\+");
            }
            else {
                posParts = new String [1];
                posParts[0] = bits[1];
            }
            double [] posStoich = new double [posParts.length];

            //check for stoichiometry different from 1
            for(int k = 0;k < posParts.length;k++) {
                String [] chop;
                if(bits[1].lastIndexOf("*") > 0) {
                    chop = posParts[k].split("\\*");
                }
                else {
                    chop = new String [1];
                    chop[0] = posParts[k];
                }
                if(chop.length > 1) {
                    posParts[k] = chop[1].trim();
                    posStoich[k] = Double.valueOf(chop[0].trim());
                }
                else {
                    posParts[k] = chop[0].trim();
                    posStoich[k] = 1.0;
                }

                boolean nothingFound = true;

                //now find them...
                for(int i = 0;i < noCompounds;i++) {
                    if(posParts[k].equals(compoundNames[i])) {
                        S[i][j] = posStoich[k];
                        nothingFound = false;
                        usedCompounds[i] = true;
                    }
                    else {
                        if(posParts[k].equals("")) {
                            nothingFound = false;
                        }
                    }

                }
            }

        }
        for(int i = 0;i < noCompounds;i++) {
            if(!usedCompounds[i]) {
                System.out.println("Compound no " + i + ", which is " + compoundNames[i] + " has not been used");
            }
        }




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
                System.out.println(compoundNames[biomassNo[k]] + " is in with value " + biomassComp[k]);

                S[biomassNo[k]][noReactions-2] = biomassComp[k];

                if(Double.valueOf(biomassComp[k])< 0) {
                    if(firstN) {
                        neg+=String.valueOf(-1*biomassComp[k]) + "*" + compoundNames[biomassNo[k]];
                        firstN = false;
                    }
                    else {
                        neg+="+" + String.valueOf(-1*biomassComp[k]) + "*" + compoundNames[biomassNo[k]];
                    }
                }
                else if(Double.valueOf(biomassComp[k]) > 0) {
                    if(firstP) {
                        pos+=String.valueOf(biomassComp[k]) + "*" + compoundNames[biomassNo[k]];
                        firstP = false;
                    }
                    else {
                        pos+="+" + String.valueOf(biomassComp[k]) + "*" + compoundNames[biomassNo[k]];
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

        compoundNames[noCompounds-1] = "Biomass";

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
    	/*
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
        }*/
    	
    	Writer writer = null;

    	try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFileName), "utf-8"));
		    for(int j=0; j<noReactions; j++){
		    	writer.write(reactions[j].name + ": " + answer[j] + "\n"); 
		    }
    	} catch (IOException ex) {
    		ex.printStackTrace();
    	} finally {
    		writer.close();	
    	}

    }

    public static void main(String[] args) {

        FBAreader fbaReader = new FBAreader(args[0],args[1]);

        fbaReader.createSmatrix();
        FBA fba = new FBA(fbaReader.noReactions,fbaReader.noCompounds);

        fba.loadS(fbaReader.S);

        fba.loadReactions(fbaReader.reactions);
        fba.loadNames(fbaReader.compoundNames);


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

