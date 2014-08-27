package model;

import java.io.*;
import java.io.File;
import jxl.*;
import jxl.write.*;
import java.lang.Math;
import java.util.*;
import java.awt.Color;

public class FBAreader {
    String inputFileName;
    String outputFileName;

    Workbook input;

    Sheet reactionList;
    Sheet compoundList;
    Sheet biomassList;

    int noReactions;
    int noCompounds;
    int noBiomass;

    String [] compoundNames;

    String [] reactionNames;
    String [] reactions;
    boolean [] lb;
    boolean [] ub;
    double [] a;
    double [] b;
    double [] c;

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

        compoundList = input.getSheet(0);
        reactionList = input.getSheet(1);
        biomassList = input.getSheet(2);

        noCompounds = compoundList.getRows();
        noReactions = reactionList.getRows();
        noBiomass = biomassList.getRows();

        compoundNames = new String [noCompounds+1];
        usedCompounds = new boolean [noCompounds+1];
        reactionNames = new String [noReactions+2];
        reactions = new String [noReactions+2];
        lb = new boolean [noReactions+2];
        ub = new boolean [noReactions+2];
        a = new double [noReactions+2];
        b = new double [noReactions+2];
        c = new double [noReactions+2];

        biomassNames = new String [noBiomass];
        biomassComp = new double [noBiomass];
        biomassIn = new boolean [noBiomass];

        for(int i = 0;i < noCompounds;i++) {
            usedCompounds[i] = false;
            compoundNames[i] = compoundList.getCell(0,i).getContents();
        }

        usedCompounds[noCompounds]=false;

        System.out.println("Compounds in");

        for(int j = 0;j < noReactions;j++) {
            reactionNames[j] = reactionList.getCell(0,j).getContents();
            reactions[j] = reactionList.getCell(1,j).getContents();

            Cell lq = reactionList.getCell(2,j);

            if(lq.getType() == CellType.LABEL) {
                String lc = lq.getContents();
                if(lc.compareTo("-infty") == 0) {
                    lb[j] = false;
                }
            }
            else {
                NumberCell nc = (NumberCell) lq;
                a[j] = nc.getValue();
                lb[j] = true;
            }


            Cell uq = reactionList.getCell(3,j);

            if(uq.getType() == CellType.LABEL) {
                String uc = uq.getContents();
                if(uc.compareTo("infty") == 0) {
                    ub[j] = false;
                }
            }
            else {
                NumberCell nc = (NumberCell) uq;
                b[j] = nc.getValue();
                ub[j] = true;
            }


            c[j] = Double.valueOf(reactionList.getCell(4,j).getContents());

            System.out.println(j + ": " + reactionNames[j] + " " + reactions[j] + " LB: " + a[j] + " UB: " + b[j]);

        }

        System.out.println("Reactions in");

        for(int k = 0;k < noBiomass;k++) {
            biomassNames[k] = biomassList.getCell(0,k).getContents();
            //System.out.println(biomassNames[k]);
            Cell nc1 = biomassList.getCell(1,k);
            //System.out.println();
            NumberCell nc = (NumberCell) nc1;
            biomassComp[k] = Double.valueOf(nc.getValue());
            //System.out.println(biomassComp[k]);
            biomassIn[k] = false;
            if(Double.valueOf(biomassList.getCell(2,k).getContents())>0) {
                biomassIn[k] = true;
            }

            System.out.println(k + ": " + biomassNames[k] + " " + biomassComp[k]);
        }

        biomassNo = new int [noBiomass];

        for(int i = 0;i < noCompounds;i++) {
            for(int k = 0;k < noBiomass;k++) {
                if(biomassNames[k].equals(compoundNames[i])) {
                    biomassNo[k] = i;
                    //System.out.println(biomassNames[k] + " is in, its compounts #" + i + " with stoich " + biomassComp[k]);
                }
            }
        }

        System.out.println("Biomass in");

        input.close();

    }


    public void createSmatrix() {
        S = new double [noCompounds+1][noReactions+2];
        //+2, +1 for the biomass reaction, growth and biomass to be included later.

        //this should work, but no error checking, woe betide!!!


        for(int j = 0;j < noReactions;j++) {
	    /*
	    System.out.print(j + " : " + reactions[j] + " " + a[j] + " " + b[j] + " " + c[j]);
	    if(lb[j]) {
		System.out.print(" Lower bound in place ");
	    }
	    else {
		System.out.print(" Lower bound absent ");
	    }
	    if(ub[j]) {
		System.out.print(" Upper bound in place ");
	    }
	    else {
		System.out.print(" Upper bound absent ");
	    }

	    System.out.print("\n");
	    */
            //pad to make sure bits splits it
            reactions[j] = reactions[j] + " ";

            //split the reactions in the two bits either side of the operator
            //note reversibility as we go
            String [] bits = new String [2];


            if(reactions[j].lastIndexOf("<==>") < 0) {
                bits = reactions[j].split("-->");
            }
            else {
                bits = reactions[j].split("<==>");

                //This is causing the problem, its overwriting what we need...
		/*
		if(a[j] == 0) {		    
		    lb[j] = false;
		}
		*/
            }


            //then break up the parts in the neg side
            String [] negParts;

            //System.out.println(bits[0].lastIndexOf("+"));

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

                //System.out.println(negParts[k]);

                boolean nothingFound = true;

                //now find them...
                for(int i = 0;i < noCompounds;i++) {
                    //System.out.print(negParts[k] + "," + compoundNames[i] + ":");
                    if(negParts[k].equals(compoundNames[i])) {
                        S[i][j] = -1*negStoich[k];
                        nothingFound = false;
                        usedCompounds[i] = true;
                        //System.out.print("FOUND " + S[i][j] + " inserted");
                    }
                    else {
                        //System.out.print("NUFFIN there ");
                        if(negParts[k].equals("")) {
                            nothingFound = false;
                        }
                    }
                    //System.out.print("\n");
                }
                if(nothingFound) {
                    System.out.println(negParts[k] + " was not found in reaction no " + j + ", " + reactionNames[j]);
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

                //System.out.println(posParts[k]);

                boolean nothingFound = true;

                //now find them...
                for(int i = 0;i < noCompounds;i++) {
                    if(posParts[k].equals(compoundNames[i])) {
                        S[i][j] = posStoich[k];
                        nothingFound = false;
                        usedCompounds[i] = true;
                    }
                    else {
                        //System.out.print("NUFFIN there ");
                        if(posParts[k].equals("")) {
                            nothingFound = false;
                        }
                    }

                }
                if(nothingFound) {
                    System.out.println(posParts[k] + " was not found in reaction no " + j + ", " + reactionNames[j]);
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

        reactionNames[noReactions-2] = "Biomass reaction";
        reactionNames[noReactions-1] = "Growth";

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
                //System.out.println(neg + "-->" + pos);
            }
            else {
                S[biomassNo[k]][noReactions-2] = 0.0;
            }
            //System.out.println(compoundNames[biomassNo[k]] + " is in with value " + biomassComp[k] + " should be the same as " + S[biomassNo[k]][noReactions-2]);
        }

        reactions[noReactions-2] = neg + "-->" + pos;
        reactions[noReactions-1] = "Biomass-->";


        noCompounds++; //one for the biomass

        compoundNames[noCompounds-1] = "Biomass";

        S[noCompounds-1][noReactions-2] = 1;
        S[noCompounds-1][noReactions-1] = -1;

        a[noReactions-2] = 0;
        b[noReactions-2] = 999;
        c[noReactions-2] = 0;
        lb[noReactions-2] = true;
        ub[noReactions-2] = false;

        a[noReactions-1] = 0;
        b[noReactions-1] = 999;
        c[noReactions-1] = 1;
        lb[noReactions-1] = true;
        ub[noReactions-1] = false;


        //Create new line in S matrix



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


            //write out a copy of the S Matrix to another sheet.
	    /*
	    for(int i = 0;i < noCompounds;i++) {	
		for(int j = 0;j < noReactions;j++) {
		    jxl.write.Number number = new jxl.write.Number(i+1,j+1,S[i][j]);
		    sheet.addCell(number);
		}
	    }
	    */

            WritableSheet sheet2 = copy.getSheet(1);
            WritableSheet sheet3 = copy.getSheet(2);
            jxl.write.Label label;
            jxl.write.Number number;

            label = new jxl.write.Label(0,noReactions-2, reactionNames[noReactions-2]);
            sheet2.addCell(label);
            label = new jxl.write.Label(1,noReactions-2, reactions[noReactions-2]);
            sheet2.addCell(label);
            number = new jxl.write.Number(2,noReactions-2, a[noReactions-2]);
            sheet2.addCell(number);
            number = new jxl.write.Number(3,noReactions-2, b[noReactions-2]);
            sheet2.addCell(number);
            number = new jxl.write.Number(4,noReactions-2, c[noReactions-2]);
            sheet2.addCell(number);
            label = new jxl.write.Label(0,noReactions-1, reactionNames[noReactions-1]);
            sheet2.addCell(label);
            label = new jxl.write.Label(1,noReactions-1, reactions[noReactions-1]);
            sheet2.addCell(label);
            number = new jxl.write.Number(2,noReactions-1, a[noReactions-1]);
            sheet2.addCell(number);
            number = new jxl.write.Number(3,noReactions-1, b[noReactions-1]);
            sheet2.addCell(number);
            number = new jxl.write.Number(4,noReactions-1, c[noReactions-1]);
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

        FBAreader R = new FBAreader(args[0],args[1]);
	/*
	for(int i = 0;i < R.noCompounds;i++) {
	    System.out.println(R.compoundNames[i]);
	}
	
	for(int j = 0;j < R.noReactions;j++) {
	    System.out.println(R.reactionNames[j] + " " + R.reactions[j]);
	}
	*/
        R.createSmatrix();
        System.out.print("Blank,");
        for(int i = 0;i < R.noCompounds;i++) {
            System.out.print(R.compoundNames[i] + ",");
        }
        System.out.print("\n");

        for(int j = 0;j < R.noReactions;j++) {
            System.out.print(R.reactionNames[j] + ",");
            for(int i = 0;i < R.noCompounds;i++) {
                System.out.print(R.S[i][j] + ",");
            }
            System.out.print("\n");
        }
	/*
	try {
	    R.writeSmatrix();
	}
	catch (IOException e) {
	    e.printStackTrace();
	}
	*/

        fba F = new fba(R.noReactions,R.noCompounds);

        F.loadS(R.S);

        System.out.println(R.a[933]);

        F.loadVectors(R.a,R.b,R.c,R.lb,R.ub);
        F.loadNames(R.compoundNames,R.reactionNames);


        double [] answer = new double[R.noReactions];
        answer = F.optimise();
        try {
            R.writeSmatrix(answer);
        }
        catch (IOException e) {
            e.printStackTrace();
        }


        PrintWriter pwout;
        FileWriter out;

        try {
            out = new FileWriter("sol.csv");
            pwout = new PrintWriter(out);

            for(int i = 0;i < R.noReactions;i++) {
                pwout.write(R.reactionNames[i] + "," + R.reactions[i]+ "," + answer[i] + "\n");
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

