package model;

import java.lang.Math;
import java.util.*;
import java.io.*;
import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.GlpkException;
import org.gnu.glpk.SWIGTYPE_p_double;
import org.gnu.glpk.SWIGTYPE_p_int;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;
import org.gnu.glpk.glp_iocp;

public class fba {

    int m; // #compounds
    int n; // #reactions

    double [][] S; // Stoichiometry matrix
    double [] v; //Flux vector
    double [] a; // lower bounds
    double [] b; // upper bounds
    double [] c; // biomass vector;
    boolean [] lowerBound;
    boolean [] upperBound;
    String [] reactionNames;
    String [] reactantNames;
    String [] reactionSpec;
    int cValue; //biomass indicator

    public fba(int n, int m) {
        this.n = n;  //number of reactions
        this.m = m;  //number of compounds

        S = new double [m][n];
        v = new double [n];
        a = new double [n];
        lowerBound = new boolean [n];
        b = new double [n];
        upperBound = new boolean [n];
        c = new double [n];

        reactionNames = new String [n];
        reactantNames = new String [m];

    }

    public void loadS(double [][] S) {
        for(int i = 0;i < m;i++) {
            for(int j = 0;j < n;j++) {
                this.S[i][j] = S[i][j];
            }
        }

    }

    public void loadVectors(double [] a, double [] b, double [] c, boolean [] lowerBound, boolean [] upperBound) {
        for(int j = 0;j < n;j++) {
            this.a[j] = a[j];
            this.b[j] = b[j];
            this.c[j] = c[j];
            this.lowerBound[j] = lowerBound[j];
            this.upperBound[j] = upperBound[j];
            if(c[j] > 0) {
                cValue = j;
            }
        }
    }

    public void loadNames(String [] cn, String [] rn) {
        for(int j = 0;j < n;j++) {
            reactionNames[j] = rn[j];
        }
        for(int i = 0;i < m;i++) {
            reactantNames[i] = cn[i];
        }
    }

    public void dump() {
        for(int j = 0;j < n;j++) {
            System.out.print(reactionNames[j] + " : ");
            boolean firstone = true;
            for(int i = 0;i < m;i++) {

                if(S[i][j]<0) {
                    if(!firstone) {
                        System.out.print("+ ");
                    }
                    else {
                        firstone = false;
                    }

                    if(S[i][j]!=-1) {
                        System.out.print((-1*S[i][j]) + "*");
                    }
                    System.out.print(reactantNames[i]+" ");

                }
            }

            if(a[j] == 0) {
                System.out.print("--> ");
            }
            else {
                System.out.print("<==> ");
            }
            firstone = true;
            for(int i = 0;i < m;i++) {
                if(S[i][j]>0) {
                    if(!firstone) {
                        System.out.print("+ ");
                    }
                    else {
                        firstone = false;
                    }


                    if(S[i][j]!=1) {
                        System.out.print(S[i][j] + "*");
                    }
                    System.out.print(reactantNames[i] + " ");
                }
            }
            System.out.print("\n");
        }
    }

    public double [] optimise() {
        //create an augmented matrix for the simplex algorithm.
        glp_prob lp;
        glp_smcp parm;
        glp_iocp parm2;
        SWIGTYPE_p_int ind;
        SWIGTYPE_p_double val;
        int ret;
        double [] retArray = new double [n];

        // describe the optimization problem
        try {
            // Create problem
            lp = GLPK.glp_create_prob();
            System.out.println("Problem created");
            GLPK.glp_set_prob_name(lp, "myProblem");

            GLPK.glp_add_rows(lp, m);

            // Define columns
            GLPK.glp_add_cols(lp, n);//noOfReactions

            for(int k = 0;k < n;k++) {

                GLPK.glp_set_col_name(lp, k+1, reactionNames[k]);
                GLPK.glp_set_col_kind(lp, k+1, GLPKConstants.GLP_CV);

                if(lowerBound[k]) {
                    if(upperBound[k]) {
                        if(a[k] < b[k]) {
                            GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_DB, a[k], b[k]);
                        }
                        else {
                            GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_FX, a[k], b[k]);
                        }
                    }
                    else {
                        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_LO, a[k], b[k]);
                    }
                }
                else {
                    if(upperBound[k]) {
                        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_UP, a[k], b[k]);
                    }
                    else {
                        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_FR, a[k], b[k]);
                    }
                }
            }


            // Create constraints
            //first we count the number of entries;
            int ne = 0;
            for(int i = 0;i < m;i++) {
                for(int j = 0;j < n;j++) {
                    if(S[i][j]!=0) {
                        ne++;
                    }
                }
                GLPK.glp_set_row_name(lp, i+1, reactantNames[i]);
                GLPK.glp_set_row_bnds(lp, i+1, GLPKConstants.GLP_FX, 0.0, 0.0);
            }

            //System.out.println(ne);

            //then we create the parts of the triplet
            SWIGTYPE_p_int rowI = GLPK.new_intArray(ne+1);
            SWIGTYPE_p_int colI = GLPK.new_intArray(ne+1);
            SWIGTYPE_p_double values = GLPK.new_doubleArray(ne+1);

            //then we populate
            int counter = 1;
            for(int i = 0;i < m;i++) {
                for(int j = 0;j < n;j++) {
                    if(S[i][j]!=0) {
                        GLPK.intArray_setitem(rowI, counter, i+1);
                        GLPK.intArray_setitem(colI, counter, j+1);
                        GLPK.doubleArray_setitem(values, counter, S[i][j]);

                        counter++;
                    }
                }
            }
	    /*
	    for(int bob = 1;bob <= ne;bob++) {
		double a = GLPK.doubleArray_getitem(values,bob);
		double x = GLPK.intArray_getitem(rowI,bob);
		double y = GLPK.intArray_getitem(colI,bob);
		System.out.println(x + "  " + y + " is " + a+ " ");
	    }
	    */

            //add the matrix in one big glob
            GLPK.glp_load_matrix(lp,ne,rowI,colI,values);



            // Define objective
            GLPK.glp_set_obj_name(lp, "z");
            GLPK.glp_set_obj_dir(lp, GLPKConstants.GLP_MAX);
            GLPK.glp_set_obj_coef(lp, cValue+1, 1.);

            // Solve model
            parm = new glp_smcp();

            GLPK.glp_init_smcp(parm);
            parm.setPresolve(GLPKConstants.GLP_ON);
            parm.setMeth(GLPKConstants.GLP_DUALP);
            parm.setMsg_lev(GLPKConstants.GLP_MSG_OFF);
            //parm.setOut_frq(1);

            ret = GLPK.glp_simplex(lp, parm);

            parm2 = new glp_iocp();
            GLPK.glp_init_iocp(parm2);
            parm2.setMsg_lev(GLPKConstants.GLP_MSG_OFF);

            ret = GLPK.glp_intopt(lp, null);
            GLPK.glp_write_lp(lp, null, "problem.out");


            //GLPK.glp_write_lp(lp, null, "problem.out");

            // Retrieve solution
            if (ret == 0) {
                GLPK.glp_print_mip(lp, "data.out");
            } else {
                System.out.println("The problem could not be solved");
                GLPK.glp_print_mip(lp, "data.out");
            }



            for(int j = 0;j < n;j++) {

                retArray[j] = GLPK.glp_get_col_prim(lp, j+1);

            }
            // Free memory
            GLPK.glp_delete_prob(lp);
        } catch (GlpkException ex) {
            ex.printStackTrace();
        }

        return retArray;

    }

    /*
    public static void main(String[] args) {
	int noReactions = Integer.valueOf(args[0]);
	int noReagents = Integer.valueOf(args[1]);

	fba example = new fba(noReactions,noReagents,args[2]);
	example.input(args[2]);
	//example.dump();
	double [] answer = new double[noReactions];
	answer = example.optimise();
	
	PrintWriter pwout;
	FileWriter out;
	
	try {
	    out = new FileWriter(args[2] + "sol.csv");
	    pwout = new PrintWriter(out);

	    for(int i = 0;i < noReactions;i++) {
		pwout.write(example.reactionNames[i] + "," + example.reactionSpec[i]+ "," + answer[i] + "\n");
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
    */
}
