package model;

import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.GlpkException;
import org.gnu.glpk.SWIGTYPE_p_double;
import org.gnu.glpk.SWIGTYPE_p_int;
import org.gnu.glpk.glp_iocp;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;

public class FBA {

    int noCompounds; // #compounds
    int noReactions; // #reactions

    double [][] S; // Stoichiometry matrix
    double [] v; //Flux vector
    Reaction [] reactions;
    Compound [] compounds;
    String [] reactionSpec;
    int biomassEquationIndex; //biomass indicator

    public FBA(int noReactions, int noCompounds) {
        this.noReactions = noReactions;
        this.noCompounds = noCompounds;

        S = new double [noCompounds][noReactions];
        v = new double [noReactions];

        compounds = new Compound [noCompounds];

    }

    public void loadS(double [][] S) {
        for(int i = 0;i < noCompounds;i++) {
            for(int j = 0;j < noReactions;j++) {
                this.S[i][j] = S[i][j];
            }
        }
    }

    public void loadReactions(Reaction[] reactions) {
    	this.reactions = reactions;
        for(int j = 0;j < noReactions;j++) {
            if(reactions[j].optimisationCoefficient > 0) {
                biomassEquationIndex = j;
            }
        }
    }

    public void loadNames(Compound [] compounds) {
        for(int i = 0;i < noCompounds;i++) {
            this.compounds[i] = compounds[i];
        }
    }

    public void dump() {
        for(int j = 0;j < noReactions;j++) {
            System.out.print(reactions[j].name + " : ");
            boolean firstone = true;
            for(int i = 0;i < noCompounds;i++) {

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
                    System.out.print(compounds[i]+" ");

                }
            }

            if(reactions[j].lowerBound == 0) {
                System.out.print("--> ");
            }
            else {
                System.out.print("<==> ");
            }
            firstone = true;
            for(int i = 0;i < noCompounds;i++) {
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
                    System.out.print(compounds[i] + " ");
                }
            }
            System.out.print("\n");
        }
    }

    public double [] optimise() throws RuntimeException {

        double [] optimalFluxes = new double [noReactions];

        glp_prob problem = GLPK.glp_create_prob();
        GLPK.glp_set_prob_name(problem, "FBA");

        setReactionBounds(problem);
        setCompoundsToBeConserved(problem);
        setSMatrix(problem);
        setObjective(problem);

        presolveProblemWithSimplex(problem);
        solveWithIntegerOptimisation(problem);
        
        GLPK.glp_write_lp(problem, null, "problem.out");
        GLPK.glp_print_mip(problem, "data.out");

        for(int j=0; j < noReactions; j++) {
            optimalFluxes[j] = GLPK.glp_get_col_prim(problem, j+1);
        }
        GLPK.glp_delete_prob(problem);

        return optimalFluxes;

    }

	private void solveWithIntegerOptimisation(glp_prob problem) {
		glp_iocp integerOptimiserParameters = new glp_iocp();
        GLPK.glp_init_iocp(integerOptimiserParameters);
        integerOptimiserParameters.setMsg_lev(GLPKConstants.GLP_MSG_OFF);
        int returnCode = GLPK.glp_intopt(problem, null);
        if(returnCode == GLPK.GLP_ENODFS){
        	throw new GlpkException("Unable to start the search, your bounds aren't restrictive enough!");
        } else if(returnCode == GLPK.GLP_ENOPFS){
        	throw new GlpkException("Unable to start the search, your bounds are too restrictive!");
        } else if(returnCode != 0){
        	throw new GlpkException("Unable to solve the problem - look up error code and put an explanatory message for this code here.");
        }
	}

	private void presolveProblemWithSimplex(glp_prob problem) {
		glp_smcp simplexParameters = new glp_smcp();

        GLPK.glp_init_smcp(simplexParameters);
        simplexParameters.setPresolve(GLPKConstants.GLP_ON);
        simplexParameters.setMeth(GLPKConstants.GLP_DUALP);
        simplexParameters.setMsg_lev(GLPKConstants.GLP_MSG_OFF);

        int returnCode = GLPK.glp_simplex(problem, simplexParameters);
        
        if(returnCode == GLPK.GLP_ENODFS){
    		throw new GlpkException("Unable to start the search, because the bounds are relaxed enough to allow unlimited growth");
        } else if(returnCode == GLPK.GLP_ENOPFS){
        	throw new GlpkException("Unable to start the search, because the bounds are too restrictive");
        } else if(returnCode != 0){
			throw new GlpkException("Unable to solve the problem - look up error code and put an explanatory message for this code here.");
        }
	}

	private void setObjective(glp_prob problem) {
		GLPK.glp_set_obj_name(problem, "z");
        GLPK.glp_set_obj_dir(problem, GLPKConstants.GLP_MAX);
        GLPK.glp_set_obj_coef(problem, biomassEquationIndex+1, 1.);
	}

	private void setSMatrix(glp_prob lp) {
		int entryCount = 0;
		for(int i = 0; i < noCompounds; i++){
			for(int j = 0; j < noReactions; j++){
				if(S[i][j] != 0){
					entryCount++;
				}
			}
		}
		
		SWIGTYPE_p_int rowIndices = GLPK.new_intArray(entryCount+1);
		SWIGTYPE_p_int colIndices = GLPK.new_intArray(entryCount+1);
		SWIGTYPE_p_double values = GLPK.new_doubleArray(entryCount+1);

		int counter = 1;
		for(int i = 0;i < noCompounds;i++) {
		    for(int j = 0;j < noReactions;j++) {
		        if(S[i][j]!=0) {
		            GLPK.intArray_setitem(rowIndices, counter, i+1);
		            GLPK.intArray_setitem(colIndices, counter, j+1);
		            GLPK.doubleArray_setitem(values, counter, S[i][j]);
		            counter++;
		        }
		    }
		}

		GLPK.glp_load_matrix(lp,entryCount,rowIndices,colIndices,values);
	}

	private void setCompoundsToBeConserved(glp_prob lp) {
		GLPK.glp_add_rows(lp, noCompounds);
		for(int i = 0;i < noCompounds;i++) {
		    GLPK.glp_set_row_name(lp, i+1, compounds[i].name);
		    GLPK.glp_set_row_bnds(lp, i+1, GLPKConstants.GLP_FX, 0.0, 0.0);
		}
	}

	private void setReactionBounds(glp_prob lp) throws RuntimeException {
		GLPK.glp_add_cols(lp, noReactions);

		for(int k=0; k < noReactions; k++) {
		    GLPK.glp_set_col_name(lp, k+1, reactions[k].name);
		    GLPK.glp_set_col_kind(lp, k+1, GLPKConstants.GLP_CV); //continuous variable

		    if(reactions[k].lowerBound != null && reactions[k].upperBound != null) {
		        if(reactions[k].lowerBound < reactions[k].upperBound) {
		            GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_DB, reactions[k].lowerBound, reactions[k].upperBound);
		        } else if(Math.abs(reactions[k].lowerBound - reactions[k].upperBound) < 0.0000001 ){
		        	// bounds are close enough to be treated the same: floating point arithmetic issues
		            GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_FX, reactions[k].lowerBound, reactions[k].lowerBound);
		        } else {
		        	throw new RuntimeException("Reaction " + reactions[k].name + " has a lower bound (" + reactions[k].lowerBound + ") greater than its upper bound (" + reactions[k].upperBound + ")");
		        }
		    } else if(reactions[k].lowerBound != null && reactions[k].upperBound == null){
		        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_LO, reactions[k].lowerBound, 999.0);
		    } else if(reactions[k].lowerBound == null && reactions[k].upperBound != null){
		        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_UP, -999.0, reactions[k].upperBound);
		    } else {
		        GLPK.glp_set_col_bnds(lp, k+1, GLPKConstants.GLP_FR,-999.0, 999.0);
		    }
		}
	}
}
