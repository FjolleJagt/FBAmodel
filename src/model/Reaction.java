package model;

public class Reaction {
	String name, equation;
	boolean isLowerBoundActive, isUpperBoundActive;
    Double lowerBound, upperBound, optimisationCoefficient;
    //unsure about optimisationCoefficient's name, but I think all are set to 0 except biomass so this fits
}
