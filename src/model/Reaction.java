package model;

public class Reaction {
	public String name, equation;
	public boolean isLowerBoundActive, isUpperBoundActive;
    public Double lowerBound, upperBound, optimisationCoefficient;
    //unsure about optimisationCoefficient's name, but I think all are set to 0 except biomass so this fits
}
