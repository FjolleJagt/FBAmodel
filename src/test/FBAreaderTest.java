package test;

import static org.junit.Assert.*;
import model.Compound;
import model.FBAreader;

import org.junit.Before;
import org.junit.Test;


public class FBAreaderTest {
	FBAreader fbaReader;
	
	@Before
	public void createFBAReader(){
		fbaReader = new FBAreader("resources/in/original.xls","resources/out.xls");
	}

	@Test
	public void testGetStoichiometryFromEquation1() {
		double[] stoichiometry = fbaReader.getStoichiometryFromEquation("C00079_[cyt]+C00026_[cyt]<==>C00166_[cyt]+C00025_[cyt]");
		Compound cyt79 = fbaReader.compounds[90];
		Compound cyt26 = fbaReader.compounds[44];
		Compound cyt166 = fbaReader.compounds[148];
		Compound cyt25 = fbaReader.compounds[42];
		assertTrue("(bad test) Row 91 is " + cyt79.name + " not C00079_[cyt]",cyt79.name.equals("C00079_[cyt]"));
		assertTrue("(bad test) Row 45 is " + cyt26.name + " not C00026_[cyt]",cyt26.name.equals("C00026_[cyt]"));
		assertTrue("(bad test) Row 149 is " + cyt166.name + " not C00166_[cyt]",cyt166.name.equals("C00166_[cyt]"));
		assertTrue("(bad test) Row 43 is " + cyt25.name + " not C00025_[cyt]",cyt25.name.equals("C00025_[cyt]"));
		assertTrue(cyt79.name + " doesn't have coefficient -1 (actually " + stoichiometry[90] + ")",stoichiometry[90] == -1);
		assertTrue(cyt26.name + " doesn't have coefficient -1 (actually " + stoichiometry[40] + ")",stoichiometry[44] == -1);
		assertTrue(cyt166.name + " doesn't have coefficient 1 (actually " + stoichiometry[148] + ")",stoichiometry[148] == 1);
		assertTrue(cyt25.name + " doesn't have coefficient 1 (actually " + stoichiometry[42] + ")",stoichiometry[42] == 1);
	}

	@Test
	public void testGetStoichiometryFromEquation101(){
		double[] stoichiometry = fbaReader.getStoichiometryFromEquation("C02094_[cyt]+C00005_[cyt]+C00080_[cyt]+C00007_[cyt]-->C08592_[cyt]+C00006_[cyt]+C00001_[cyt]+C00080_[cyt]");
		Compound cyt2094 = fbaReader.compounds[334];
		Compound cyt80 = fbaReader.compounds[92];
		assertTrue("(bad test) Row 335 is " + cyt2094.name + " not C02094_[cyt]",cyt2094.name.equals("C02094_[cyt]"));
		assertTrue("(bad test) Row 93 is " + cyt80.name + " not C00080_[cyt]",cyt80.name.equals("C00080_[cyt]"));
		assertTrue(cyt2094.name + " doesn't have coefficient -1 (actually " + stoichiometry[334] + ")",stoichiometry[334] == -1);
		assertTrue(cyt80.name + " doesn't have coefficient 0 (actually " + stoichiometry[92] + ")",stoichiometry[92] == 0);
	}
}
