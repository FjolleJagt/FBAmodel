package model;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.lang.Math;
import java.lang.Number;
import java.awt.image.*;
import java.io.*;
import java.util.ArrayList;

// import java.lang.Object;

public class mas implements ActionListener {

	// All the Swing objects are initialised here

	JFrame controlFrame, interactionFrame;
	JPanel controlPanel, interactionPanel;

	JLabel textinLabel, textoutLabel, div1Label,div2Label, growthLabel, answerLabel, lookLabel;
	JLabel imLabel,barLabel;
	JButton startSim, okayIT;
	JTextField textinField,textoutField;
	JSpinner lookSpin;

	JMenuBar menuBar;
	JMenu Moptions;
	JCheckBoxMenuItem cbOut;
	JMenuItem cbIV;

	// Internal variables

	boolean Outputstate;
	int wdg,hdg;

	// Dependent objects

	FBAreader fbaReader; // MAIN dependent object

	public mas() {

		//Create and set up the window.
		controlFrame = new JFrame("Control Panel");
		controlFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		controlFrame.setSize(new Dimension(800, 600));

		//Create and set up the panel.
		controlPanel = new JPanel(new GridLayout(7, 2));

		controlFrame.setJMenuBar(createMenuBar());

		//Add the widgets.
		addWidgets1();
		//addWidgets2();

		//Set the default button.
		controlFrame.getRootPane().setDefaultButton(startSim);

		//Add the panel to the window.
		controlFrame.getContentPane().add(controlPanel, BorderLayout.WEST);

		//Display the window.
		controlFrame.pack();
		controlFrame.setVisible(true);

		Outputstate =false;

	}

	/*
	 * Create and add the widgets.
	 */
	private void addWidgets1() {
		//Create widgets.

		//Make the Insolation spinner and add it
		textinField = new JTextField("iHK677c.xls");
		textinLabel = new JLabel("input filename", SwingConstants.LEFT);	
		controlPanel.add(textinLabel);
		controlPanel.add(textinField);

		textoutField = new JTextField("output.xls");
		textoutLabel = new JLabel("output filename", SwingConstants.LEFT);	
		controlPanel.add(textoutLabel);
		controlPanel.add(textoutField);

		SpinnerModel lookModel = new SpinnerNumberModel(1, 1, 10000, 1);
		lookLabel = new JLabel("Look at", SwingConstants.LEFT);
		lookSpin = new JSpinner(lookModel);		
		controlPanel.add(lookLabel);
		controlPanel.add(lookSpin);
		lookLabel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

		growthLabel = new JLabel("Growth", SwingConstants.LEFT);
		answerLabel = new JLabel("");
		controlPanel.add(growthLabel);
		controlPanel.add(answerLabel);

		textinLabel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		textoutLabel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

		div1Label = new JLabel("--------", SwingConstants.CENTER);
		controlPanel.add(div1Label);
		div2Label = new JLabel("--------", SwingConstants.CENTER);
		controlPanel.add(div2Label);

		//Make the buttons
		startSim = new JButton("Compute");
		startSim.setActionCommand("Compute");

		//Listen to events from the buttons.
		startSim.addActionListener(this);

		//Add these widgets to the container
		controlPanel.add(startSim);

	}

	public void actionPerformed(ActionEvent event) {

		if("wakeup".equals(event.getActionCommand())){
			interactionFrame.setVisible(true);
		}
		else if("ok".equals(event.getActionCommand())) {
			interactionFrame.setVisible(false);
		}
		else {

			if ("Compute".equals(event.getActionCommand())) {
				fbaReader = new FBAreader("resources/in.xls","resources/out.xls");
				
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
			
				//To here
			
				System.out.println("Growth is " + answer[fbaReader.noReactions-1]);
			
				for(int k = 0;k < inBiomass.length;k++) {
					fbaReader.biomassIn[k] = inBiomass[k];		   
				}
			
				try {		    
					fbaReader.writeSmatrix(answer);
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
	}


	public JMenuBar createMenuBar() {

		//Create the menu bar.
		menuBar = new JMenuBar();

		//Build the first menu.
		Moptions = new JMenu("Options");
		Moptions.setMnemonic(KeyEvent.VK_A);
		Moptions.getAccessibleContext().setAccessibleDescription(
				"The only menu in this program that has menu items");
		menuBar.add(Moptions);

		//a group of check box menu items

		cbOut = new JCheckBoxMenuItem("Output Monitor",true);
		cbOut.setMnemonic(KeyEvent.VK_O);
		Moptions.add(cbOut);
		Moptions.addSeparator();
		//a popup box

		cbIV = new JMenuItem("Change Interaction");
		cbIV.setMnemonic(KeyEvent.VK_I);
		cbIV.setActionCommand("wakeup");
		cbIV.addActionListener(this);
		Moptions.add(cbIV);

		return menuBar;
	}

	/**
	 * Create the GUI and show it.  For thread safety,
	 * this method should be invoked from the
	 * event-dispatching thread.
	 */

	private static void createAndShowGUI() {
		//Make sure we have nice window decorations.
		JFrame.setDefaultLookAndFeelDecorated(true);

		mas control = new mas();

	}

	public static void main(String[] args) {
		//Schedule a job for the event-dispatching thread:
		//creating and showing this application's GUI.

		/*javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				createAndShowGUI();
			}
		});*/
		TestRig testRig = new TestRig();
		testRig.test();
	}
}
