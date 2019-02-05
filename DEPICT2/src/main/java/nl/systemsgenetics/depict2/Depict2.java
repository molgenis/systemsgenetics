/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import org.apache.commons.cli.ParseException;

/**
 *
 * @author patri
 */
public class Depict2 {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		
		Depict2Options options;

		if (args.length == 0) {
			Depict2Options.printHelp();
			return;
		}

		try {
			options = new Depict2Options(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			Depict2Options.printHelp();
			return;
		}

		options.printOptions();
		
		
	}
	
}
