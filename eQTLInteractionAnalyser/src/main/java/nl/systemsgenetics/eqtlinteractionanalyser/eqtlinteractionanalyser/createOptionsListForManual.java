/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 *
 * @author patri
 */
public class createOptionsListForManual {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		
		Options options = EQTLInteractionAnalyser.OPTIONS;
		
		System.out.println("| Short | Long | Description |");
		System.out.println("|-------|------|-------------|");
		
		for(Object optionO : options.getOptions()){
			Option option = (Option) optionO;
			System.out.print("| -");
			System.out.print(option.getOpt());
			System.out.print(" | --");
			System.out.print(option.getLongOpt());
			System.out.print(" | ");
			System.out.print(option.getDescription());
			System.out.println(" | ");
		}
		
	}
	
}
