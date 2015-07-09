/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import java.io.IOException;

/**
 *
 * @author lude
 */
public class EQTLInteractionAnalyser {


	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		// TODO code application logic here

        new TestEQTLDatasetForInteractions(args[0], args[1], args[2]);
    }
		
	}
