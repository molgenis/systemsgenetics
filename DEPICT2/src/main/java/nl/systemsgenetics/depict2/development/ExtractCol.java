/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.development;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ExtractCol {
	
	public static void extract(String inputFile, String column, String output) throws Exception{
		
		DoubleMatrixDataset<String, String> input = DoubleMatrixDataset.loadDoubleTextData(inputFile, '\t');
		
		DoubleMatrixDataset<String, String> outputMatrix = input.viewColSelection(column);
		
		outputMatrix.save(output);
		
	}
	
}
