/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 *
 * @author patri
 */
public class First1000qtl {

	/**
	 * @param options
	 * @throws java.io.IOException
	 */
	public static void printFirst1000(DownstreamerOptionsDeprecated options) throws IOException, Exception {
		
		System.out.println("");
		System.out.println("WARNING: This is an undocumented function used during development");
		System.out.println("");
		
		DoubleMatrixDatasetFastSubsetLoader loader = new DoubleMatrixDatasetFastSubsetLoader(options.getGwasZscoreMatrixPath());
		
		ArrayList<String> rowsToOutput = new ArrayList<>();
		
		Set<String> rowMap = loader.getOriginalRowMap(); 
		
		int i = 0;
		for(String row : rowMap){
			if(i++ > 1000){
				break;
			}
			rowsToOutput.add(row);
		}
		
		DoubleMatrixDataset<String, String> subset = loader.loadSubsetOfRowsBinaryDoubleData(rowsToOutput);
		
		subset.save(options.getOutputBasePath());
		
		
	}
	
}
