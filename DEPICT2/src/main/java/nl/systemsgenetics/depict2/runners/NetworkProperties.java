/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.runners;

import java.io.IOException;
import java.util.LinkedHashMap;

import nl.systemsgenetics.depict2.Depict2Options;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class NetworkProperties {
	
	public static void investigateNetwork(Depict2Options options) throws IOException {
		
		DoubleMatrixDataset<String, String> network = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		
		System.out.println("Loaded network");
		
		double[] thresholds = new double[]{1,2,3,5,10,20,30};
		
		LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<>();
		for(int i = 0 ; i < thresholds.length ; i++){
			hashCols.put(Double.toString(thresholds[i]), i);
		}
		
		DoubleMatrixDataset<String, String> connectivity = new DoubleMatrixDataset<>(network.getHashRows(), hashCols);
		
		for(int geneI = 0 ; geneI < network.rows() ; ++geneI){
			
			int[] geneThresholdCount = new int[thresholds.length];
			
			for(int geneJ = 0 ; geneJ < network.columns() ; ++geneJ){
				
				if(geneI == geneJ){
					continue;
				}
				
				double zScore = Math.abs(network.getElementQuick(geneI, geneJ));

				for(int thresholdI = 0; thresholdI < thresholds.length ; ++thresholdI){
					if(zScore >= thresholds[thresholdI]){
						geneThresholdCount[thresholdI]++;
					}
				}
				
			}
			
			for(int thresholdI = 0; thresholdI < thresholds.length ; ++thresholdI){
				connectivity.setElementQuick(geneI, thresholdI, geneThresholdCount[thresholdI]);
			}
			
		}
		
		connectivity.save(options.getOutputBasePath() + ".txt");
		
	}
	
}
