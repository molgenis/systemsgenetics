/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 *
 * @author patri
 */
public class NetworkProperties {
	
	public static void investigateNetwork(DownstreamerOptions options) throws IOException, Exception {
		
		
		List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
		
		
		DoubleMatrixDatasetFastSubsetLoader networkLoader = new DoubleMatrixDatasetFastSubsetLoader(options.getGwasZscoreMatrixPath());

		Set<String> genesInNetwork = networkLoader.getOriginalRowMap().keySet();
		
        ArrayList<String> geneNames =  new ArrayList<>(genes.size());
		for(Gene gene : genes){
			if(genesInNetwork.contains(gene.getGene())){
				geneNames.add(gene.getGene());
			}
		}
		
		DoubleMatrixDataset<String, String> network = networkLoader.loadSubsetOfRowsBinaryDoubleData(geneNames);
		network = network.viewColSelection(geneNames);
		
		System.out.println("Loaded network");
		
		double[] thresholds = new double[]{1,2,3,4,5,6,8,10,20,30};
		
		LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<>();
		for(int i = 0 ; i < thresholds.length ; i++){
			hashCols.put(Double.toString(thresholds[i]), i);
		}
		hashCols.put("sumChi2", thresholds.length);
		
		
		DoubleMatrixDataset<String, String> connectivity = new DoubleMatrixDataset<>(network.getHashRows(), hashCols);
		
		for(int geneI = 0 ; geneI < network.rows() ; ++geneI){
			
			int[] geneThresholdCount = new int[thresholds.length];
			double sumChi2 = 0;
			
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
				
				sumChi2 += zScore * zScore;
				
			}
			
			for(int thresholdI = 0; thresholdI < thresholds.length ; ++thresholdI){
				connectivity.setElementQuick(geneI, thresholdI, geneThresholdCount[thresholdI]);
			}
			
			connectivity.setElementQuick(geneI, thresholds.length, sumChi2);
			
		}
		
		connectivity.save(options.getOutputBasePath() + ".txt");
		
	}
	
}
