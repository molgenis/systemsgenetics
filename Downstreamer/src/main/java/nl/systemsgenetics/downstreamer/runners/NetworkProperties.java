/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import java.io.IOException;
import java.util.*;

import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class NetworkProperties {

	private static final Logger LOGGER = Logger.getLogger(NetworkProperties.class);
	
	public static void investigateNetworkPatrick(DownstreamerOptions options) throws IOException, Exception {

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


	/**
	 * Take a co-regulation matrix and for each gene determine the number of connections that gene has. Reports at
	 * bonferoni significant z-scores, as well as the sumchisqr for each gene.
	 * A bit redundant with above, but didnt realise it was implemented already :'(
	 */
	public static void investigateNetwork(DownstreamerOptions options) throws Exception {
		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");

		DoubleMatrixDataset<String, String> corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		List<String> genesToKeep = corMatrix.getRowObjects();
		LOGGER.info("Read " + genesToKeep.size() + " genes in correlation matrix");
		genesToKeep.retainAll(genes.keySet());
		LOGGER.info("Retained " + genesToKeep.size() + " genes that overlap with --genes file");
		corMatrix = corMatrix.viewSelection(genesToKeep, genesToKeep);

		if (!corMatrix.getHashRows().keySet().containsAll(corMatrix.getHashCols().keySet())) {
			throw new Exception("Co-expression matrix is not squared with same row and col names");
		}

		if (!genes.keySet().containsAll(corMatrix.getHashRows().keySet())) {
			throw new Exception("Not all genes Co-expression matrix are found in gene mapping file");
		}

		final int genesInMatrix = corMatrix.rows();
		Set<String> colnames = new HashSet<>();
		colnames.add("sum_chi_sqr");
		colnames.add("z_nominal");
		colnames.add("z_bonf_sig");

		double zNominal = Math.abs(ZScores.pToZTwoTailed(0.05));
		double zBonfSig = Math.abs(ZScores.pToZTwoTailed(0.05/genesInMatrix));

		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>(corMatrix.getRowObjects(), colnames);

		for (int i = 0; i < genesInMatrix; ++i) {

			double curSumChiSqr = 0;
			double curZNominal = 0;
			double curZBonfSig = 0;

			for(int j =0 ; j < genesInMatrix; j++) {

				double curVal = corMatrix.getElementQuick(i, j);
				curSumChiSqr += (curVal * curVal);

				if (Math.abs(curVal) > zNominal) {
					curZNominal ++;
				}

				if (Math.abs(curVal) > zBonfSig) {
					curZBonfSig ++;
				}

			}

			output.setElementQuick(i, 0, curSumChiSqr);
			output.setElementQuick(i, 1, curZNominal);
			output.setElementQuick(i, 2, curZBonfSig);
		}

		LOGGER.info("Done, saving output");
		output.save(options.getOutputBasePath() + ".degree.tsv");

	}

}
