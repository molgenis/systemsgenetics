/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer;

import java.io.File;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedWriter;

/**
 *
 * This class is only used to preprocces the test gen-gen correlation matrix
 * 
 * @author patri
 */
public class PrepareRandomCorrelationsPerChunk {

	public static void main(String[] args) throws URISyntaxException, Exception {
		new PrepareRandomCorrelationsPerChunk().run();
	}

	public void run() throws URISyntaxException, Exception {
		
		File corFile = new File(this.getClass().getResource("/random/random_true_correlation_matrix.txt.gz").toURI());
		File genesFile = new File(this.getClass().getResource("/random/genes.txt").toURI());
		File outputFolder = new File(corFile.getParentFile(),"corPerArm");
		
		System.out.println(outputFolder.getAbsolutePath());
		
		final LinkedHashMap<String, List<Gene>> genes = IoUtils.readGenesAsChrArmMap(genesFile);
		
		DoubleMatrixDataset<String, String> fullCorMatrix = DoubleMatrixDataset.loadDoubleData(corFile.getAbsolutePath());
		
		for(Map.Entry<String, List<Gene>> chrArmGenes : genes.entrySet()){
			
			String arm = chrArmGenes.getKey();
			List<Gene> armGenes = chrArmGenes.getValue();
			
			ArrayList<String> genesWithCor = new ArrayList<>();
			
			for(Gene gene : armGenes){
				if(fullCorMatrix.containsCol(gene.getGene())){
					genesWithCor.add(gene.getGene());
				}
			}
			
			if(genesWithCor.isEmpty()){
				continue;
			}
			
			
			DoubleMatrixDataset<String, String> armCorMatrix = fullCorMatrix.viewColSelection(genesWithCor).viewRowSelection(genesWithCor);
			
			DoubleMatrixDatasetRowCompressedWriter.saveDataset(
							outputFolder.getAbsolutePath() + "/genecor_chr_" + arm + "_correlations",
							armCorMatrix,
							"Gene-Gene correlation matrix of chr " + arm, "Genes", "Genes");
					
			
			
			
		}

	}

}
