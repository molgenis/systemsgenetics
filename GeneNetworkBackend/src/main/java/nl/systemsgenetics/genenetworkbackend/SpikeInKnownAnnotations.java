/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend;

import java.io.File;
import java.io.IOException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class SpikeInKnownAnnotations {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		
		final File hpoPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions.txt.gz");
		final File hpoAnnotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File hpoPredictionSpikedMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_spiked.txt");
		final int minZknownAnnotations = 3;
		
		DoubleMatrixDataset<String, String> hpoPredictionMatrix = DoubleMatrixDataset.loadDoubleData(hpoPredictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> hpoAnnotationMatrix = DoubleMatrixDataset.loadDoubleData(hpoAnnotationMatrixFile.getAbsolutePath());
		
		if (!hpoPredictionMatrix.getColObjects().equals(hpoAnnotationMatrix.getColObjects())) {
			System.err.println("Differnce in terms");
			return;
		}

		if (!hpoPredictionMatrix.getRowObjects().equals(hpoAnnotationMatrix.getRowObjects())) {
			System.err.println("Differnce in genes");
			return;
		}
		
		for(int i = 0; i < hpoAnnotationMatrix.rows() ; ++i){
			for(int j = 0; j < hpoAnnotationMatrix.columns(); ++j){
				if(hpoAnnotationMatrix.getElementQuick(i, j) > 0 && hpoPredictionMatrix.getElementQuick(i, j) < minZknownAnnotations){
					hpoPredictionMatrix.setElementQuick(i, j, minZknownAnnotations);
				}
			}
		}
		
		hpoPredictionMatrix.save(hpoPredictionSpikedMatrixFile);
		
	}
	
}
