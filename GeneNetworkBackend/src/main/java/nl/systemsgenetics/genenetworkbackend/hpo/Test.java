/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import java.io.File;
import java.io.IOException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class Test {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		final File updatedPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_improved.txt.gz");
		
		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleData(updatedPredictionMatrixFile.getAbsolutePath());
		
		System.out.println(matrix.getElement("ENSG00000165917", "HP:0001324"));
		
	}
	
}
