/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.functionenrichmentoftransqtls;

import java.io.File;
import java.io.IOException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CorrelateSumChi2ToPathways {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		final File pathwayMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File sumChi2MatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\");

		DoubleMatrixDataset pathwayMatrix = DoubleMatrixDataset.loadDoubleData(pathwayMatrixFile.getAbsolutePath());
		DoubleMatrixDataset sumChi2Matrix = DoubleMatrixDataset.loadDoubleData(sumChi2MatrixFile.getAbsolutePath());
		
		

	}

}
