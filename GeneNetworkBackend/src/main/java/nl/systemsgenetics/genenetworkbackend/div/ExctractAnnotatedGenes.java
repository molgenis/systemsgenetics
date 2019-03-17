/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ExctractAnnotatedGenes {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Extracts");

		HashSet<String> termsToExtract = new HashSet<>();
		termsToExtract.add("HP:0001644");
		termsToExtract.add("HP:0001638");

		DoubleMatrixDataset<String, String> annotationMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(annotationMatrixFile.getAbsolutePath(), '\t', null, termsToExtract);

		ArrayList<String> rowNames = annotationMatrix.getRowObjects();

		for (String term : termsToExtract) {

			System.out.println(term);
			
			DoubleMatrix1D termAnnoations = annotationMatrix.getCol(term);
			ArrayList<String> termGenes = new ArrayList<>();

			for (int r = 0; r < termAnnoations.size(); ++r) {
				if (termAnnoations.getQuick(r) > 0) {
					termGenes.add(rowNames.get(r));
				}
			}

			CSVWriter writer = new CSVWriter(new FileWriter(new File(outputFolder, term.replace(':', '_') + ".txt")), '\t', '\0', '\0', "\n");
			String[] outputLine = new String[1];
			int c = 0;
			outputLine[c++] = "Gene";
			writer.writeNext(outputLine);

			for (String gene : termGenes) {
				c = 0;
				outputLine[c++] = gene;
				writer.writeNext(outputLine);
			}
			
			writer.close();

		}

	}

}
