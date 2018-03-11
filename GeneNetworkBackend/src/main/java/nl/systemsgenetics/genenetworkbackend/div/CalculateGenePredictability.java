/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class CalculateGenePredictability {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions.txt.gz");
		File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Ensembl2Reactome_All_Levels.txt_matrix.txt.gz");
		File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions_bonSigTerms.txt");
		File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\reactome_predictions_genePredictability.txt");

		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		DoubleMatrixDataset<String, String> predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrix = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());

		DoubleMatrixDataset<String, String> predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);
		DoubleMatrixDataset<String, String> annotationMatrixSignificant = annotationMatrix.viewColSelection(significantTerms);

		if (!predictionMatrixSignificant.getColObjects().equals(annotationMatrixSignificant.getColObjects())) {
			System.err.println("Differnce in terms");
			return;
		}

		if (!predictionMatrixSignificant.getRowObjects().equals(annotationMatrixSignificant.getRowObjects())) {
			System.err.println("Differnce in genes");
			return;
		}

		MannWhitneyUTest2 uTest = new MannWhitneyUTest2();

		double[] genePredictabilityZscores = new double[predictionMatrixSignificant.rows()];

		for (int g = 0; g < predictionMatrixSignificant.rows(); g++) {

			DoubleMatrix1D geneAnnotations = annotationMatrixSignificant.getRow(g);

			int geneAnnotationCount = geneAnnotations.cardinality();

			if (geneAnnotationCount >= 10) {

				double[] zScoresAnnotatedPathways = new double[geneAnnotationCount];
				double[] zScoresOtherPathways = new double[annotationMatrixSignificant.columns() - geneAnnotationCount];

				int x = 0;
				int y = 0;

				for (int p = 0; p < geneAnnotations.size(); p++) {
					if (geneAnnotations.getQuick(p) != 0) {
						zScoresAnnotatedPathways[x++] = predictionMatrixSignificant.getElementQuick(g, p);
					} else {
						zScoresOtherPathways[y++] = predictionMatrixSignificant.getElementQuick(g, p);
					}
				}

				uTest.setData(zScoresOtherPathways, zScoresAnnotatedPathways);
				genePredictabilityZscores[g] = uTest.getZ();

			} else {
				genePredictabilityZscores[g] = Double.NaN;
			}

		}

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "Gene";
		outputLine[c++] = "Z-score";
		writer.writeNext(outputLine);

		ArrayList<String> geneNames = predictionMatrixSignificant.getRowObjects();
		for (int g = 0; g < predictionMatrixSignificant.rows(); g++) {
			c = 0;
			outputLine[c++] = geneNames.get(g);
			outputLine[c++] = String.valueOf(genePredictabilityZscores[g]);
			writer.writeNext(outputLine);
		}
		
		writer.close();

	}

	private static LinkedHashSet<String> loadSignificantTerms(File significantTermsFile) throws IOException {

		LinkedHashSet<String> significantTerms = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(significantTermsFile));

		String line;
		while ((line = reader.readLine()) != null) {
			significantTerms.add(line);
		}

		return significantTerms;

	}

}
