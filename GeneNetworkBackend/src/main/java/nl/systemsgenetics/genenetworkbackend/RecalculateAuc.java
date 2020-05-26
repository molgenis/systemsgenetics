/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class RecalculateAuc {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions.txt.gz");
		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File predictedHpoTermFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_auc_bonferroni.txt");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_auc_proteinCoding");
		final File ensgBiotypeFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgV83.txt");

		Map<String, String> ensgToBiotype = loadEnsgToBiotype(ensgBiotypeFile);

		final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();

		LinkedHashSet<String> predictedHpoTerms = readPredictedHpoTermFileBonCutoff(predictedHpoTermFile, 0.05);

		DoubleMatrixDataset<String, String> predictionMatrixFull = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrixFull = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());

		DoubleMatrixDataset<String, String> predictionMatrixPredicted = predictionMatrixFull.viewColSelection(predictedHpoTerms);
		DoubleMatrixDataset<String, String> annotationMatrixPredicted = annotationMatrixFull.viewColSelection(predictedHpoTerms);

		LinkedHashSet<String> proteinCodingGenes = new LinkedHashSet<>();
		for (String gene : predictionMatrixPredicted.getRowObjects()) {

			String geneBiotype = ensgToBiotype.get(gene);

			if (geneBiotype != null && geneBiotype.equals("protein_coding")) {

				proteinCodingGenes.add(gene);

			}

		}

		System.out.println("Protein coding in matrix: " + proteinCodingGenes.size());

		DoubleMatrixDataset<String, String> predictionMatrixPredictedProteinCoding = predictionMatrixPredicted.viewRowSelection(proteinCodingGenes);
		DoubleMatrixDataset<String, String> annotationMatrixPredictedProteinCoding = annotationMatrixPredicted.viewRowSelection(proteinCodingGenes);

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");
		String[] outputLine = new String[3];
		int c = 0;
		outputLine[c++] = "HPO";
		outputLine[c++] = "Origanl_AUC";
		outputLine[c++] = "ProteinCoding_AUC";
		writer.writeNext(outputLine);
		
		for (String term : predictionMatrixPredicted.getColObjects()) {

			DoubleMatrix1D hpoAnnotations = annotationMatrixPredicted.getCol(term);
			DoubleMatrix1D hpoPredictions = predictionMatrixPredicted.getCol(term);

			int hpoAnnotationCount = hpoAnnotations.cardinality();

			double[] zScoresAnnotatedGenes = new double[hpoAnnotationCount];
			double[] zScoresOtherGenes = new double[annotationMatrixPredicted.rows() - hpoAnnotationCount];

			int x = 0;
			int y = 0;

			for (int g = 0; g < hpoAnnotations.size(); g++) {
				if (hpoAnnotations.getQuick(g) != 0) {
					zScoresAnnotatedGenes[x++] = hpoPredictions.getQuick(g);
				} else {
					zScoresOtherGenes[y++] = hpoPredictions.getQuick(g);
				}
			}

			uTest.setData(zScoresAnnotatedGenes, zScoresOtherGenes);

			double originalAuc = uTest.getAuc();

			hpoAnnotations = annotationMatrixPredictedProteinCoding.getCol(term);
			hpoPredictions = predictionMatrixPredictedProteinCoding.getCol(term);

			hpoAnnotationCount = hpoAnnotations.cardinality();

			zScoresAnnotatedGenes = new double[hpoAnnotationCount];
			zScoresOtherGenes = new double[annotationMatrixPredicted.rows() - hpoAnnotationCount];

			x = 0;
			y = 0;

			for (int g = 0; g < hpoAnnotations.size(); g++) {
				if (hpoAnnotations.getQuick(g) != 0) {
					zScoresAnnotatedGenes[x++] = hpoPredictions.getQuick(g);
				} else {
					zScoresOtherGenes[y++] = hpoPredictions.getQuick(g);
				}
			}

			uTest.setData(zScoresAnnotatedGenes, zScoresOtherGenes);

			double proteinCodingAuc = uTest.getAuc();

			c = 0;
			outputLine[c++] = term;
			outputLine[c++] = String.valueOf(originalAuc);
			outputLine[c++] = String.valueOf(proteinCodingAuc);
			writer.writeNext(outputLine);

		}

	}

	private static Map<String, String> loadEnsgToBiotype(File ensgFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(ensgFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[0], nextLine[4]);

		}

		return Collections.unmodifiableMap(mapping);

	}

	private static LinkedHashSet<String> readPredictedHpoTermFileBonCutoff(File predictedHpoTermFile, double cutoff) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(predictedHpoTermFile))).withSkipLines(1).withCSVParser(parser).build();

		LinkedHashSet<String> hpos = new LinkedHashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if(Double.parseDouble(nextLine[4]) <= cutoff){
				hpos.add(nextLine[0]);
			}
			
		}

		reader.close();

		return hpos;

	}

}
