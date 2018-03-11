/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class HpoGenePrioritisation {

	public static void main(String[] args) throws IOException, ParseException {

		final File hpoPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions.txt");
		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations");
		final File ensgSymbolMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\selectedHpo2.txt");

		Map<String, String> ensgSymbolMapping = loadEnsgToHgnc(ensgSymbolMappingFile);

		HashMap<String, LinkedHashSet<String>> caseHpo = loadCaseHpo(caseHpoFile);

		DoubleMatrixDataset<String, String> hpoPredictionMatrix = DoubleMatrixDataset.loadDoubleData(hpoPredictionMatrixFile.getAbsolutePath());
		ArrayList<String> genes = hpoPredictionMatrix.getRowObjects();


		System.out.println("Done loading data");

		BufferedWriter sampleFileWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "samples.txt")));

		for (Map.Entry<String, LinkedHashSet<String>> caseHpoEntry : caseHpo.entrySet()) {

			String caseId = caseHpoEntry.getKey();
			LinkedHashSet<String> hpo = caseHpoEntry.getValue();
			double zSum = 0;

			sampleFileWriter.append(caseId);
			sampleFileWriter.append('\n');

			System.out.println("Processing: " + caseId);

			for (String term : hpo) {
				if (!hpoPredictionMatrix.containsCol(term)) {
					System.err.println("Missing HPO: " + term);
				}
			}

			DoubleMatrixDataset<String, String> predictionCaseTerms = hpoPredictionMatrix.viewColSelection(hpo);
			DoubleMatrix2D predictionCaseTermsMatrix = predictionCaseTerms.getMatrix();
			GenePrioritisationResult[] geneResults = new GenePrioritisationResult[genes.size()];

			double denominator = Math.sqrt(hpo.size());

			for (int g = 0; g < predictionCaseTermsMatrix.rows(); ++g) {
				String gene = genes.get(g);
				String symbol = ensgSymbolMapping.get(gene);
				if (symbol == null) {
					symbol = "";
				}
				double geneScore = predictionCaseTermsMatrix.viewRow(g).zSum() / denominator;

				geneResults[g] = new GenePrioritisationResult(gene, symbol, geneScore);

			}

			Arrays.sort(geneResults);

			CSVWriter writer = new CSVWriter(new FileWriter(new File(outputFolder, caseId + ".txt")), '\t', '\0', '\0', "\n");

			String[] outputLine = new String[4 + hpo.size()];
			int c = 0;
			outputLine[c++] = "Ensg";
			outputLine[c++] = "Hgnc";
			outputLine[c++] = "Rank";
			outputLine[c++] = "Zscore";
			for (String term : hpo) {
				outputLine[c++] = term;
			}
			writer.writeNext(outputLine);

			int rank = 1;
			for (GenePrioritisationResult geneResult : geneResults) {
				c = 0;
				outputLine[c++] = geneResult.getEnsg();
				outputLine[c++] = geneResult.getSymbol();
				outputLine[c++] = String.valueOf(rank++);
				outputLine[c++] = String.valueOf(geneResult.getGeneScore());
				DoubleMatrix1D geneZs = predictionCaseTerms.viewRow(geneResult.getEnsg());
				for (int i = 0; i < hpo.size(); i++) {
					outputLine[c++] = String.valueOf(geneZs.getQuick(i));
				}
				writer.writeNext(outputLine);
			}

			writer.close();

		}

		sampleFileWriter.close();

	}

	private static HashMap<String, LinkedHashSet<String>> loadCaseHpo(File caseHpoFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(caseHpoFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, LinkedHashSet<String>> caseHpo = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[5].isEmpty()) {

				LinkedHashSet<String> hpo = caseHpo.get(nextLine[0]);
				if (hpo == null) {
					hpo = new LinkedHashSet<>();
					caseHpo.put(nextLine[0], hpo);
				}
				hpo.add(nextLine[1]);

			}

		}

		return caseHpo;

	}

	private static Map<String, String> loadEnsgToHgnc(File mappingFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mappingFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[0], nextLine[1]);

		}

		return Collections.unmodifiableMap(mapping);

	}

}
