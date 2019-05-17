/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import org.apache.log4j.Logger;
import org.biojava.nbio.ontology.Term;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class HpoGenePrioritisation {

	private static final NumberFormat Z_FORMAT = new DecimalFormat("#0.0##");
	private static final Logger LOGGER = Logger.getLogger(HpoGenePrioritisation.class);

	public static void prioritize(GadoOptions options) throws IOException, ParseException, Exception {

		final File hpoPredictionMatrixFile = options.getPredictionMatrixFile();
		final File ensgSymbolMappingFile = options.getGenesFile();
		//final File hpoOboFile = options.getHpoOboFile();

		final File caseHpoFile = options.getProcessedCaseHpoFile();
		final File outputFolder = new File(options.getOutputBasePath());

		outputFolder.mkdirs();

		Map<String, String> ensgSymbolMapping = loadEnsgToHgnc(ensgSymbolMappingFile);

		HashMap<String, LinkedHashSet<String>> caseHpo = loadCaseHpo(caseHpoFile);

		//HpoOntology hpoOntology = new HpoOntology(hpoOboFile);

		DoubleMatrixDataset<String, String> hpoPredictionMatrix = DoubleMatrixDataset.loadDoubleBinaryData(hpoPredictionMatrixFile.getAbsolutePath());
		ArrayList<String> genes = hpoPredictionMatrix.getRowObjects();

		LOGGER.info("Loaded HPO prediction matrix");

		BufferedWriter sampleFileWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "samples.txt")));

		for (Map.Entry<String, LinkedHashSet<String>> caseHpoEntry : caseHpo.entrySet()) {

			String caseId = caseHpoEntry.getKey();
			LinkedHashSet<String> hpo = caseHpoEntry.getValue();

			sampleFileWriter.append(caseId);
			sampleFileWriter.append('\n');

			for (String term : hpo) {
				if (!hpoPredictionMatrix.containsCol(term)) {
					throw new RuntimeException("Missing HPO: " + term);
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
			//outputLine[c++] = "PossibleAdditionalHpos";
			writer.writeNext(outputLine);

			int rank = 1;
			for (GenePrioritisationResult geneResult : geneResults) {
				c = 0;
				outputLine[c++] = geneResult.getEnsg();
				outputLine[c++] = geneResult.getSymbol();
				outputLine[c++] = String.valueOf(rank++);
				outputLine[c++] = Z_FORMAT.format(geneResult.getGeneScore());
				DoubleMatrix1D geneZs = predictionCaseTerms.viewRow(geneResult.getEnsg());
				for (int i = 0; i < hpo.size(); i++) {
					outputLine[c++] = Z_FORMAT.format(geneZs.getQuick(i));
				}

				//outputLine[c++] = getOtherPossibleHpoTerms(geneResult.getEnsg(), hpo, hpoPredictionMatrixSignificant, minZscoreOtherCandidates, hpoOntology);

				writer.writeNext(outputLine);
			}

			writer.close();
			
			LOGGER.info("Finished gene prioritization for " + caseId);

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

	private static String getOtherPossibleHpoTerms(String ensg, LinkedHashSet<String> annotatedHpo, DoubleMatrixDataset<String, String> hpoPredictionMatrixSignificant, double minZscoreOtherCandidates, HpoOntology hpoOntology) {

		StringBuilder result = new StringBuilder();
		boolean notFirstResult = false;

		DoubleMatrix1D geneZscores = hpoPredictionMatrixSignificant.getRow(ensg);
		ArrayList<String> hpos = hpoPredictionMatrixSignificant.getColObjects();

		ArrayList<Term> patientHpoTerms = new ArrayList<>();
		for (String patientHpo : annotatedHpo) {
			patientHpoTerms.add(hpoOntology.nameToTerm(patientHpo));
		}

		for (int p = 0; p < hpos.size(); p++) {
			if (geneZscores.get(p) >= minZscoreOtherCandidates) {

				String potentialOtherHpo = hpos.get(p);

				if (potentialOtherHpo.equals(hpos)) {
					continue;
				}

				Term potentialOtherHpoTerm = hpoOntology.nameToTerm(potentialOtherHpo);

				//Only report if not parent or child of patient hpo
				boolean childOrParent = false;
				for (Term patientTerm : patientHpoTerms) {

					if (hpoOntology.isTermABelowTermB(patientTerm, potentialOtherHpoTerm)) {
						childOrParent = true;
						break;
					}

					if (hpoOntology.isTermABelowTermB(potentialOtherHpoTerm, patientTerm)) {
						childOrParent = true;
						break;
					}
				}

				if (!childOrParent) {

					//here with have a term that is not child or parent of already annotated HPO terms but has strong Z-score
					String potentialOtherHpoDescription = potentialOtherHpoTerm.getDescription();

					if (notFirstResult) {
						result.append('\t');
					}
					notFirstResult = true;

					result.append(Z_FORMAT.format(geneZscores.get(p)));
					result.append(';');
					result.append(potentialOtherHpo);
					result.append(';');
					result.append(potentialOtherHpoDescription);

				}

			}
		}

		return result.toString();

	}

}
