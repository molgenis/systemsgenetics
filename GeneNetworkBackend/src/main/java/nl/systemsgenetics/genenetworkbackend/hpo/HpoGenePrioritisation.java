/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import umcg.genetica.io.hpo.HpoOntology;
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
import org.biojava.nbio.ontology.Term;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class HpoGenePrioritisation {

	private static final NumberFormat Z_FORMAT = new DecimalFormat("#0.0##");

	public static void main(String[] args) throws IOException, ParseException, Exception {

		//final File hpoPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions.txt.gz");
		//spliked
		final File hpoPredictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_spiked.txt");
		final File ensgSymbolMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_bonSigTerms.txt");
		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final double minZscoreOtherCandidates = 5;

//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\selectedHpo2.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\HPO_termen_Linda_processed.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\Prioritisations");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\hpos_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\hpoSolvedCases_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations3Spiked");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\hpo_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\Prioritisations");

//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\extraUnsolved\\hpo_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\extraUnsolved\\Prioritisations");

		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Vibe\\hpo_selected.txt");
		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Vibe\\Prioritisations");


//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioMieke\\hpo_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioMieke\\Prioritisations\\");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioEdgar\\hpo_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioEdgar\\Prioritisations\\");
//		final File caseHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\orginalcasehpo_EditedSipko3_JanJongbloed_selected.txt");
//		final File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations2");
		outputFolder.mkdirs();

		Map<String, String> ensgSymbolMapping = loadEnsgToHgnc(ensgSymbolMappingFile);

		HashMap<String, LinkedHashSet<String>> caseHpo = loadCaseHpo(caseHpoFile);

		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		HpoOntology hpoOntology = new HpoOntology(hpoOboFile);

		DoubleMatrixDataset<String, String> hpoPredictionMatrix = DoubleMatrixDataset.loadDoubleData(hpoPredictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> hpoPredictionMatrixSignificant = hpoPredictionMatrix.viewColSelection(significantTerms);
		ArrayList<String> genes = hpoPredictionMatrixSignificant.getRowObjects();

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
				if (!hpoPredictionMatrixSignificant.containsCol(term)) {
					System.err.println("Missing HPO: " + term);
				}
			}

			DoubleMatrixDataset<String, String> predictionCaseTerms = hpoPredictionMatrixSignificant.viewColSelection(hpo);
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

			String[] outputLine = new String[5 + hpo.size()];
			int c = 0;
			outputLine[c++] = "Ensg";
			outputLine[c++] = "Hgnc";
			outputLine[c++] = "Rank";
			outputLine[c++] = "Zscore";
			for (String term : hpo) {
				outputLine[c++] = term;
			}
			outputLine[c++] = "PossibleAdditionalHpos";
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

				outputLine[c++] = getOtherPossibleHpoTerms(geneResult.getEnsg(), hpo, hpoPredictionMatrixSignificant, minZscoreOtherCandidates, hpoOntology);

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

	public static LinkedHashSet<String> loadSignificantTerms(File significantTermsFile) throws IOException {

		LinkedHashSet<String> significantTerms = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(significantTermsFile));

		String line;
		while ((line = reader.readLine()) != null) {
			significantTerms.add(line);
		}

		return significantTerms;

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
