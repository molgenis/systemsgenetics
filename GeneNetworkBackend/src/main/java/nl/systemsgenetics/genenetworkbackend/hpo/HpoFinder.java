/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.io.OboParser;

/**
 *
 * @author patri
 */
public class HpoFinder {

	final Ontology hpoOntology;
	final Map<String, PredictionInfo> predictionInfo;
	final Term is_a;
	final Term is_obsolete = null;
	final Term trueValue = null;
	final Term replaced_by = null;

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, ParseException {

		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final File hpoPredictionInfoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_auc_fdr.txt");
		final File queryFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\originalHpo.txt");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\selectedHpo.txt");
		final double fdr = 0.05;

		Map<String, PredictionInfo> predictionInfo = loadPredictionInfo(hpoPredictionInfoFile);

		OboParser parser = new OboParser();

		BufferedReader oboFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(hpoOboFile)));
		Ontology hpoOntology = parser.parseOBO(oboFileReader, "HPO", "HPO");

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		int c = 0;
		String[] outputLine = new String[7];
		outputLine[c++] = "originalHPO";
		outputLine[c++] = "originalDescription";
		outputLine[c++] = "matchHPO";
		outputLine[c++] = "matchDescription";
		outputLine[c++] = "matchPvalue";
		outputLine[c++] = "matchauc";
		outputLine[c++] = "multiMatch";

		writer.writeNext(outputLine);

		BufferedReader queryReader = new BufferedReader(new FileReader(queryFile));
		String queryHpo;

		while ((queryHpo = queryReader.readLine()) != null) {

			if (hpoOntology.containsTerm(queryHpo)) {
				Term queryHpoTerm = hpoOntology.getTerm(queryHpo);

				List<Term> alternativeTerms = hpoFinder.getPredictableTerms(queryHpoTerm, fdr);

				for (Term alternativeTerm : alternativeTerms) {

					PredictionInfo info = predictionInfo.get(alternativeTerm.getName());

					c = 0;
					outputLine[c++] = queryHpo;
					outputLine[c++] = queryHpoTerm.getDescription();
					outputLine[c++] = alternativeTerm.getName();
					outputLine[c++] = alternativeTerm.getDescription();
					outputLine[c++] = String.valueOf(info.getpValue());
					outputLine[c++] = String.valueOf(info.getAuc());
					outputLine[c++] = alternativeTerms.size() > 1 ? "x" : "-";
					writer.writeNext(outputLine);

					//System.out.println(alternativeTerm.getName() + " P-value: " + info.getpValue() + " AUC: " + info.getAuc() + " " + alternativeTerm.getDescription());
				}

				if (alternativeTerms.isEmpty()) {
					c = 0;
					outputLine[c++] = queryHpo;
					outputLine[c++] = queryHpoTerm.getDescription();
					outputLine[c++] = "NA";
					outputLine[c++] = "NA";
					outputLine[c++] = "NA";
					outputLine[c++] = "NA";
					outputLine[c++] = "NA";
					writer.writeNext(outputLine);
				}

			} else {
				c = 0;
				outputLine[c++] = queryHpo;
				outputLine[c++] = "NA";
				outputLine[c++] = "NA";
				outputLine[c++] = "NA";
				outputLine[c++] = "NA";
				outputLine[c++] = "NA";
				outputLine[c++] = "NA";
				writer.writeNext(outputLine);
			}

		}

		writer.close();

	}

	public HpoFinder(Ontology hpoOntology, Map<String, PredictionInfo> predictionInfo) {
		this.hpoOntology = hpoOntology;
		this.predictionInfo = predictionInfo;

		is_a = hpoOntology.getTerm("is_a");
		//is_obsolete = hpoOntology.getTerm("is_obsolete");
		//trueValue = hpoOntology.getTerm("true");
		//replaced_by = hpoOntology.getTerm("replaced_by");

	}

	private static Map<String, PredictionInfo> loadPredictionInfo(File hpoPredictionInfoFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoPredictionInfoFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, PredictionInfo> predictionInfo = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			PredictionInfo info = new PredictionInfo(nextLine[0], Double.parseDouble(nextLine[1]), Double.parseDouble(nextLine[2]), Double.parseDouble(nextLine[3]));
			predictionInfo.put(info.getHpo(), info);

		}

		return Collections.unmodifiableMap(predictionInfo);

	}

	public List<Term> getPredictableTerms(Term queryTerm, double fdr) {

		List<Term> result;

		if (queryTerm.getAnnotation().containsProperty("is_obsolete")) {
			System.err.println("Obsolete term cannot be handeld. Please find alternative manually: " + queryTerm.getName());
			return Collections.emptyList();
		}

		PredictionInfo info = predictionInfo.get(queryTerm.getName());

		if (info == null || info.getFdr() > fdr) {

			if (info == null) {
				System.out.println("No predictions for: " + queryTerm.getName());
			} else {
				System.out.println("Bad predictions for: " + queryTerm.getName() + " P-value: " + info.getpValue() + " AUC: " + info.getAuc() + " FDR: " + info.getFdr());
			}

			result = new ArrayList<>();

			hpoOntology.getTriples(queryTerm, null, is_a).forEach((parentTriple) -> {
				result.addAll(getPredictableTerms(parentTriple.getObject(), fdr));
			});

		} else {
			result = new ArrayList<>(1);
			result.add(queryTerm);
		}

		return result;

	}

}
