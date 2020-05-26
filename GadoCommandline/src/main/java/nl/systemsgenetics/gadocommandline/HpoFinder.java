/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

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

	public static Ontology loadHpoOntology(final File hpoOboFile) throws FileNotFoundException, IOException, ParseException {
		OboParser parser = new OboParser();
		BufferedReader oboFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(hpoOboFile)));
		Ontology hpoOntology = parser.parseOBO(oboFileReader, "HPO", "HPO");
		return hpoOntology;
	}

	public HpoFinder(Ontology hpoOntology, Map<String, PredictionInfo> predictionInfo) {
		this.hpoOntology = hpoOntology;
		this.predictionInfo = predictionInfo;

		is_a = hpoOntology.getTerm("is_a");
		//is_obsolete = hpoOntology.getTerm("is_obsolete");
		//trueValue = hpoOntology.getTerm("true");
		//replaced_by = hpoOntology.getTerm("replaced_by");

	}

	public static Map<String, PredictionInfo> loadPredictionInfo(File hpoPredictionInfoFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoPredictionInfoFile))).withSkipLines(1).withCSVParser(parser).build();

		HashMap<String, PredictionInfo> predictionInfo = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			PredictionInfo info = new PredictionInfo(nextLine[0], Double.parseDouble(nextLine[2]), Double.parseDouble(nextLine[3]), Double.parseDouble(nextLine[4]));
			predictionInfo.put(info.getHpo(), info);

		}

		return Collections.unmodifiableMap(predictionInfo);

	}

	public List<Term> getPredictableTerms(String queryName, double correctedPCutoff) {
		return getPredictableTerms(hpoOntology.getTerm(queryName), correctedPCutoff);
	}

	public List<Term> getPredictableTerms(Term queryTerm, double correctedPCutoff) {

		List<Term> result;

		if (queryTerm.getAnnotation().containsProperty("is_obsolete")) {
			throw new RuntimeException("Obsolete term cannot be handeld. Please find alternative manually: " + queryTerm.getName());
		}

		PredictionInfo info = predictionInfo.get(queryTerm.getName());

		if (info == null || info.getCorrectedP() > correctedPCutoff) {

			if (info == null) {
				//System.out.println("No predictions for: " + queryTerm.getName());
			} else {
				//System.out.println("Bad predictions for: " + queryTerm.getName() + " P-value: " + info.getpValue() + " AUC: " + info.getAuc() + " FDR: " + info.getFdr());
			}

			result = new ArrayList<>();

			hpoOntology.getTriples(queryTerm, null, is_a).forEach((parentTriple) -> {
				result.addAll(getPredictableTerms(parentTriple.getObject(), correctedPCutoff));
			});

		} else {
			result = new ArrayList<>(1);
			result.add(queryTerm);
		}

		return result;

	}
	
	public List<String> getTermsToNames(List<Term> terms){
		ArrayList<String> termNames = new ArrayList<>();
		
		terms.forEach((term) -> {
			termNames.add(term.getName());
		});
		
		return termNames;
		
	}

}
