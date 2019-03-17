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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;

/**
 *
 * @author patri
 */
public class processCaseHpo {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException, ParseException {

		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final File hpoPredictionInfoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_auc_bonferroni.txt");
		final File updatedIdFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\updatedHpoId.txt");
//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\orginalCaseHpo.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\selectedHpo.txt");
//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\unsolvedDcmHpo.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\unsolvedDcmHpoSelected.txt");
//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\hpos.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\hpos_selected.txt");

//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\hpoSolvedCases.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\hpoSolvedCases_selected.txt");

//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioMieke\\hpos.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioMieke\\hpo_selected.txt");

//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\hpo5gpm.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\hpo_selected.txt");

		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\extraUnsolved\\hpo.txt");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\extraUnsolved\\hpo_selected.txt");

		

//		final File caseHpo = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\HPO_5GPM.txt");
//		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\HPO_5GPM_processed.txt");
		final double correctedPCutoff = 0.05;

		Map<String, PredictionInfo> predictionInfo = HpoFinder.loadPredictionInfo(hpoPredictionInfoFile);

		Ontology hpoOntology = HpoFinder.loadHpoOntology(hpoOboFile);

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		Map<String, String> updatedHpoId = loadUpdatedIds(updatedIdFile);

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(caseHpo))).withCSVParser(parser).build();

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[6];
		int c = 0;
		outputLine[c++] = "Sample";
		outputLine[c++] = "SelectedHpo";
		outputLine[c++] = "SelectedHpoDescription";
		outputLine[c++] = "OriginalHpo";
		outputLine[c++] = "OriginalHpoDescription";
		outputLine[c++] = "ExcludeFromPrioritisation";
		writer.writeNext(outputLine);

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String sampleId = nextLine[0];
			HashSet<String> sampleHpo = new HashSet<>();

			for (int i = 1; i < nextLine.length; i++) {

				String hpo = nextLine[i];
				
				if(hpo.length() == 0){
					continue;
				}

				if (updatedHpoId.containsKey(hpo)) {
					hpo = updatedHpoId.get(hpo);
				}

				if (sampleHpo.add(hpo)) {

					Term hpoTerm = hpoOntology.getTerm(hpo);
					PredictionInfo info = predictionInfo.get(hpo);

					if (info == null || info.getCorrectedP() > correctedPCutoff) {
						//in case of no prediction or bad prediction

						List<Term> alternativeTerms = hpoFinder.getPredictableTerms(hpoTerm, correctedPCutoff);

						if(alternativeTerms.isEmpty()){
							System.out.println("Warning no alternative found for: " + hpo);
						}
						
						for (Term alternativeTerm : alternativeTerms) {
							c = 0;
							outputLine[c++] = sampleId;
							outputLine[c++] = alternativeTerm.getName();
							outputLine[c++] = alternativeTerm.getDescription();
							outputLine[c++] = hpo;
							outputLine[c++] = hpoTerm.getDescription();
							outputLine[c++] = "";
							writer.writeNext(outputLine);
						}

					} else {

						c = 0;
						outputLine[c++] = sampleId;
						outputLine[c++] = hpo;
						outputLine[c++] = hpoTerm.getDescription();
						outputLine[c++] = "";
						outputLine[c++] = "";
						outputLine[c++] = "";
						writer.writeNext(outputLine);

					}

				}

			}

		}
		
		writer.close();

	}

	private static Map<String, String> loadUpdatedIds(File updatedIdFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(updatedIdFile))).withCSVParser(parser).build();

		HashMap<String, String> updates = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			updates.put(nextLine[0], nextLine[1]);

		}

		return Collections.unmodifiableMap(updates);

	}

}
