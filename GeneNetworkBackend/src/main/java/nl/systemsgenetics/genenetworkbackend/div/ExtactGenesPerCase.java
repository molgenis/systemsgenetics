/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author patri
 */
public class ExtactGenesPerCase {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		File[] folders = new File[]{
			new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\rankingCandidateGenes"),
			new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsCardioEdgar\\"),
			new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\extraUnsolved\\rankingCandidateGenes"),
			new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations3\\rankingCandidateGenes")};

		File sampleMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\sampleMappings.txt");
		File outputFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\sampleGenesForManuscript\\");

		outputFolder.mkdirs();

		HashMap<String, String> dnasToPsuedos = new HashMap<>();

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleMappingFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			dnasToPsuedos.put(nextLine[0], nextLine[1]);
		}
		reader.close();

		int samplesProcessed = 0;

		for (Map.Entry<String, String> sample : dnasToPsuedos.entrySet()) {

			String dnaNumber = sample.getKey();
			String pseudoId = sample.getValue();

			File sampleFile = null;

			for (File folder : folders) {

				sampleFile = new File(folder, dnaNumber + ".txt");

				if (sampleFile.exists()) {
					break;
				}

			}

			if (!sampleFile.exists()) {
				System.err.println("Sample not found: " + dnaNumber);
				continue;
			}

			File outputFile = new File(outputFolder, pseudoId + ".txt");

			final CSVReader reader2 = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(0).withCSVParser(parser).build();

			CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

			String[] outputLine = new String[2];
			int c;

			while ((nextLine = reader2.readNext()) != null) {
				c = 0;
				outputLine[c++] = nextLine[0];
				outputLine[c++] = nextLine[1];
				writer.writeNext(outputLine);
			}

			reader2.close();
			writer.close();

			++samplesProcessed;

		}

		System.out.println(samplesProcessed + " out of " + dnasToPsuedos.size() + " done");

	}

}
