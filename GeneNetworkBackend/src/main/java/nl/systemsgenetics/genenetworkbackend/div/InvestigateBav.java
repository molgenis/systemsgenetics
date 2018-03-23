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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author patri
 */
public class InvestigateBav {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		final File rankFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\rankingCandidateGenes");
		final File samplesFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\BavSamples.txt");

		ArrayList<String> samples = loadLines(samplesFile);
		
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		
		
		for(String sample : samples){
			
			File sampleFile = new File(rankFolder, sample + ".txt");
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(1).withCSVParser(parser).build();
			
			
			
			
		}
		
		
	}
	
	
	public static ArrayList<String> loadLines(File file) throws IOException {

		ArrayList<String> lines = new ArrayList<>();

		BufferedReader reader = new BufferedReader(new FileReader(file));

		String line;
		while ((line = reader.readLine()) != null) {
			lines.add(line);
		}

		return lines;

	}

}
