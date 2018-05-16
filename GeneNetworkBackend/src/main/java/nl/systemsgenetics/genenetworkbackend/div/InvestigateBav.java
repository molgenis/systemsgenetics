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
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

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
		final File samplesFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\AtaxiaSamples.txt");
		final File mutationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\MutationMatrixAtaxia.txt");

		ArrayList<String> samples = loadLines(samplesFile);

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		TObjectIntHashMap<String> geneCount = new TObjectIntHashMap();
		TObjectDoubleHashMap<String> geneZscore = new TObjectDoubleHashMap<>();
		HashMap<String, String> geneOtherHpo = new HashMap<>();
		HashMap<String, HashSet<String>> geneMutatedSamples = new HashMap<>();

		for (String sample : samples) {

			File sampleFile = new File(rankFolder, sample + ".txt");
			
			if(!sampleFile.exists()){
				System.out.println("Skipping: " + sample );
				continue;
			}
			
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(1).withCSVParser(parser).build();

			for (int i = 0; i < 20; i++) {

				String[] line = reader.readNext();
				String gene = line[1];
				double zscore = Double.parseDouble(line[3]);
				String otherHpos = line[5];

				if(zscore > 3){
					System.out.println(sample + " " + gene + " " + zscore);
				}
				
				geneZscore.put(gene, zscore);
				geneCount.adjustOrPutValue(gene, 1, 1);
				geneOtherHpo.put(gene, otherHpos);

				geneMutatedSamples.putIfAbsent(gene, new HashSet<>());
				geneMutatedSamples.get(gene).add(sample);

			}

			reader.close();

		}

		ArrayList<String> recurrentGenes = new ArrayList<>();
		for (String gene : geneCount.keySet()) {
			int count = geneCount.get(gene);

			if (count > 1) {

				if (geneZscore.get(gene) > 3) {
					recurrentGenes.add(gene);
					System.out.println(gene + "\t" + count + "\t" + geneZscore.get(gene) + "\t" + geneOtherHpo.get(gene) + "\t" + String.join(";", geneMutatedSamples.get(gene)));
				}

				
			}
		}

		DoubleMatrixDataset<String, String> mutationMatrix = new DoubleMatrixDataset<>(recurrentGenes, samples);

		for (String sample : samples) {

			File sampleFile = new File(rankFolder, sample + ".txt");
			
			if(!sampleFile.exists()){
				System.out.println("Skipping: " + sample );
				continue;
			}
			
			
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(1).withCSVParser(parser).build();

			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				String gene = nextLine[1];

				if (mutationMatrix.containsRow(gene)) {

					mutationMatrix.setElement(gene, sample, 1);

				}

			}

		}

		mutationMatrix.save(mutationMatrixFile);

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
