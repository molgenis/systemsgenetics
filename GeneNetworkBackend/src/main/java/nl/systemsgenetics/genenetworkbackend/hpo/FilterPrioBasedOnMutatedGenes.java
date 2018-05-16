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
import java.util.HashSet;

/**
 *
 * @author patri
 */
public class FilterPrioBasedOnMutatedGenes {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
//
//		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\samplesWithGeno.txt");
//		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\gavinRes\\");
//		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations");
//		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\rankingCandidateGenes");

//		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations2\\samplesWithGeno.txt");
//		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\gavinRes\\");
//		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations2");
//		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations2\\rankingCandidateGenes");

		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\samplesWithGeno.txt");
		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\gavinRes\\");
		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm");
		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\PrioritisationsDcm\\rankingCandidateGenes");


		resultFolder.mkdirs();
		
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader sampleFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(0).withCSVParser(parser).build();

		String[] nextLine;
		while ((nextLine = sampleFileReader.readNext()) != null) {

			String sample = nextLine[0];

			String genoSampleName = new File(nextLine[1]).getName();
			if (!genoSampleName.endsWith(".txt")) {
				genoSampleName += ".txt";
			}

			File genoFile = new File(genoFolder, genoSampleName);
			File prioFile = new File(prioFolder, sample + ".txt");
			File rankingFile = new File(resultFolder, sample + ".txt");

			System.out.println("------------------------------------------------------------------");
			System.out.println("Sample: " + sample);
			System.out.println("Geno: " + genoFile.getAbsolutePath());
			System.out.println("Prio: " + prioFile.getAbsolutePath());
			System.out.println("Ranking: " + rankingFile.getAbsolutePath());

			HashSet<String> genesWithMutation = getMutatedGenes(genoFile);

			final CSVReader prioFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(prioFile))).withSkipLines(0).withCSVParser(parser).build();

			CSVWriter writer = new CSVWriter(new FileWriter(rankingFile), '\t', '\0', '\0', "\n");

			String[] outputLine = prioFileReader.readNext();
			writer.writeNext(outputLine);

			while ((outputLine = prioFileReader.readNext()) != null) {

				if (genesWithMutation.contains(outputLine[1])) {
					writer.writeNext(outputLine);
				}

			}

			writer.close();
			prioFileReader.close();

		}

	}

	private static HashSet<String> getMutatedGenes(File genoFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(genoFile))).withSkipLines(1).withCSVParser(parser).build();

		HashSet<String> genes = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(nextLine[8]);

		}

		reader.close();

		return genes;

	}

}
